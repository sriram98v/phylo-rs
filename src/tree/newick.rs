//! Iterative Newick-format parser.
//!
//! Parses the input in a single left-to-right pass, building directly into the
//! tree arena. Nesting is tracked with an explicit `Vec` of ancestor node ids
//! rather than by recursion, so parsing a tree nested `d` levels deep uses
//! `O(d)` heap and no call-stack depth — a pectinate ("caterpillar") tree of
//! `n` tips parses without overflowing the stack. There is no intermediate AST
//! or token buffer, so peak memory is just the tree plus the ancestor stack.
//!
//! # Accepted grammar
//!
//! The parser aims for the common ground of Newick as emitted by popular
//! phylogenetics software (RAxML, IQ-TREE, RevBayes, BEAST, ...), restricted to
//! *trees* (network / extended-Newick `#H` hybrid nodes are out of scope):
//!
//! * **Topology** via `(`, `)`, `,`, terminated by `;` (a missing terminator at
//!   end of input is tolerated; anything after the first `;` is ignored, so a
//!   multi-tree file yields its first tree).
//! * **Quoted labels** `'...'`, where `''` is an escaped single quote and all
//!   other characters — including whitespace, `:`, `;`, `,` — are preserved.
//! * **Unquoted labels**: any run of characters that are not whitespace and not
//!   one of `( ) [ ] , : ; '`. Underscores are kept **literal** (not converted
//!   to spaces), matching modern tools and keeping identifiers such as
//!   `GB_GCA_015163815.1` intact. Support values such as IQ-TREE's `95.3/98`
//!   are ordinary unquoted labels.
//! * **Branch lengths** `:<number>` in decimal or scientific notation. A token
//!   that does not parse as a weight is stored as no weight rather than being
//!   rejected.
//! * **Comments** `[...]` anywhere whitespace is allowed (NHX `[&&NHX:...]`,
//!   BEAST `[&...]`, RevBayes `[&index=..]`, rooting markers `[&R]`/`[&U]`).
//!   Their contents are ignored.
//!
//! Malformed input yields a [`NewickError`] carrying a byte offset; the parser
//! never panics.

use std::sync::Arc;

use crate::error::NewickError;
use crate::node::Node;
use crate::prelude::*;

/// Returns true if `c` terminates an unquoted label or a branch-length token.
/// This is exactly the set of characters the Newick standard forbids in an
/// unquoted label.
fn is_delimiter(c: char) -> bool {
    c.is_whitespace() || matches!(c, '(' | ')' | '[' | ']' | ',' | ':' | ';' | '\'')
}

/// A cursor over the source string tracking a byte offset.
struct Scanner<'a> {
    src: &'a str,
    pos: usize,
}

impl<'a> Scanner<'a> {
    fn new(src: &'a str) -> Self {
        Scanner { src, pos: 0 }
    }

    fn rest(&self) -> &'a str {
        &self.src[self.pos..]
    }

    fn peek(&self) -> Option<char> {
        self.rest().chars().next()
    }

    /// Advances past the current character and returns it.
    fn bump(&mut self) -> Option<char> {
        let c = self.peek()?;
        self.pos += c.len_utf8();
        Some(c)
    }

    /// Skips whitespace and `[...]` comments. Each comment's raw text (brackets
    /// included) is appended to `sink` so the caller can preserve it as a node
    /// annotation. Returns an error on an unterminated comment.
    fn skip_trivia(&mut self, sink: &mut String) -> Result<(), NewickError> {
        loop {
            let trimmed = self.rest().trim_start();
            self.pos = self.src.len() - trimmed.len();
            if self.rest().starts_with('[') {
                let start = self.pos;
                match self.rest().find(']') {
                    Some(rel) => {
                        let end = self.pos + rel + ']'.len_utf8();
                        sink.push_str(&self.src[start..end]);
                        self.pos = end;
                    }
                    None => return Err(NewickError::UnterminatedComment { idx: start }),
                }
            } else {
                return Ok(());
            }
        }
    }

    /// Reads a run of non-delimiter characters, returned as a borrowed slice.
    fn read_token(&mut self) -> &'a str {
        let start = self.pos;
        while let Some(c) = self.peek() {
            if is_delimiter(c) {
                break;
            }
            self.pos += c.len_utf8();
        }
        &self.src[start..self.pos]
    }

    /// Reads a `'...'` quoted label, unescaping `''` to a single `'`. Assumes
    /// the current character is the opening quote.
    fn read_quoted(&mut self) -> Result<String, NewickError> {
        let start = self.pos;
        self.bump(); // opening quote
        let mut out = String::new();
        loop {
            match self.bump() {
                None => return Err(NewickError::UnterminatedQuote { idx: start }),
                Some('\'') => {
                    if self.peek() == Some('\'') {
                        self.bump();
                        out.push('\'');
                    } else {
                        return Ok(out);
                    }
                }
                Some(c) => out.push(c),
            }
        }
    }
}

/// Commits the raw comment text gathered for node `id` via `handler`, then
/// clears the buffer. The handler decides whether and in what form the
/// annotation is stored; if a node is committed more than once (comments both
/// before and after its subtree) the stored values are concatenated.
fn commit_annotation<T, W, Z, H>(
    tree: &mut SimpleRootedTree<T, W, Z>,
    id: TreeNodeID<SimpleRootedTree<T, W, Z>>,
    raw: &mut String,
    handler: &H,
) where
    T: NodeTaxa,
    W: EdgeWeight,
    Z: NodeWeight,
    H: AnnotationHandler,
{
    if raw.is_empty() {
        return;
    }
    if let Some(stored) = handler.handle(raw) {
        if let Some(node) = tree.get_node_mut(id) {
            let combined = match node.get_annotation() {
                Some(existing) => Arc::from(format!("{existing}{stored}")),
                None => stored,
            };
            node.set_annotation(Some(combined));
        }
    }
    raw.clear();
}

/// Parses a Newick string into a [`SimpleRootedTree`], retaining node `[...]`
/// annotations as decided by `handler`.
///
/// Only the first `;`-terminated tree is read. See the module documentation for
/// the accepted grammar.
pub(crate) fn parse_newick<T, W, Z, H>(
    src: &str,
    handler: &H,
) -> Result<SimpleRootedTree<T, W, Z>, NewickError>
where
    T: NodeTaxa,
    W: EdgeWeight,
    Z: NodeWeight,
    H: AnnotationHandler,
{
    type Id<T, W, Z> = TreeNodeID<SimpleRootedTree<T, W, Z>>;

    let mut sc = Scanner::new(src);
    let mut tree = SimpleRootedTree::new(0);
    // `current` is the node that a label or branch length applies to; `stack`
    // holds its ancestors, deepest last.
    let mut current: Id<T, W, Z> = tree.get_root_id();
    let mut stack: Vec<Id<T, W, Z>> = Vec::new();
    let mut saw_content = false;
    // Raw `[...]` comment text gathered while `current` is unchanged. Committed
    // (and cleared) whenever a topology token moves to a different node, so the
    // handler sees all of a node's comments together.
    let mut raw = String::new();

    loop {
        sc.skip_trivia(&mut raw)?;
        let c = match sc.peek() {
            // End of input without a terminator is tolerated.
            None => break,
            // First `;` ends the first tree; the rest of the input is ignored.
            Some(';') => break,
            Some(c) => c,
        };

        match c {
            '(' => {
                sc.bump();
                // Descend: `current` becomes a parent, open its first child.
                commit_annotation(&mut tree, current, &mut raw, handler);
                stack.push(current);
                let child = tree.next_id();
                tree.set_node(Node::new(child));
                tree.set_child(current, child);
                current = child;
                saw_content = true;
            }
            ',' => {
                sc.bump();
                // Start a sibling under the same parent.
                commit_annotation(&mut tree, current, &mut raw, handler);
                let parent = *stack
                    .last()
                    .ok_or(NewickError::UnbalancedParens { idx: sc.pos })?;
                let sibling = tree.next_id();
                tree.set_node(Node::new(sibling));
                tree.set_child(parent, sibling);
                current = sibling;
            }
            ')' => {
                sc.bump();
                // Ascend: subsequent label/length annotate the parent node.
                commit_annotation(&mut tree, current, &mut raw, handler);
                current = stack
                    .pop()
                    .ok_or(NewickError::UnbalancedParens { idx: sc.pos })?;
            }
            ':' => {
                sc.bump();
                sc.skip_trivia(&mut raw)?;
                let weight = sc
                    .read_token()
                    .parse::<TreeNodeWeight<SimpleRootedTree<T, W, Z>>>();
                if let Some(node) = tree.get_node_mut(current) {
                    node.set_weight(weight.ok());
                }
            }
            '\'' => {
                let label = sc.read_quoted()?;
                tree.set_node_taxa(current, T::from_str(&label).ok());
                saw_content = true;
            }
            _ => {
                let label = sc.read_token();
                if label.is_empty() {
                    // A delimiter with no meaning here (e.g. a stray `]`).
                    return Err(NewickError::InvalidCharacter { idx: sc.pos });
                }
                tree.set_node_taxa(current, T::from_str(label).ok());
                saw_content = true;
            }
        }
    }
    // Commit any comments trailing the final node.
    commit_annotation(&mut tree, current, &mut raw, handler);

    if !stack.is_empty() {
        return Err(NewickError::UnbalancedParens { idx: sc.pos });
    }
    if !saw_content {
        return Err(NewickError::Empty);
    }
    Ok(tree)
}
