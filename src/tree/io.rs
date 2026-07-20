use itertools::Itertools;
use std::borrow::Cow;
use std::ffi::OsStr;
use std::fmt::Display;
use std::path::Path;
use std::sync::Arc;
use std::{fs, io};

use crate::prelude::*;

/// Enum to track block of Nexus file. This enum can be extended in the future to include new blocks for different use cases.
pub enum NexusBlock {
    /// Tree block
    TREE,
    /// Miscellaneous block to be ignored
    NONE,
}

/// Decides whether and how a node's raw Newick `[...]` annotation is retained
/// while parsing.
///
/// The parser gathers every `[...]` comment attached to a node (brackets
/// included, concatenated in source order) and passes the raw text to
/// [`AnnotationHandler::handle`]. Returning `Some` stores that value on the
/// node; returning `None` drops it. Any closure `Fn(&str) -> Option<Arc<str>>`
/// is a handler, so callers can transform annotations however they like — strip
/// brackets, keep only certain keys, normalise NHX, etc.
pub trait AnnotationHandler {
    /// Maps the raw annotation text for one node to the value to store, or
    /// `None` to discard it.
    fn handle(&self, raw: &str) -> Option<Arc<str>>;
}

/// Retains each node's annotation exactly as written (lossless). This is the
/// default used by [`Newick::from_newick`].
#[derive(Debug, Clone, Copy, Default)]
pub struct KeepRawAnnotations;

impl AnnotationHandler for KeepRawAnnotations {
    fn handle(&self, raw: &str) -> Option<Arc<str>> {
        Some(Arc::from(raw))
    }
}

/// Drops all annotations, leaving every node's annotation empty.
#[derive(Debug, Clone, Copy, Default)]
pub struct DiscardAnnotations;

impl AnnotationHandler for DiscardAnnotations {
    fn handle(&self, _raw: &str) -> Option<Arc<str>> {
        None
    }
}

impl<F> AnnotationHandler for F
where
    F: Fn(&str) -> Option<Arc<str>>,
{
    fn handle(&self, raw: &str) -> Option<Arc<str>> {
        self(raw)
    }
}

/// Decides whether and how a node's stored annotation is written back out by
/// [`Newick::to_newick_with`] / [`Newick::subtree_to_newick_with`].
///
/// For each annotated node the stored annotation text is passed to
/// [`AnnotationWriter::render`]. Returning `Some` emits that text (right after
/// the node label); returning `None` omits the annotation. Any closure
/// `Fn(&str) -> Option<String>` is a writer, so callers can filter or rewrite
/// what ends up in the output.
pub trait AnnotationWriter {
    /// Maps a node's stored annotation to the text to emit, or `None` to omit
    /// it. Borrowing the input (via `Cow::Borrowed`) avoids allocating when the
    /// annotation is emitted unchanged.
    fn render<'a>(&self, annotation: &'a str) -> Option<Cow<'a, str>>;
}

impl AnnotationWriter for KeepRawAnnotations {
    fn render<'a>(&self, annotation: &'a str) -> Option<Cow<'a, str>> {
        Some(Cow::Borrowed(annotation))
    }
}

impl AnnotationWriter for DiscardAnnotations {
    fn render<'a>(&self, _annotation: &'a str) -> Option<Cow<'a, str>> {
        None
    }
}

impl<F> AnnotationWriter for F
where
    F: Fn(&str) -> Option<String>,
{
    fn render<'a>(&self, annotation: &'a str) -> Option<Cow<'a, str>> {
        self(annotation).map(Cow::Owned)
    }
}

/// A trait descibing Newick encoding of a tree.
pub trait Newick: RootedTree {
    /// Creates a new tree from a Newick string, retaining node `[...]`
    /// annotations according to `annotations`.
    ///
    /// Pass [`KeepRawAnnotations`] to keep them verbatim, [`DiscardAnnotations`]
    /// to drop them, or any `Fn(&str) -> Option<Arc<str>>` closure to transform
    /// them.
    fn from_newick_with<H: AnnotationHandler>(
        newick_str: &[u8],
        annotations: H,
    ) -> std::io::Result<Self>;

    /// Creates a new tree using a Newick string, keeping node annotations
    /// verbatim. Equivalent to [`Newick::from_newick_with`] with
    /// [`KeepRawAnnotations`].
    fn from_newick(newick_str: &[u8]) -> std::io::Result<Self> {
        Self::from_newick_with(newick_str, KeepRawAnnotations)
    }

    /// Encodes a subtree starting from a node as a Newick string, emitting node
    /// annotations as decided by `annotations`.
    ///
    /// Pass [`KeepRawAnnotations`] to write them verbatim, [`DiscardAnnotations`]
    /// to omit them, or any `Fn(&str) -> Option<String>` closure to filter or
    /// rewrite them.
    fn subtree_to_newick_with<H: AnnotationWriter>(
        &self,
        node_id: TreeNodeID<Self>,
        annotations: H,
    ) -> impl Display;

    /// Encodes a subtree starting from a node as a Newick string, writing node
    /// annotations verbatim.
    fn subtree_to_newick(&self, node_id: TreeNodeID<Self>) -> impl Display {
        self.subtree_to_newick_with(node_id, KeepRawAnnotations)
    }

    /// Encodes a tree as a Newick string, emitting node annotations as decided
    /// by `annotations`.
    fn to_newick_with<H: AnnotationWriter>(&self, annotations: H) -> impl Display {
        format!(
            "{};",
            self.subtree_to_newick_with(self.get_root_id(), annotations)
        )
    }

    /// Encodes a tree as a Newick string, writing node annotations verbatim.
    fn to_newick(&self) -> impl Display {
        format!("{};", self.subtree_to_newick(self.get_root_id()))
    }

    /// Writes Newick String to file
    fn to_file(&self, p: &Path) -> io::Result<()> {
        assert!(p.extension() == Some(OsStr::new("nwk")));
        fs::write(p, self.to_newick().to_string().as_bytes())
    }

    /// Reads Newick String to file
    /// Note: this attempts to read only the first tree in the file
    fn from_file(p: &Path) -> io::Result<Self> {
        assert!(p.extension() == Some(OsStr::new("nwk")));
        let nwk_string = fs::read_to_string(p)?
            .as_bytes()
            .iter()
            .copied()
            .take_while(|x| *x != b';')
            .collect_vec();

        Self::from_newick(nwk_string.as_slice())
    }
}

/// A trait for reading and writing Nexus files
pub trait Nexus: Newick {
    /// Creates tree from Nexus string
    /// Note: this attempts to read only the first tree in the file
    fn from_nexus(p: String) -> std::io::Result<Self> {
        let file_lines = p.lines().collect_vec();
        // Reject anything that doesn't open with the `#NEXUS` marker -- this
        // also covers empty input, where there is no first line to inspect.
        if file_lines.first() != Some(&"#NEXUS") {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                NexusError::InvalidHeader,
            ));
        }
        let mut tree_block = String::new();
        let mut curr_block = NexusBlock::NONE;
        for line in file_lines {
            let line_words = line
                .split_ascii_whitespace()
                .map(|x| x.to_ascii_lowercase())
                .collect_vec();
            match line_words.first().map(String::as_str) {
                None => continue,
                Some("begin") => {
                    // A bare `begin` with no block name is not a trees block.
                    curr_block = match line_words.get(1).map(String::as_str) {
                        Some("trees;") => NexusBlock::TREE,
                        _ => NexusBlock::NONE,
                    };
                }
                Some("end;") => curr_block = NexusBlock::NONE,
                Some(_) => {
                    if matches!(curr_block, NexusBlock::TREE) {
                        tree_block.push_str(line);
                    }
                }
            }
        }
        // Take the newick after the first `=` of the first `;`-terminated
        // definition. Malformed or absent tree blocks yield an error rather
        // than panicking on an out-of-bounds index.
        let first_tree = tree_block
            .split(';')
            .next()
            .and_then(|def| def.split_once('='))
            .map(|(_, newick)| newick.split_whitespace().collect::<String>())
            .filter(|newick| !newick.is_empty())
            .ok_or_else(|| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    NexusError::MissingTreeBlock,
                )
            })?;
        Self::from_newick(format!("{first_tree};").as_bytes())
    }

    /// Creates tree from Nexus file
    fn from_nexus_file(p: &Path) -> std::io::Result<Self> {
        assert!(p.extension() == Some(OsStr::new("nex")));
        let file_data = fs::read_to_string(p)?;
        Self::from_nexus(file_data)
    }

    /// Writes Newick String to file
    fn to_nexus(&self) -> io::Result<String> {
        Ok(format!(
            "#NEXUS\n\nBEGIN TREES;\n\tTree tree={}\nEND;",
            self.to_newick()
        ))
    }

    /// Writes Newick String to file
    fn to_nexus_file(&self, p: &Path) -> io::Result<()> {
        assert!(p.extension() == Some(OsStr::new("nex")));
        fs::write(p, self.to_nexus()?.as_bytes())
    }
}
