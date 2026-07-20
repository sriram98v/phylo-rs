#[cfg(feature = "non_crypto_hash")]
use fxhash::{FxHashMap as HashMap, FxHashSet as HashSet};
#[cfg(feature = "parallel")]
use rayon::prelude::*;
#[cfg(not(feature = "non_crypto_hash"))]
use std::collections::{HashMap, HashSet};

use itertools::Itertools;
use phylo::node::PhyloNode;
use phylo::prelude::*;
use phylo::tree::PhyloTree;
use std::fs::{read_to_string, File};
use std::io::Write;

#[test]
fn pd() {
    let paths: HashMap<_, _> = std::fs::read_dir("examples/phylogenetic-diversity/trees")
        .unwrap()
        .map(|x| {
            (
                x.as_ref().unwrap().file_name().into_string().unwrap(),
                std::fs::read_dir(x.unwrap().path())
                    .unwrap()
                    .map(|f| {
                        (
                            f.as_ref()
                                .unwrap()
                                .file_name()
                                .into_string()
                                .unwrap()
                                .split("-")
                                .map(|x| x.to_string())
                                .collect_vec()[0]
                                .clone(),
                            PhyloTree::from_newick(
                                read_to_string(f.unwrap().path()).unwrap().as_bytes(),
                            )
                            .unwrap(),
                        )
                    })
                    .collect::<HashMap<_, _>>(),
            )
        })
        .collect();

    let mut output_file = File::create("examples/phylogenetic-diversity/pds.out").unwrap();

    for (clade, trees) in paths.iter() {
        println!("{}", clade);
        let mut pds = vec![];
        for year in 2015..2023 {
            let tree = trees.get(&year.to_string());
            match tree {
                Some(t) => {
                    println!(
                        "{}: {}",
                        year,
                        t.get_nodes()
                            .map(|n| n.get_weight().unwrap_or(0.0))
                            .sum::<f32>()
                    );
                    pds.push(
                        t.get_nodes()
                            .map(|n| n.get_weight().unwrap_or(0.0))
                            .sum::<f32>(),
                    );
                }
                _ => {
                    println!("{}: {}", year, 0.0);
                    pds.push(0.0);
                }
            };
        }
        println!("{:?}", pds);
        let out = format!(
            "{}: {}\n",
            clade,
            pds.iter().map(|x| x.to_string()).join(",")
        );
        output_file.write_all(out.as_bytes()).unwrap()
    }
}

#[test]
fn distance_matrix() {
    let input_str = String::from("((A:0.1,B:0.2):0.3,C:0.6);");
    let tree = PhyloTree::from_newick(input_str.as_bytes()).unwrap();
    let matrix = tree.matrix();
    dbg!(&matrix);
}

#[test]
fn build_small_tree() {
    let mut tree = PhyloTree::new(1);
    dbg!(&tree);
    let new_node = PhyloNode::new(2);
    tree.add_child(tree.get_root_id(), new_node);
    let new_node = PhyloNode::new(3);
    tree.add_child(tree.get_root_id(), new_node);
    let new_node = PhyloNode::new(4);
    tree.add_child(2, new_node);
    let new_node = PhyloNode::new(5);
    tree.add_child(2, new_node);
    dbg!(&tree, tree.get_node(1).unwrap().get_children());
    dbg!(RootedTree::get_node_depth(&tree, 2));
    dbg!(&tree.to_newick().to_string());
    tree.clear();
    dbg!(&tree);
}

#[test]
fn tree_iter() {
    let mut tree = PhyloTree::new(1);
    let new_node = PhyloNode::new(2);
    tree.add_child(tree.get_root_id(), new_node);
    let new_node = PhyloNode::new(5);
    tree.add_child(tree.get_root_id(), new_node);
    let new_node = PhyloNode::new(3);
    tree.add_child(2, new_node);
    let new_node = PhyloNode::new(4);
    tree.add_child(2, new_node);
    let new_node = PhyloNode::new(6);
    tree.add_child(5, new_node);
    let new_node = PhyloNode::new(7);
    tree.add_child(5, new_node);
    dbg!(&tree.get_node(1).unwrap().get_children());
    dbg!(&tree.dfs(tree.get_root_id()).collect_vec());
    dbg!(&tree.bfs_ids(tree.get_root_id()).collect_vec());
    dbg!(&tree.postord_ids(tree.get_root_id()).collect_vec());
    dbg!(&tree.euler_walk_ids(tree.get_root_id()).collect_vec());
    dbg!(&tree.dfs(tree.get_root_id()).collect_vec());
    dbg!(&tree.node_to_root(5).collect_vec());
    dbg!(&tree.root_to_node(5).collect_vec());
}
#[test]
fn read_small_tree() {
    let input_str = String::from("((A,B),C);");
    let tree = PhyloTree::from_newick(input_str.as_bytes()).unwrap();
    dbg!(&tree.euler_walk_ids(tree.get_root_id()).collect_vec());
    let input_str = String::from("((A:0.1,B:0.2),C:0.6);");
    let tree = PhyloTree::from_newick(input_str.as_bytes()).unwrap();
    dbg!(&tree.euler_walk_ids(tree.get_root_id()).collect_vec());
    dbg!(format!("{}", &tree.to_newick()));
    assert_eq!(
        &tree.get_taxa_space().collect::<HashSet<&String>>(),
        &vec![&"A".to_string(), &"B".to_string(), &"C".to_string()]
            .into_iter()
            .collect()
    );
    let input_str = String::from("((A:1e-3,B:2e-3),C:6e-3);");
    let tree = PhyloTree::from_newick(input_str.as_bytes()).unwrap();
    dbg!(format!("{}", &tree.to_newick()));
    for node in tree.postord_nodes(tree.get_root_id()) {
        dbg!(node.get_weight());
    }
}
#[test]
fn newick_quoted_labels_and_escaping() {
    // Regression for issue #32: an unquoted `;`/`:`/space inside a *quoted*
    // label must not terminate or corrupt parsing (GTDB trees hit this).
    let input_str = String::from(
        "((GB_GCA_929200395.1:0.01646,GB_GCA_929200515.1:0.0207)\
         '100:f__CAKMZJ01; g__CAKMZJ01':0.39515,X:0.1);",
    );
    let tree = PhyloTree::from_newick(input_str.as_bytes()).unwrap();

    let taxa = tree.get_taxa_space().collect::<HashSet<&String>>();
    // Underscores and dots in identifiers stay literal.
    assert!(taxa.contains(&"GB_GCA_929200395.1".to_string()));
    assert!(taxa.contains(&"GB_GCA_929200515.1".to_string()));
    assert!(taxa.contains(&"X".to_string()));
    assert_eq!(tree.get_leaves().count(), 3);
    // The quoted internal label is preserved verbatim, whitespace included.
    assert!(taxa.contains(&"100:f__CAKMZJ01; g__CAKMZJ01".to_string()));

    // A doubled quote inside a quoted label is an escaped single quote.
    let escaped = PhyloTree::from_newick("('a''b':0.1,C:0.2);".as_bytes()).unwrap();
    assert!(escaped.get_taxa_space().any(|t| t == &"a'b".to_string()));
}

#[test]
fn newick_tool_variations() {
    // Comments (rooting marker, NHX, BEAST-style with commas/braces inside the
    // brackets) are skipped; topology and labels are unaffected.
    for input in [
        "[&R] (A,B);",
        "(A[&&NHX:S=human]:0.1,B:0.2);",
        "(A[&rate=0.5,height={1,2}]:0.1,B:0.2);",
        "(\n  A:0.1,\n  B:0.2\n);",
    ] {
        let tree = PhyloTree::from_newick(input.as_bytes()).unwrap();
        let taxa = tree.get_taxa_space().collect::<HashSet<&String>>();
        assert!(taxa.contains(&"A".to_string()) && taxa.contains(&"B".to_string()));
        assert_eq!(tree.get_leaves().count(), 2);
    }

    // IQ-TREE concatenated support (`SH-aLRT/UFBoot`) is one unquoted label.
    let iqtree = PhyloTree::from_newick("((A,B)95.3/98:0.1,C);".as_bytes()).unwrap();
    assert!(iqtree.get_taxa_space().any(|t| t == &"95.3/98".to_string()));

    // Scientific and negative branch lengths parse.
    let weighted = PhyloTree::from_newick("(A:1.5e-3,B:-0.2);".as_bytes()).unwrap();
    assert_eq!(weighted.get_node(2).unwrap().get_weight(), Some(-0.2));
    assert!(weighted.get_node(1).unwrap().get_weight().unwrap() - 1.5e-3 < 1e-9);
}

#[test]
fn newick_malformed_input_does_not_panic() {
    // Each of these must return an error rather than panic.
    for bad in [
        "".as_bytes(),
        ";".as_bytes(),
        "('A,B);".as_bytes(),         // unterminated quote
        "(A[unclosed,B);".as_bytes(), // unterminated comment
        "(A,B));".as_bytes(),         // extra ')'
        "((A,B);".as_bytes(),         // unclosed '('
        "(A]B);".as_bytes(),          // stray ']'
        &[0xff, 0x28],                // invalid UTF-8
    ] {
        assert!(PhyloTree::from_newick(bad).is_err());
    }
    // A trailing root branch length with no terminator is tolerated.
    assert!(PhyloTree::from_newick("(A,B):0.5".as_bytes()).is_ok());
}

#[test]
fn newick_preserves_annotations() {
    // NHX / BEAST-style `[...]` metadata is captured verbatim per node and
    // survives a round trip, so downstream analysis loses no information.
    let input = "(A[&&NHX:S=human]:0.1,B:0.2)[&rate=0.5]:0.3;";
    let tree = PhyloTree::from_newick(input.as_bytes()).unwrap();

    // Leaf A carries its NHX comment; the internal (root) node carries its own.
    let a = tree
        .get_nodes()
        .find(|n| n.get_taxa() == Some(&"A".to_string()))
        .unwrap();
    assert_eq!(a.get_annotation(), Some("[&&NHX:S=human]"));
    let root = tree.get_node(tree.get_root_id()).unwrap();
    assert_eq!(root.get_annotation(), Some("[&rate=0.5]"));
    // A node with no comment has no annotation.
    let b = tree
        .get_nodes()
        .find(|n| n.get_taxa() == Some(&"B".to_string()))
        .unwrap();
    assert_eq!(b.get_annotation(), None);

    // Round trip reproduces the annotations verbatim.
    assert_eq!(tree.to_newick().to_string(), input);
}

#[test]
fn newick_annotation_handlers() {
    use std::sync::Arc;
    let input = "(A[&&NHX:S=human]:0.1,B:0.2);";

    let taxon = |tree: &PhyloTree, name: &str| -> Option<Arc<str>> {
        tree.get_nodes()
            .find(|n| n.get_taxa() == Some(&name.to_string()))
            .unwrap()
            .get_annotation()
            .map(Arc::from)
    };

    // Default keeps the annotation verbatim.
    let kept = PhyloTree::from_newick(input.as_bytes()).unwrap();
    assert_eq!(taxon(&kept, "A").as_deref(), Some("[&&NHX:S=human]"));

    // Discard drops every annotation.
    let dropped = PhyloTree::from_newick_with(input.as_bytes(), DiscardAnnotations).unwrap();
    assert_eq!(taxon(&dropped, "A"), None);

    // A closure can transform the raw text — here, strip the brackets.
    let stripped = PhyloTree::from_newick_with(input.as_bytes(), |raw: &str| {
        Some(Arc::from(raw.trim_start_matches('[').trim_end_matches(']')))
    })
    .unwrap();
    assert_eq!(taxon(&stripped, "A").as_deref(), Some("&&NHX:S=human"));
}

#[test]
fn newick_to_newick_annotation_writers() {
    // Annotations are retained in memory; the caller decides what to write out.
    let input = "(A[&&NHX:S=human]:0.1,B[&keep]:0.2);";
    let tree = PhyloTree::from_newick(input.as_bytes()).unwrap();

    // Default writes everything verbatim (round trip).
    assert_eq!(tree.to_newick().to_string(), input);

    // Discard omits all annotations from the output (tree in memory unchanged).
    assert_eq!(
        tree.to_newick_with(DiscardAnnotations).to_string(),
        "(A:0.1,B:0.2);"
    );

    // A closure can include only some annotations — here, keep `[&keep]` and
    // drop the NHX comment.
    let filtered = tree
        .to_newick_with(|ann: &str| {
            if ann.contains("keep") {
                Some(ann.to_string())
            } else {
                None
            }
        })
        .to_string();
    assert_eq!(filtered, "(A:0.1,B[&keep]:0.2);");

    // The in-memory annotations are untouched by any of the above.
    assert_eq!(tree.to_newick().to_string(), input);
}

#[test]
fn newick_deep_tree_does_not_overflow_stack() {
    // A pectinate tree nested `n` deep would overflow a recursive parser or
    // serialiser; both are iterative, so parsing, serialising, and re-parsing
    // all handle it. Round-tripping preserves the (binary) topology exactly.
    let n = 100_000;
    let mut s = String::with_capacity(n * 10);
    for i in 0..n {
        s.push('(');
        s.push_str(&format!("t{i},"));
    }
    s.push_str(&format!("t{n}"));
    for _ in 0..n {
        s.push(')');
    }
    s.push(';');

    let tree = PhyloTree::from_newick(s.as_bytes()).unwrap();
    assert_eq!(tree.get_leaves().count(), n + 1);

    let serialised = tree.to_newick().to_string();
    assert_eq!(serialised, s);
    let reparsed = PhyloTree::from_newick(serialised.as_bytes()).unwrap();
    assert_eq!(reparsed.get_leaves().count(), n + 1);
}

#[test]
fn read_nexus() {
    let input_str =
        String::from("#NEXUS\n\nBEGIN TREES;\n\tTree tree=((A:1,B:1):1,(C:1,D:1):1);\nEND;");
    let tree = PhyloTree::from_nexus(input_str.clone()).unwrap();
    assert_eq!(
        tree.to_newick().to_string(),
        "((A:1,B:1):1,(C:1,D:1):1);".to_string()
    );
    assert_eq!(tree.to_nexus().unwrap(), input_str);
}
#[test]
fn tree_spr() {
    let input_str = String::from("((A,B),C);");
    let mut tree = PhyloTree::from_newick(input_str.as_bytes()).unwrap();
    dbg!(format!("{}", &tree.to_newick()));
    dbg!(tree.get_nodes().collect_vec());
    let p_tree = tree.prune(1).unwrap();
    dbg!(format!("{}", &tree.to_newick()));
    dbg!(format!("{}", &p_tree.to_newick()));
    tree.graft(p_tree, (0, 4)).unwrap();
    tree.clean();
    dbg!(format!("{}", &tree.to_newick()));
    dbg!(&tree.get_node_parent(4));
    tree.spr((1, 2), (5, 4)).unwrap();
    dbg!(format!("{}", &tree.to_newick()));
}

#[test]
fn tree_nni() {
    let input_str = String::from("(A,(B,(C,D)));");
    let mut tree = PhyloTree::from_newick(input_str.as_bytes()).unwrap();
    dbg!(format!("{}", &tree.to_newick()));
    assert!(tree.nni(4, true).is_ok());
    dbg!(format!("True: {}", &tree.to_newick()));

    let input_str = String::from("(A,(B,(C,D)));");
    let mut tree = PhyloTree::from_newick(input_str.as_bytes()).unwrap();
    dbg!(format!("{}", &tree.to_newick()));
    assert!(tree.nni(4, false).is_ok());
    dbg!(format!("False: {}", &tree.to_newick()));
}

#[test]
fn tree_cluster() {
    let input_str: String = String::from("((A,B),C);");
    let tree = PhyloTree::from_newick(input_str.as_bytes()).unwrap();
    dbg!(&tree.get_cluster(0).collect_vec());
    dbg!(&tree.get_cluster(1).collect_vec());
    let bp = tree.get_bipartition((0, 1));
    dbg!(&bp.0.collect_vec());
    dbg!(&bp.1.collect_vec());
}
#[test]
fn balance_tree() {
    let input_str: String = String::from("(((A,B),C),D);");
    let mut tree = PhyloTree::from_newick(input_str.as_bytes()).unwrap();
    tree.balance_subtree().unwrap();
    dbg!(format!("{}", &tree.to_newick()));
    let input_str: String = String::from("(D,(C,(A,B)));");
    let mut tree = PhyloTree::from_newick(input_str.as_bytes()).unwrap();
    tree.balance_subtree().unwrap();
    dbg!(format!("{}", &tree.to_newick()));
    let input_str: String = String::from("(D,(A,(C,B)));");
    let mut tree = PhyloTree::from_newick(input_str.as_bytes()).unwrap();
    tree.balance_subtree().unwrap();
    dbg!(format!("{}", &tree.to_newick()));
    dbg!(tree.get_nodes().collect_vec());
    dbg!(tree.get_root_id());
}
#[test]
fn induce_tree() {
    let input_str: String = String::from("(((A,B),C),D);");
    let tree = PhyloTree::from_newick(input_str.as_bytes()).unwrap();
    dbg!(format!("{}", &tree.to_newick()));
    let mut x = tree.induce_tree(vec![3, 5, 6]).unwrap();
    x.clean();
    dbg!(x.get_root().get_children());
    dbg!(x.get_nodes().collect_vec());
    dbg!(format!("{}", &x.to_newick()));
}
#[test]
fn median_node() {
    let input_str: String = String::from("(((A,B),C),D);");
    let tree = PhyloTree::from_newick(input_str.as_bytes()).unwrap();
    dbg!(format!("{}", &tree.to_newick()));
    dbg!(tree.get_cluster(tree.get_median_node_id()).collect_vec());
}

#[test]
fn yule() {
    let tree1 = PhyloTree::yule(20);
    dbg!(format!("{}", &tree1.to_newick()));
}

#[test]
fn uniform() {
    let tree1 = PhyloTree::unif(20);
    dbg!(format!("{}", &tree1.to_newick()));
}

#[test]
fn const_lca() {
    let tree = PhyloTree::yule(20);
    dbg!(format!("{}", tree.lca().get_lca_id(vec![1, 10].as_slice())));
}

#[test]
fn contract_tree() {
    fn depth(tree: &PhyloTree, node_id: usize) -> f32 {
        tree.depth(node_id) as f32
    }
    let mut tree = PhyloTree::yule(10);
    dbg!(&tree);
    tree.set_zeta(depth).unwrap();
    println!("{}", tree.to_newick());
    let taxa_subset = vec![
        "1".to_string(),
        "4".to_string(),
        "3".to_string(),
        "7".to_string(),
    ]
    .into_iter()
    .map(|x| tree.get_taxa_node_id(&x).unwrap())
    .collect_vec();
    let new_tree = tree.contract_tree(taxa_subset.as_slice()).unwrap();
    println!("{}", new_tree.to_newick());

    let input_str: String = String::from("(1:1.13,((0:0.93,3:1.40):0.58,(2:1.14,4:1.04)):0.11);");
    let tree: SimpleRootedTree<String, f32, f32> =
        SimpleRootedTree::from_newick(input_str.as_bytes()).unwrap();
    dbg!(&tree);
    let taxa_subset = vec![
        "1".to_string(),
        "0".to_string(),
        "4".to_string(),
        // "7".to_string(),
    ]
    .into_iter()
    .map(|x| tree.get_taxa_node_id(&x).unwrap())
    .collect_vec();

    let new_tree = tree.contract_tree(taxa_subset.as_slice()).unwrap();
    println!("{}", new_tree.to_newick());
}

#[test]
fn cophenetic_dist() {
    fn depth(tree: &PhyloTree, node_id: usize) -> f32 {
        tree.depth(node_id) as f32
    }
    let t1_input_str: String = String::from("((A,B),C);");
    let t2_input_str: String = String::from("(A,(B,C));");
    let mut t1 = PhyloTree::from_newick(t1_input_str.as_bytes()).unwrap();
    let mut t2 = PhyloTree::from_newick(t2_input_str.as_bytes()).unwrap();

    t1.set_zeta(depth).unwrap();
    t2.set_zeta(depth).unwrap();

    assert_eq!(t1.cophen_dist(&t2, 1), 4_f32);

    dbg!(t1.cophen_dist(&t2, 0));
}

#[test]
fn suppress_tree_node() {
    let input_str: String = String::from("(((A,B),C),D);");
    let mut tree = PhyloTree::from_newick(input_str.as_bytes()).unwrap();
    tree.supress_node(2).expect("node id should be valid");
}

#[test]
fn robinson_foulds() {
    let input_str: String = String::from("(((A,B),C),D);");
    let t1 = PhyloTree::from_newick(input_str.as_bytes()).unwrap();
    let input_str: String = String::from("(A,(B,(C,D)));");
    let t2 = PhyloTree::from_newick(input_str.as_bytes()).unwrap();
    assert_eq!(t1.rf(&t2), 0);

    let input_str: String = String::from("(((A,B),C),D);");
    let t1 = PhyloTree::from_newick(input_str.as_bytes()).unwrap();
    let input_str: String = String::from("(A,(D,(C,B)));");
    let t2 = PhyloTree::from_newick(input_str.as_bytes()).unwrap();
    assert_eq!(t1.rf(&t2), 2);

    let input_str: String = String::from("((A:0.1,B:0.2):0.6,(C:0.3,D:0.4):0.5);");
    let t1 = PhyloTree::from_newick(input_str.as_bytes()).unwrap();
    let input_str: String = String::from("((A:0.3,C:0.4):0.5,(B:0.2,D:0.1):0.6);");
    let t2 = PhyloTree::from_newick(input_str.as_bytes()).unwrap();
    assert_eq!(t1.rf(&t2), 2);

    let input_str: String = String::from("(A, ((B, (C, (D, E))), ((F, G), (H, I))));");
    let t1 = PhyloTree::from_newick(input_str.as_bytes()).unwrap();
    let input_str: String = String::from("(A, ((B, (C, (D, (H, I)))), ((F, G), E)));");
    let t2 = PhyloTree::from_newick(input_str.as_bytes()).unwrap();
    assert_eq!(t1.rf(&t2), 8);
}

#[test]
fn cluster_affinity() {
    let input_str: String = String::from("(((A,B),C),D);");
    let t1 = PhyloTree::from_newick(input_str.as_bytes()).unwrap();
    let input_str: String = String::from("(A,(B,(C,D)));");
    let t2 = PhyloTree::from_newick(input_str.as_bytes()).unwrap();
    assert_eq!(t1.ca(&t2), 2);

    let input_str: String = String::from("(((A,B),C),D);");
    let t1 = PhyloTree::from_newick(input_str.as_bytes()).unwrap();
    let input_str: String = String::from("(((A,B),C),D);");
    let t2 = PhyloTree::from_newick(input_str.as_bytes()).unwrap();
    assert_eq!(t1.ca(&t2), 0);

    let input_str: String = String::from("(((A,B),C),D);");
    let t1 = PhyloTree::from_newick(input_str.as_bytes()).unwrap();
    let input_str: String = String::from("(A,(D,(C,B)));");
    let t2 = PhyloTree::from_newick(input_str.as_bytes()).unwrap();
    assert_eq!(t1.ca(&t2), 2);
}

#[test]
fn bipartitions() {
    let input_str: String = String::from("(((A,B),C),D);");
    let t1 = PhyloTree::from_newick(input_str.as_bytes()).unwrap();
    let _bps = t1
        .get_bipartitions_ids()
        .map(|(p1, p2)| {
            (
                p1.map(|x| t1.get_node_taxa(x).cloned().unwrap())
                    .collect_vec(),
                p2.map(|x| t1.get_node_taxa(x).cloned().unwrap())
                    .collect_vec(),
            )
        })
        .collect_vec();

    let input_str: String = String::from("(A, (B, (C, (D, (E, (F, (G, H)))))));");
    let t1 = PhyloTree::from_newick(input_str.as_bytes()).unwrap();
    let _bps = t1
        .get_bipartitions_ids()
        .map(|(p1, p2)| {
            (
                p1.map(|x| t1.get_node_taxa(x).cloned().unwrap())
                    .collect_vec(),
                p2.map(|x| t1.get_node_taxa(x).cloned().unwrap())
                    .collect_vec(),
            )
        })
        .collect_vec();
}

#[test]
#[cfg(feature = "parallel")]
fn compute_norm_parallel() {
    for norm in 1..10 {
        let x = (1..1000).map(|x| x as f32).collect_vec();
        let y = x.clone();
        assert!(
            (PhyloTree::compute_norm(x.into_iter(), norm) - PhyloTree::compute_norm_par(y, norm))
                .abs()
                < 0.1
        );
    }

    let x = (1..3).combinations_with_replacement(2).collect_vec();
    let y = (1..3)
        .combinations_with_replacement(2)
        .par_bridge()
        .map(|x| x[0] + x[1])
        .collect::<Vec<_>>();

    dbg!(x, y);
}

#[test]
#[cfg(feature = "serde")]
fn serde_round_trip() {
    let input_str = String::from("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);");
    let tree = PhyloTree::from_newick(input_str.as_bytes()).unwrap();

    let json = serde_json::to_string(&tree).unwrap();
    let tree2: PhyloTree = serde_json::from_str(&json).unwrap();

    // Same newick output
    assert_eq!(tree.to_newick().to_string(), tree2.to_newick().to_string());

    // Same taxa space
    let mut taxa1: Vec<_> = tree.get_taxa_space().cloned().collect();
    let mut taxa2: Vec<_> = tree2.get_taxa_space().cloned().collect();
    taxa1.sort();
    taxa2.sort();
    assert_eq!(taxa1, taxa2);

    // Same node count
    assert_eq!(tree.num_nodes(), tree2.num_nodes());

    // Taxa lookup still works after deserialization
    assert!(tree2.get_taxa_node(&"A".to_string()).is_some());
    assert!(tree2.get_taxa_node(&"D".to_string()).is_some());
}

#[test]
#[cfg(feature = "parallel")]
fn cophenetic_dist_par() {
    fn depth(tree: &PhyloTree, node_id: usize) -> f32 {
        tree.depth(node_id) as f32
    }
    let t1_input_str: String = String::from("((A,B),C);");
    let t2_input_str: String = String::from("(A,(B,C));");
    let mut t1 = PhyloTree::from_newick(t1_input_str.as_bytes()).unwrap();
    let mut t2 = PhyloTree::from_newick(t2_input_str.as_bytes()).unwrap();

    t1.set_zeta(depth).unwrap();
    t2.set_zeta(depth).unwrap();

    assert_eq!(t1.cophen_dist_par(&t2, 1), 4_f32);
}

// ── OLA tests ──────────────────────────────────────────────────────────────

/// Helper: parse a Newick string, call to_vec, return (taxa, indices).
fn ola_encode(newick: &str) -> (Vec<String>, Vec<i64>) {
    let tree = PhyloTree::from_newick(newick.as_bytes()).unwrap();
    let ola = tree.to_vec();
    (ola.taxa, ola.indices)
}

/// Helper: build an OLATree, decode with from_vec, then re-encode and return
/// (taxa, indices) so we can compare against the original encoding.
fn ola_roundtrip(taxa: Vec<&str>, indices: Vec<i64>) -> (Vec<String>, Vec<i64>) {
    let ola_in = OLATree {
        taxa: taxa.iter().map(|s| s.to_string()).collect(),
        indices,
    };
    let decoded: PhyloTree = PhyloTree::from_vec(ola_in);
    let ola_out = decoded.to_vec();
    (ola_out.taxa, ola_out.indices)
}

// ── to_vec ──────────────────────────────────────────────────────────────────

/// Right-leaning comb (A,(B,C)): sibling of each new leaf is always the
/// immediately preceding leaf, so all entries are non-negative.
#[test]
fn ola_to_vec_right_comb() {
    let (taxa, indices) = ola_encode("(A,(B,C));");
    assert_eq!(taxa, vec!["A", "B", "C"]);
    // i=1: sibling of B = A(0)
    // i=2: sibling of C = B(1)
    assert_eq!(indices, vec![0, 1]);
}

/// Left-leaning comb ((A,B),C): when the sibling is an internal node the
/// entry is negative.
#[test]
fn ola_to_vec_left_comb() {
    let (taxa, indices) = ola_encode("((A,B),C);");
    assert_eq!(taxa, vec!["A", "B", "C"]);
    // i=1: sibling of B = A(0)        → 0
    // i=2: sibling of C = (A,B) node  → -max(μ(A)=0, μ(B)=1) = -1
    assert_eq!(indices, vec![0, -1]);
}

/// Balanced 4-leaf tree ((A,B),(C,D)).
#[test]
fn ola_to_vec_balanced_4() {
    let (taxa, indices) = ola_encode("((A,B),(C,D));");
    assert_eq!(taxa, vec!["A", "B", "C", "D"]);
    // i=1: sibling of B = A(0)           → 0
    // i=2: sibling of C = (A,B) node     → -max(0,1) = -1
    // i=3: sibling of D = C(2)           → 2
    assert_eq!(indices, vec![0, -1, 2]);
}

/// Mixed tree ((A,B),(C,(D,E))): tests a deeper right subtree.
#[test]
fn ola_to_vec_mixed_5() {
    let (taxa, indices) = ola_encode("((A,B),(C,(D,E)));");
    assert_eq!(taxa, vec!["A", "B", "C", "D", "E"]);
    // i=1: sibling B → A(0)                         → 0
    // i=2: sibling C → (A,B): -max(0,1)             → -1
    // i=3: sibling D → C(2)                         → 2
    // i=4: sibling E → D(3)                         → 3
    assert_eq!(indices, vec![0, -1, 2, 3]);
}

/// Balanced 6-leaf tree (((A,B),(C,D)),(E,F)): requires a 2-level internal
/// node index computation.
#[test]
fn ola_to_vec_balanced_6() {
    let (taxa, indices) = ola_encode("(((A,B),(C,D)),(E,F));");
    assert_eq!(taxa, vec!["A", "B", "C", "D", "E", "F"]);
    // i=1: B → A(0)                           → 0
    // i=2: C → (A,B): -max(0,1)              → -1
    // i=3: D → C(2)                           → 2
    // i=4: E → ((A,B),(C,D)): sib=root of ABCD subtree
    //       split by children of ABCD: AB side min=0, CD side min=2
    //       -max(0,2)                          → -2
    // i=5: F → E(4)                           → 4
    assert_eq!(indices, vec![0, -1, 2, -2, 4]);
}

// ── from_vec ────────────────────────────────────────────────────────────────

/// Decode a right-comb OLA vector and verify the resulting tree encodes back
/// to the same OLA representation.
#[test]
fn ola_from_vec_right_comb() {
    let (taxa_out, indices_out) = ola_roundtrip(vec!["A", "B", "C"], vec![0, 1]);
    assert_eq!(taxa_out, vec!["A", "B", "C"]);
    assert_eq!(indices_out, vec![0, 1]);
}

/// Decode a left-comb OLA vector and verify the roundtrip.
#[test]
fn ola_from_vec_left_comb() {
    let (taxa_out, indices_out) = ola_roundtrip(vec!["A", "B", "C"], vec![0, -1]);
    assert_eq!(taxa_out, vec!["A", "B", "C"]);
    assert_eq!(indices_out, vec![0, -1]);
}

/// Decode a balanced 4-leaf tree and verify the roundtrip.
#[test]
fn ola_from_vec_balanced_4() {
    let (taxa_out, indices_out) = ola_roundtrip(vec!["A", "B", "C", "D"], vec![0, -1, 2]);
    assert_eq!(taxa_out, vec!["A", "B", "C", "D"]);
    assert_eq!(indices_out, vec![0, -1, 2]);
}

/// Decode a 6-leaf tree with a 2-level internal node index and verify the
/// roundtrip.
#[test]
fn ola_from_vec_balanced_6() {
    let (taxa_out, indices_out) =
        ola_roundtrip(vec!["A", "B", "C", "D", "E", "F"], vec![0, -1, 2, -2, 4]);
    assert_eq!(taxa_out, vec!["A", "B", "C", "D", "E", "F"]);
    assert_eq!(indices_out, vec![0, -1, 2, -2, 4]);
}

// ── roundtrip (Newick → to_vec → from_vec → to_vec) ────────────────────────

/// For any tree T, to_vec(from_vec(to_vec(T))) must equal to_vec(T).
#[test]
fn ola_roundtrip_newick_3() {
    for newick in &["((A,B),C);", "(A,(B,C));"] {
        let tree = PhyloTree::from_newick(newick.as_bytes()).unwrap();
        let ola1 = tree.to_vec();
        let decoded: PhyloTree = PhyloTree::from_vec(OLATree {
            taxa: ola1.taxa.clone(),
            indices: ola1.indices.clone(),
        });
        let ola2 = decoded.to_vec();
        assert_eq!(ola1.taxa, ola2.taxa, "taxa mismatch for {newick}");
        assert_eq!(ola1.indices, ola2.indices, "indices mismatch for {newick}");
    }
}

#[test]
fn ola_roundtrip_newick_5() {
    for newick in &[
        "((A,B),(C,(D,E)));",
        "(((A,B),C),(D,E));",
        "(A,(B,(C,(D,E))));",
    ] {
        let tree = PhyloTree::from_newick(newick.as_bytes()).unwrap();
        let ola1 = tree.to_vec();
        let decoded: PhyloTree = PhyloTree::from_vec(OLATree {
            taxa: ola1.taxa.clone(),
            indices: ola1.indices.clone(),
        });
        let ola2 = decoded.to_vec();
        assert_eq!(ola1.taxa, ola2.taxa, "taxa mismatch for {newick}");
        assert_eq!(ola1.indices, ola2.indices, "indices mismatch for {newick}");
    }
}

/// LCA of two nodes computed straight from parent links, with no index and no
/// euler walk. Deliberately naive: it is the oracle the precomputed path is
/// checked against, so it must share nothing with it.
fn oracle_lca(tree: &PhyloTree, x: usize, y: usize) -> usize {
    let x_ancestors = tree.node_to_root_ids(x).collect_vec();
    tree.node_to_root_ids(y)
        .find(|a| x_ancestors.contains(a))
        .expect("nodes of one tree share a root")
}

#[test]
fn lca_oracle_rebuilds_correctly_after_mutation() {
    let mut tree = PhyloTree::yule(32);
    let leaves = tree.get_leaf_ids().collect_vec();

    // A fresh oracle agrees with the naive walk before anything moves.
    {
        let oracle = tree.lca();
        for pair in leaves.iter().combinations(2) {
            let (x, y) = (*pair[0], *pair[1]);
            assert_eq!(oracle.get_lca_id(&[x, y]), oracle_lca(&tree, x, y));
        }
    } // The oracle's borrow of `tree` ends here, releasing it for mutation.

    // Move a subtree. The borrow checker guarantees no oracle is alive across
    // this mutation -- a stale index is unrepresentable now, which is why the
    // old runtime `invalidate_lca_index` machinery is gone. (Uncommenting a
    // query on an oracle built before this line would fail to compile.)
    let moved = leaves[0];
    let target = leaves[leaves.len() - 1];
    let from = (tree.get_node_parent_id(moved).unwrap(), moved);
    let to = (tree.get_node_parent_id(target).unwrap(), target);
    tree.spr(from, to).expect("spr on a yule tree");

    // A newly built oracle describes the tree as it is now, not as it was.
    let oracle = tree.lca();
    for pair in leaves.iter().combinations(2) {
        let (x, y) = (*pair[0], *pair[1]);
        if !tree.contains_node(x) || !tree.contains_node(y) {
            continue;
        }
        assert_eq!(
            oracle.get_lca_id(&[x, y]),
            oracle_lca(&tree, x, y),
            "rebuilt LCA for ({x}, {y}) after the topology changed"
        );
    }
}

#[test]
fn subtree_extraction_leaves_no_phantom_nodes() {
    let tree = PhyloTree::yule(20);
    // A non-root internal node: its extracted subtree excludes the original
    // root, so the old `Self::new()` placeholder would linger unreachable.
    let root = tree
        .get_node_ids()
        .find(|&id| id != tree.get_root_id() && !tree.is_leaf(id))
        .expect("a yule tree has a non-root internal node");
    let child = tree.subtree(root).unwrap();

    // The arena must hold exactly the nodes reachable from the subtree root.
    let reachable = child.dfs(child.get_root_id()).count();
    assert_eq!(
        child.get_node_ids().count(),
        reachable,
        "extracted subtree arena holds unreachable phantom nodes"
    );
}

#[test]
fn restricted_subtree_oracle_matches_a_fresh_build() {
    let tree = PhyloTree::yule(30);
    let parent_oracle = tree.lca();

    // Every internal node roots a real subtree; restrict at each of them.
    let internal = tree
        .get_node_ids()
        .filter(|&id| !tree.is_leaf(id))
        .collect_vec();
    assert!(!internal.is_empty());

    for root in internal {
        let child = tree.subtree(root).unwrap();
        let restricted = parent_oracle.restrict_to_subtree(&child);
        let fresh = child.lca();

        let child_nodes = child.get_node_ids().collect_vec();
        let child_leaves = child.get_leaf_ids().collect_vec();
        for pair in child_leaves.iter().combinations(2) {
            let q = [*pair[0], *pair[1]];
            let sliced = restricted.get_lca_id(&q);
            assert_eq!(
                sliced,
                fresh.get_lca_id(&q),
                "sliced oracle disagrees with a fresh build (subtree root {root})"
            );
            assert_eq!(sliced, oracle_lca(&child, q[0], q[1]));
            // LCA is invariant under subtree extraction, so the parent oracle
            // must give the same answer for nodes retained in the child.
            assert_eq!(sliced, parent_oracle.get_lca_id(&q));
        }

        // Depths must be re-based so the subtree root sits at depth 0.
        for id in child_nodes {
            assert_eq!(
                restricted.get_node_depth(id),
                RootedTree::get_node_depth(&child, id),
                "sliced depth mismatch for node {id} (subtree root {root})"
            );
        }
    }
}

#[test]
fn removing_a_node_removes_its_taxon() {
    let mut tree = PhyloTree::from_newick("((A,B),(C,D));".as_bytes()).unwrap();
    let before = tree.num_taxa();
    let a = tree.get_taxa_node_id(&"A".to_string()).unwrap();

    tree.remove_node(a);

    assert_eq!(
        tree.num_taxa(),
        before - 1,
        "taxa count must drop when a labelled node is removed"
    );
    assert_eq!(
        tree.get_taxa_node_id(&"A".to_string()),
        None,
        "a removed node's taxon must not still resolve"
    );
}

#[test]
fn from_nodes_rebuilds_the_taxa_map() {
    let tree = PhyloTree::from_newick("((A,B),(C,D));".as_bytes()).unwrap();

    // Reassemble an identical arena from the tree's own nodes and hand it back
    // to from_nodes -- the taxa map must come out populated with no manual
    // re-registration pass.
    let max_id = tree.get_node_ids().max().unwrap();
    let mut nodes: Vec<Option<PhyloNode>> = (0..=max_id).map(|_| None).collect();
    for id in tree.get_node_ids() {
        nodes[id] = Some(tree.get_node(id).unwrap().clone());
    }

    let rebuilt = PhyloTree::from_nodes(nodes, tree.get_root_id());

    assert_eq!(
        rebuilt.num_taxa(),
        tree.num_taxa(),
        "from_nodes must rebuild the taxa map, not leave it empty"
    );
    for taxon in ["A", "B", "C", "D"] {
        assert_eq!(
            rebuilt.get_taxa_node_id(&taxon.to_string()),
            tree.get_taxa_node_id(&taxon.to_string()),
            "taxon {taxon} must resolve to the same node after from_nodes"
        );
    }
}

#[test]
fn relabelling_a_node_does_not_grow_the_taxa_map() {
    let mut tree = PhyloTree::from_newick("((A,B),(C,D));".as_bytes()).unwrap();
    let before = tree.num_taxa();
    let a = tree.get_taxa_node_id(&"A".to_string()).unwrap();

    tree.set_node_taxa(a, Some("Z".to_string()));

    assert_eq!(
        tree.num_taxa(),
        before,
        "relabelling replaces a taxon, it does not add one"
    );
    assert_eq!(tree.get_taxa_node_id(&"A".to_string()), None);
    assert_eq!(tree.get_taxa_node_id(&"Z".to_string()), Some(a));
}

#[test]
fn binary_nodes_do_not_over_reserve_children() {
    // Vec's first push jumps to capacity 4, which for a bifurcating node is
    // twice what it will ever use. A yule tree is entirely leaves and
    // bifurcations, so every node should hold exactly what it uses.
    let tree = PhyloTree::yule(256);
    for id in tree.get_node_ids() {
        let node = tree.get_node(id).unwrap();
        let len = node.get_children().len();
        let capacity = node.heap_size() / std::mem::size_of::<usize>();
        assert_eq!(
            capacity, len,
            "node {id} holds capacity for {capacity} children but has {len}"
        );
    }
}

proptest::proptest! {
    /// The cursor must never hand out a slot that is already live, and must
    /// always hand out the lowest vacant one -- callers depend on both.
    #[test]
    fn next_id_is_the_lowest_vacant_slot(
        taxa in 4usize..40,
        removals in proptest::collection::vec(0usize..40, 0..12),
    ) {
        let mut tree = PhyloTree::yule(taxa);

        for r in removals {
            let live = tree.get_node_ids().collect_vec();
            if live.len() <= 2 {
                break;
            }
            // Never remove the root: it is not a normal arena slot.
            let victim = live[r % live.len()];
            if victim == tree.get_root_id() {
                continue;
            }
            tree.remove_node(victim);

            let next = tree.next_id();
            proptest::prop_assert!(
                !tree.contains_node(next),
                "next_id returned {next}, which is live"
            );
            // The scan-from-zero answer the cursor is standing in for.
            let expected = (0..)
                .find(|id| !tree.contains_node(*id))
                .unwrap();
            proptest::prop_assert_eq!(next, expected, "next_id skipped a lower vacant slot");
        }
    }

    /// The oracle's constant-time query must agree with plain parent-link
    /// walking, on any tree, for any pair.
    #[test]
    fn precomputed_lca_matches_the_naive_walk(taxa in 4usize..40) {
        let tree = PhyloTree::yule(taxa);
        let oracle = tree.lca();
        let leaves = tree.get_leaf_ids().collect_vec();
        for pair in leaves.iter().combinations(2) {
            let (x, y) = (*pair[0], *pair[1]);
            proptest::prop_assert_eq!(
                oracle.get_lca_id(&[x, y]),
                oracle_lca(&tree, x, y),
                "indexed LCA disagrees with the walk for ({}, {})", x, y
            );
        }
    }

    /// A reused oracle and the one-shot `EulerWalk::get_lca_id` fallback (which
    /// builds a throwaway oracle per call) are two paths to the same query and
    /// must not disagree.
    #[test]
    fn indexed_and_unindexed_lca_agree(taxa in 4usize..40) {
        let tree = PhyloTree::yule(taxa);
        let oracle = tree.lca();

        let leaves = tree.get_leaf_ids().collect_vec();
        for pair in leaves.iter().combinations(2) {
            let q = [*pair[0], *pair[1]];
            proptest::prop_assert_eq!(oracle.get_lca_id(&q), tree.get_lca_id(&q));
        }
    }
}

/// `unweight` must clear every edge weight even when the arena has holes. The
/// old body filtered for empty slots and unwrapped them, so it panicked on the
/// first hole and cleared nothing on a full arena.
#[test]
fn unweighting_clears_every_edge_weight() {
    let mut tree = PhyloTree::from_newick("((A:0.1,B:0.2):0.3,C:0.6);".as_bytes()).unwrap();
    let a_leaf = tree.get_taxa_node_id(&"A".to_string()).unwrap();
    tree.delete_node(a_leaf); // leave a vacant slot in the arena
    tree.unweight();
    for node in tree.get_nodes() {
        assert!(
            node.get_weight().is_none(),
            "node {} still weighted after unweight",
            node.get_id()
        );
    }
}

/// `to_nexus_file` must accept a `.nex` path; it previously asserted `.nwk`
/// (copy-pasted from the newick writer) and panicked on any nexus filename.
#[test]
fn to_nexus_file_accepts_a_nex_path() {
    let tree = PhyloTree::from_newick("((A:1,B:1):1,(C:1,D:1):1);".as_bytes()).unwrap();
    let path = std::env::temp_dir().join("phylo_nexus_roundtrip.nex");
    tree.to_nexus_file(&path).unwrap();
    let contents = read_to_string(&path).unwrap();
    assert!(contents.starts_with("#NEXUS"));
    std::fs::remove_file(&path).ok();
}

/// Malformed Nexus input must return an error, not panic on an out-of-bounds
/// index. Each case below used to crash the parser.
#[test]
fn from_nexus_rejects_malformed_input() {
    // Empty input: no first line to compare against the header.
    assert!(PhyloTree::from_nexus(String::new()).is_err());
    // Header present but no trees block at all.
    assert!(PhyloTree::from_nexus("#NEXUS\n".to_string()).is_err());
    // A trees block whose entry has no `=` to split on.
    assert!(PhyloTree::from_nexus("#NEXUS\nBEGIN TREES;\n\tTree tree;\nEND;".to_string()).is_err());
    // Wrong header is still rejected.
    assert!(PhyloTree::from_nexus("NOT A NEXUS FILE".to_string()).is_err());
}
