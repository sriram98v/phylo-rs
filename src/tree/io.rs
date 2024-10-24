#![allow(clippy::needless_lifetimes)]

use itertools::Itertools;
use std::ffi::OsStr;
use std::fmt::Display;
use std::path::Path;
use std::{fs, io};

use crate::prelude::*;

/// Enum to track block of Nexus file. This enum can be extended in the future to include new blocks for different use cases.
pub enum NexusBlock {
    /// Tree block
    TREE,
    /// Miscellaneous block to be ignored
    NONE,
}

/// A trait descibing Newick encoding of a tree.
pub trait Newick: RootedTree {
    /// Creates a new tree using a Newick string
    fn from_newick(newick_str: &[u8]) -> std::io::Result<Self>;

    /// Encodes a subtree starting from a node as a Newick string
    fn subtree_to_newick(&self, node_id: TreeNodeID<Self>) -> impl Display;

    /// Encodes a tree as a Newick string
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
        if file_lines[0] != "#NEXUS" {
            Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                NexusError::InvalidHeader,
            ))
        } else {
            let mut tree_block = "".to_string();
            let mut curr_block = NexusBlock::NONE;
            for line in file_lines {
                let line_words = line
                    .split_ascii_whitespace()
                    .map(|x| x.to_ascii_lowercase())
                    .collect_vec();
                if line_words.is_empty() {
                    continue;
                } else if line_words[0] == *"begin" {
                    curr_block = match line_words[1].as_str() {
                        "trees;" => NexusBlock::TREE,
                        _ => NexusBlock::NONE,
                    };
                } else if line_words[0] == "end;" {
                    curr_block = NexusBlock::NONE;
                } else {
                    match curr_block {
                        NexusBlock::TREE => tree_block.push_str(line),
                        NexusBlock::NONE => {}
                    };
                }
            }
            let first_tree = tree_block.split(';').collect_vec()[0]
                .split('=')
                .collect_vec()[1]
                .split_whitespace()
                .collect::<String>();
            Self::from_newick(format!("{first_tree};").as_bytes())
        }
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
        assert!(p.extension() == Some(OsStr::new("nwk")));
        fs::write(p, self.to_nexus()?.as_bytes())
    }
}
