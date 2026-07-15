#[cfg(feature = "non_crypto_hash")]
use fxhash::FxHashMap as HashMap;
#[cfg(not(feature = "non_crypto_hash"))]
use std::collections::HashMap;

use itertools::Itertools;
use phylo::prelude::*;
use std::fs::{read_to_string, File};
use std::io::Write;

fn main() {
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
        println!("Clade: {}", clade);
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
        let out = format!(
            "{}: {}\n",
            clade,
            pds.iter().map(|x| x.to_string()).join(",")
        );
        output_file.write_all(out.as_bytes()).unwrap()
    }
}
