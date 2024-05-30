extern crate clap;

use phylo::node::simple_rnode::{RootedTreeNode, RootedWeightedNode};
use phylo::tree::io::Newick;
use phylo::*;
use clap::{arg, Command};
use phylo::tree::{SimpleRootedTree, ops::CopheneticDistance};
use itertools::Itertools;
use phylo::iter::node_iter::{Ancestors, EulerWalk};
use phylo::tree::simple_rtree::{RootedTree, RootedMetaTree};
use phylo::tree::distances::PathFunction;
use tree::stats::PhylogeneticDiversity;
use std::time::{Duration, Instant};
use std::fs::File;
use std::io::Read;
use rand::prelude::IteratorRandom;

fn main(){
    let matches = Command::new("Phylogenetics Rust")
        .version("1.0")
        .author("Sriram Vijendran <vijendran.sriram@gmail.com>")
        .subcommand(Command::new("cophen-dist-repr")
            .about("Build suffix tree index from reference fasta file")
            .arg(arg!(-k --norm <NORM> "nth norm")
                .value_parser(clap::value_parser!(usize))
            )
            .arg(arg!(-n --num_trees <NUM_TREES> "number of trees")
                .required(true)
                .value_parser(clap::value_parser!(usize))
            )
            .arg(arg!(-x --num_taxa <NUM_TAXA> "number of taxa")
                .required(true)
                .value_parser(clap::value_parser!(usize))
            )
            .arg(arg!(-t --threads <THREADS> "number of threads")
                .value_parser(clap::value_parser!(usize))
            )
        )
        .subcommand(Command::new("PD")
            .about("Compute Phylogenetic Diversity")
            .subcommand(Command::new("min")
                .about("Compute minPD")
                .arg(arg!(-f --file <TREE_FILE> "Input Tree File")
                    .required(true)
                    .value_parser(clap::value_parser!(String))
                )
                .arg(arg!(-n --num_taxa <NUM_TAXA> "Input Tree File")
                    .required(true)
                    .value_parser(clap::value_parser!(usize))
                )
            )
        )
        .about("CLI tool for quick tree operations")
        .get_matches();

        match matches.subcommand(){
            Some(("cophen-dist-repr",  sub_m)) => {            
                fn depth(tree: &SimpleRootedTree, node_id: <SimpleRootedTree as RootedTree>::NodeID)->f32
                {
                    EulerWalk::get_node_depth(tree, node_id) as f32
                }
            
                let norm = sub_m.get_one::<usize>("norm").expect("required");
                let num_trees = sub_m.get_one::<usize>("num_trees").expect("required");
                let num_taxa = sub_m.get_one::<usize>("num_taxa").expect("required");
                println!("Number of trees: {}", num_trees);
                println!("Number of taxa per tree: {}", num_taxa);
                println!("Norm: {}", norm);
                let num_threads: usize = *sub_m.get_one::<usize>("threads").unwrap_or(&1);
            
                let mut t1 = SimpleRootedTree::yule(*num_taxa).unwrap();
                let mut t2 = SimpleRootedTree::yule(*num_taxa).unwrap();

                t1.precompute_constant_time_lca();
                t2.precompute_constant_time_lca();

                // dbg!(t1.get_nodes().map(|x| (x.get_parent(), x.get_id(), x.get_children().collect_vec(), x.get_weight())).collect_vec());

                t1.set_zeta(depth);
                t2.set_zeta(depth);

                // dbg!(t1.get_nodes().map(|x| (x.get_id(), x.get_weight())).collect_vec());

                println!("{}", t1.to_newick());
                println!("{}", t2.to_newick());

                // dbg!(&t1);
                // dbg!(&t2);

                println!("Computing runtime");
                let mean_dist = (0..*num_trees).map(|_| {
                        let taxa_set = t1.get_taxa_space();
                        let now = Instant::now();
                        dbg!(t1.cophen_dist_naive_by_taxa(&t2, *norm, taxa_set.clone()));
                        dbg!(t1.cophen_dist(&t2, *norm));
                        return now.elapsed();
                    }).sum::<Duration>()/(*num_trees as u32);
                
                println!("Mean time: {:?}", mean_dist);            
            },
            Some(("PD", sub_m)) => {
                match sub_m.subcommand(){
                    Some(("min", min_pd)) => {
                        let mut tree_file = File::open(min_pd.get_one::<String>("file").expect("required")).unwrap();
                        let num_taxa = min_pd.get_one::<usize>("num_taxa").expect("required");
                        let mut trees = String::new();

                        tree_file.read_to_string(&mut trees).unwrap();
                        let tree_string = trees.split("\n").collect_vec()[0];
                        let mut tree = SimpleRootedTree::from_newick(tree_string.as_bytes());
                        tree.precompute_minPDs();
                        println!("minPD: {}\nnormalized minPD: {}", tree.get_minPD(num_taxa.clone()), tree.get_norm_minPD(num_taxa.clone()));
                        // dbg!("{}", tree);
                    }
                    _ => println!("No valid PD metric chosen! Refer help page (-h flag)")
                }
            },
            _ => {
                println!("No option selected! Refer help page (-h flag)");
            }
        }
}