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
use std::time::{Duration, Instant};
use rand::prelude::IteratorRandom;

fn main(){
    let matches = Command::new("Generalized suffix tree")
        .version("1.0")
        .author("Sriram Vijendran <vijendran.sriram@gmail.com>")
        .subcommand(Command::new("cophen-dist-repr")
            .about("Build suffix tree index from reference fasta file")
            .arg(arg!(-k --norm <NORM> "nth norm")
                .value_parser(clap::value_parser!(u32))
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
        .about("CLI tool for quick tree operations")
        .get_matches();

        match matches.subcommand(){
            Some(("cophen-dist-repr",  sub_m)) => {            
                fn depth(tree: &SimpleRootedTree, node_id: <SimpleRootedTree as RootedTree>::NodeID)->f32
                {
                    EulerWalk::get_node_depth(tree, node_id) as f32
                }
            
                let norm = sub_m.get_one::<u32>("norm").expect("required");
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

                // dbg!(&t1);
                // dbg!(&t2);

                println!("Computing runtime");
                let mean_dist = (0..*num_trees).map(|_| {
                        let taxa_set = t1.get_taxa_space();
                        let now = Instant::now();
                        dbg!(t1.cophen_dist_naive_by_taxa(&t2, *norm, taxa_set.clone()));
                        return now.elapsed();
                    }).sum::<Duration>()/(*num_trees as u32);
                
                println!("Mean time: {:?}", mean_dist);            
            },
            _ => {
                println!("No option selected! Refer help page (-h flag)");
            }
        }
}