use phylo::prelude::*;
use std::{fs, env};
use std::error::Error;
use std::time::Instant;
use rand::{seq::IteratorRandom, thread_rng};



fn main() -> Result<(), Box<dyn Error>> {
    let args = env::args().collect::<Vec<_>>();

    let input_str: String = fs::read_to_string(args[1].clone())?;

    let mut rng = thread_rng();

    let mut tree = PhyloTree::from_newick(&input_str.as_bytes())?;
    let taxa_set = tree.get_leaf_ids().collect::<Vec<_>>();
    let ntaxa = taxa_set.len();
    let taxa_subset = taxa_set
        .into_iter()
        .choose_multiple(&mut rng, ((ntaxa as f32)*0.05) as usize);
    tree.precompute_constant_time_lca();

    let now = Instant::now();

    let subtree = tree.contract_tree(taxa_subset.as_slice()).unwrap();

    let elapsed = now.elapsed();
    println!("{}", &subtree.to_newick());
    println!("Internal time: {:.7?}", elapsed.as_secs_f64());

    Ok(())
}
