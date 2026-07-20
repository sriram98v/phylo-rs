use phylo::prelude::*;
use std::error::Error;
use std::time::Instant;
use std::{env, fs};

fn main() -> Result<(), Box<dyn Error>> {
    let args = env::args().collect::<Vec<_>>();

    let input_str: String = fs::read_to_string(args[1].clone())?;

    let tree = PhyloTree::from_newick(input_str.as_bytes())?;
    let oracle = tree.lca();
    let now = Instant::now();

    let lca_id = oracle.get_lca_id(vec![10, 20].as_slice());

    let elapsed = now.elapsed();
    println!("{}", &lca_id);
    println!("Internal time: {:.7?}", elapsed.as_secs_f64());

    Ok(())
}
