use phylo::prelude::*;
use std::{fs, env};
use std::error::Error;
use std::time::Instant;


fn main() -> Result<(), Box<dyn Error>> {
    let args = env::args().collect::<Vec<_>>();

    let input_str: String = fs::read_to_string(args[1].clone())?;

    let mut tree = PhyloTree::from_newick(&input_str.as_bytes())?;
    tree.precompute_constant_time_lca();
    let now = Instant::now();

    let lca_id = tree.get_lca_id(vec![10, 20].as_slice());

    let elapsed = now.elapsed();
    println!("{}", &lca_id);
    println!("Internal time: {:.7?}", elapsed.as_secs_f64());

    Ok(())
}
