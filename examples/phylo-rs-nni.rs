use phylo::prelude::*;
use std::{fs, env};
use std::error::Error;
use std::time::Instant;


fn main() -> Result<(), Box<dyn Error>> {
    let args = env::args().collect::<Vec<_>>();

    let input_str: String = fs::read_to_string(args[1].clone())?;

    let mut tree = PhyloTree::from_newick(&input_str.as_bytes())?;
    let leaf_edges = tree.get_leaf_ids().map(|l_id| (tree.get_node_parent_id(l_id).unwrap(), l_id)).collect::<Vec<_>>();
    let node_id = leaf_edges[0].0;


    
    let now = Instant::now();

    let success = tree.nni(node_id, true);

    let elapsed = now.elapsed();
    println!("{}", &success.is_ok());
    println!("Internal time: {:.7?}", elapsed.as_secs_f64());

    Ok(())
}
