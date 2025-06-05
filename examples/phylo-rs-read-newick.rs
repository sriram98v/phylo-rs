use phylo::prelude::*;
use std::{fs, env};
use std::error::Error;
use std::time::Instant;


fn main() -> Result<(), Box<dyn Error>> {
    let args = env::args().collect::<Vec<_>>();

    let input_str: String = fs::read_to_string(args[1].clone())?;

    let now = Instant::now();
    let tree = PhyloTree::from_newick(&input_str.as_bytes())?;
    
    let elapsed = now.elapsed();
    println!("{}", &tree.to_newick());
    println!("Internal time: {:.5?}", elapsed.as_secs_f64());

    Ok(())
}
