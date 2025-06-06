use phylo::prelude::*;
use std::{fs, env};
use std::error::Error;
use std::time::Instant;


fn main() -> Result<(), Box<dyn Error>> {
    let args = env::args().collect::<Vec<_>>();

    let t1_str: String = fs::read_to_string(args[1].clone())?;
    let t2_str: String = fs::read_to_string(args[2].clone())?;

    let t1 = PhyloTree::from_newick(&t1_str.as_bytes())?;
    let t2 = PhyloTree::from_newick(&t2_str.as_bytes())?;
    
    let now = Instant::now();

    let rf = t1.rf(&t2);

    let elapsed = now.elapsed();
    println!("{}", &rf);
    println!("Internal time: {:.7?}", elapsed.as_secs_f64());

    Ok(())
}
