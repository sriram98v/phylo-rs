use std::{fs, env};
use std::error::Error;
use std::time::Instant;
use phylotree::tree::Tree;


fn main() -> Result<(), Box<dyn Error>> {
    let args = env::args().collect::<Vec<_>>();

    let t1_str: String = fs::read_to_string(args[1].clone())?;
    let t2_str: String = fs::read_to_string(args[2].clone())?;

    let now = Instant::now();
    let tree1 = Tree::from_newick(&t1_str)?;
    let tree2 = Tree::from_newick(&t2_str)?;

    let rf = tree1.robinson_foulds(&tree2).unwrap();

    let elapsed = now.elapsed();

    println!("{:?}", &rf);
    
    println!("Internal time: {:.7?}", elapsed.as_secs_f64());

    Ok(())
}
