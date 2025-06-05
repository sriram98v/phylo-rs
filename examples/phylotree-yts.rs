use std::{fs, env};
use std::error::Error;
use std::time::Instant;
use phylotree::generate_yule;
use phylotree::tree::Tree;


fn main() -> Result<(), Box<dyn Error>> {
    let args = env::args().collect::<Vec<_>>();

    let input_str: String = fs::read_to_string(args[1].clone())?;

    let tree = Tree::from_newick(&input_str)?;
    let n_leaves = tree.n_leaves();

    let now = Instant::now();

    let yule_tree = generate_yule(n_leaves, false, phylotree::distr::Distr::Uniform).unwrap();

    let elapsed = now.elapsed();
    println!("{:?}", yule_tree.to_newick());
    
    println!("Internal time: {:.7?}", elapsed.as_secs_f64());

    Ok(())
}
