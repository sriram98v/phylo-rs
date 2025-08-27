use std::{fs, env};
use std::error::Error;
use std::time::Instant;
use phylotree::tree::Tree;


fn main() -> Result<(), Box<dyn Error>> {
    let args = env::args().collect::<Vec<_>>();

    let input_str: String = fs::read_to_string(args[1].clone())?;

    let tree: Tree = Tree::from_newick(&input_str)?;
    let now = Instant::now();

    let ancestor = tree.get_common_ancestor(
        &tree.get_by_name("Tip10").unwrap().id,
        &tree.get_by_name("Tip20").unwrap().id,
    ).unwrap();
    

    let elapsed = now.elapsed();
    println!("{:?}", &ancestor);
    println!("Internal time: {:.7?}", elapsed.as_secs_f64());

    Ok(())
}
