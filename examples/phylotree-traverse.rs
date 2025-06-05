use std::{fs, env};
use std::error::Error;
use std::time::Instant;
use phylotree::tree::Tree;


fn main() -> Result<(), Box<dyn Error>> {
    let args = env::args().collect::<Vec<_>>();

    let input_str: String = fs::read_to_string(args[1].clone())?;

    let now = Instant::now();
    let tree = Tree::from_newick(&input_str)?;
    let root = tree.get_root().unwrap();
    let postorder = tree.postorder(&root).unwrap()
        .iter()
        .map(|node_id| tree.get(node_id).unwrap().name.clone().unwrap_or("x".to_string()))
        .collect::<Vec<_>>();
    let elapsed = now.elapsed();
    println!("{:?}", &postorder);
    
    println!("Internal time: {:.7?}", elapsed.as_secs_f64());

    Ok(())
}
