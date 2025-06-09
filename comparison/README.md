# Comparison Study
Scripts and tests for the comparison study of Phylo-rs.

The scripts and tests provided here are used for evaluating the runtime and memory of Phylo-rs for some typical tasks, compared to other bioinformatics libraries with similar functionality.

## Test Cases

Memory Test Cases:
 - `read_newick`: Read large trees in newick format.

Runtime Test Cases:
 - `Tree traversal`: traversing a tree in postfix order.
 - `Least Common Ancestor retrieval`: Finding the LCA for a set of taxa in a rooted tree.
 - `Nearest Neighbor Interchange`: Performing an NNI operation on a rooted tree.
 - `Yule tree Simulation`: Simulating a tree under the Yule evolutionary model.
 - `Tree Contraction`: Contracting a tree to an arbitrary subset of leaves.
 - `Robinson Foulds metric`: Computing the Robinson Foulds distance for a pair of given trees.

 
## Software

In order to run the tests, the following software is required:

 - [Genesis 0.33.0](https://github.com/lczech/genesis/releases/tag/v0.33.0)
 - [ape 5.8-1](https://cran.r-project.org/web/packages/ape/index.html)
 - [DendroPy 5.0.8](https://jeetsukumaran.github.io/DendroPy/)
 - [CompactTree 1.0.0](https://github.com/niemasd/CompactTree/releases/tag/1.0.0)
 - [Gotree 0.4.5](https://github.com/evolbioinfo/gotree/releases/tag/v0.4.5)
 - [TreeSwift 1.1.45](https://github.com/niemasd/TreeSwift/releases/tag/v1.1.45)
 - [phylo 2.0.0](https://crates.io/crates/phylo)
 - [phylotree 0.1.3](https://crates.io/crates/phylotree)
 
## Instructions

To build all the executables used, run the ```build.sh``` script. Be sure to install all the dependencies for Genesis and CompactTree before running the build script. Following the build script use the ```measure-mem.sh``` and ```measure-time.sh``` scripts to run the memory and runtime analysis, respectively. Be sure to install all the software listed in the [Software](#software) section, as gotree will be required to simulate the trees used in the comparative analysis. Finally, use the ```plot-memory-scalability.py``` and ```plot-runtime-scalability.py``` scripts to generate plots. Note that the requirements in the ```requirements.txt``` file should be installed before running any of the scripts in this directory. It is recommended to use a Python virtual environment when installing these dependencies, instructions for which can be found [here](https://docs.python.org/3/library/venv.html)
