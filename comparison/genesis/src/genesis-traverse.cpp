#include "genesis/genesis.hpp"

#include <chrono>

using namespace genesis;
using namespace genesis::tree;

int main( int argc, char** argv )
{
    // Get input.
    if (argc != 3) {
        throw std::runtime_error( "Need to provide a newick tree file." );
    }
    auto const infile = std::string( argv[1] );

    auto const tree = CommonTreeNewickReader().read( utils::from_file( infile ));

    // Start the clock.
    std::cout << "Start reading" << utils::current_time() << "\n";
    auto const start = std::chrono::steady_clock::now();

    // Run, Forrest, Run!

    for( auto it : postorder( tree ) ) {
        std::cout << it.node().data<CommonNodeData>().name << " ";
    }
    std::cout << std::endl;

    // Stop the clock
    auto const duration = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now() - start
    );
    std::cout << "Finished reading " << utils::current_time() << "\n";
    double const elapsed_secs = double(duration.count()) / 1000.0;
    std::cout << "Internal time: " << elapsed_secs << "\n";

    // Check output
    std::cout << "Leaves: " << leaf_node_count( tree ) << "\n";
    return 0;
}
