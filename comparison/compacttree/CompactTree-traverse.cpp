#include "compact_tree.h"
#include <chrono>

int main(int argc, char** argv) {
    
    // Start the clock.
    compact_tree tree(argv[1]);

    auto start = std::chrono::high_resolution_clock::now();
    

    CT_NODE_T curr_node;
    compact_tree::preorder_iterator it_end = tree.preorder_end();
    for(compact_tree::preorder_iterator it = tree.preorder_begin(); it != it_end; ++it) {
        curr_node = *it;
        std::cout << "- Node " << curr_node << std::endl;
    }

    auto const duration = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - start
    );
    double const elapsed_secs = double(duration.count()) / 1000.0;
    std::cout << "Internal time: " << elapsed_secs << "\n";

    return 0;
}