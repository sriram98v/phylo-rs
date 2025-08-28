#include "compact_tree.h"
#include <chrono>

int main(int argc, char** argv) {
    
    // Start the clock.
    compact_tree tree(argv[1]);

    auto start = std::chrono::high_resolution_clock::now();
    

    CT_NODE_T curr_node;
    std::tuple<const std::string*, CT_LENGTH_T, CT_NODE_T, const std::vector<CT_NODE_T>*> curr_data;
    const std::string* curr_label_ptr;
    CT_LENGTH_T curr_length;
    CT_NODE_T curr_parent;
    const std::vector<CT_NODE_T>* curr_children_ptr;
    size_t curr_children_size;
    size_t curr_child_ind;
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