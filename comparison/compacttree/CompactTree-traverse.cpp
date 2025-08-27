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
        curr_data = tree[curr_node];
        curr_label_ptr = std::get<0>(curr_data);
        curr_length = std::get<1>(curr_data);
        curr_parent = std::get<2>(curr_data);
        curr_children_ptr = std::get<3>(curr_data);
        curr_children_size = curr_children_ptr->size();
        std::cout << "- Node " << curr_node << std::endl;
        if(!(curr_label_ptr->empty())) {
            std::cout << "  - Label: " << (*curr_label_ptr) << std::endl;
        }
        std::cout << "  - Edge Length: " << curr_length << std::endl;
        if(curr_parent != NULL_NODE) {
            std::cout << "  - Parent: Node " << curr_parent << std::endl;
        }
        if(curr_children_size != 0) {
            std::cout << "  - Children: {Node " << (*curr_children_ptr)[0];
            for(curr_child_ind = 1; curr_child_ind < curr_children_size; ++curr_child_ind) {
                std::cout << ", Node " << (*curr_children_ptr)[curr_child_ind];
            }
            std::cout << '}' << std::endl;
        }
    }

    auto const duration = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - start
    );
    double const elapsed_secs = double(duration.count()) / 1000.0;
    std::cout << "Internal time: " << elapsed_secs << "\n";


    return 0;
}