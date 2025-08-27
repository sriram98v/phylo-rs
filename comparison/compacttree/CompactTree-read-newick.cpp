#include "compact_tree.h"
#include <chrono>

int main(int argc, char** argv) {

    compact_tree tree(argv[1]);

    auto start = std::chrono::high_resolution_clock::now();

    auto const duration = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - start
    );

    double const elapsed_secs = double(duration.count()) / 1000.0;
    std::cout << "Internal time: " << elapsed_secs << "\n";

    return 0;
}