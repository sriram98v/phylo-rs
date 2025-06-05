#!/bin/bash

# Build Genesis
cd genesis
rm -f genesis/apps/*.cpp
cp src/*.cpp genesis/apps
cd genesis
make -j 4
make update -j 4
cd ..

cp genesis/bin/libgenesis.so ./
cp genesis/bin/apps/genesis-read-newick mem/
cp genesis/bin/apps/genesis-traverse time/

cd ..

# Build CompactTree

cd compacttree
g++ -o time/CompactTree-traverse CompactTree-traverse.cpp
g++ -o mem/CompactTree-read-newick CompactTree-read-newick.cpp

cd ..

# Build phylo-rs and phylotree
cd ..
cargo build --release --example phylo-*
cd comparison

for file in ../target/release/examples/phylo-rs*; do
    if [[ -x "$file" ]] then
        cp $file phylo-rs/time
    fi
done
mv phylo-rs/time/phylo-rs-read-newick phylo-rs/mem/

for file in ../target/release/examples/phylotree-*; do
    if [[ -x "$file" ]] then
        cp $file phylotree/time
    fi
done
mv phylotree/time/phylotree-read-newick phylotree/mem/

