#!/usr/bin/env python
import treeswift
import sys
import time
from random import sample



treefile=sys.argv[1]

num_sub_taxa = 200//2
subsample_list =  [f"Tip{i}" for i in sample(list(range(200)),num_sub_taxa)]

tree = treeswift.read_tree_newick(treefile)
start_time = time.time()
x = tree.mrca(subsample_list)

print(f"Internal time: {(time.time()-start_time)}\n")
