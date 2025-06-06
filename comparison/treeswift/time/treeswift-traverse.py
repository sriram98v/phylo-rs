#!/usr/bin/env python

import treeswift
import sys
import time
from random import sample



treefile=sys.argv[1]

start_time = time.time()
tree = treeswift.read_tree_newick(treefile)
x = list([node for node in tree.traverse_postorder()])

print(f"Internal time: {(time.time()-start_time)}\n")
