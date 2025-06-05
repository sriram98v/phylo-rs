#!/usr/bin/env python

import treeswift
import sys

treefile=sys.argv[1]

tree = treeswift.read_tree_newick(treefile)
print(tree)