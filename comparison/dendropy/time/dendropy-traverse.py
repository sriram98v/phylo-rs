#!/usr/bin/env python

from dendropy import TaxonNamespace, Tree
import sys
import time


treefile=sys.argv[1]

start_time = time.time()
tns = TaxonNamespace()
tree = Tree.get(path=treefile, schema="newick", taxon_namespace=tns)
x = list(tree.postorder_node_iter())
print(f"Internal time: {(time.time()-start_time)}\n")
