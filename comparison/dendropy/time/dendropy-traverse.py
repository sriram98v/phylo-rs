#!/usr/bin/env python

from dendropy import TaxonNamespace, Tree
import sys
import time


treefile=sys.argv[1]

tns = TaxonNamespace()
tree = Tree.get(path=treefile, schema="newick", taxon_namespace=tns)

start_time = time.time()
x = list(tree.postorder_node_iter())
print(f"Internal time: {(time.time()-start_time)}\n")
