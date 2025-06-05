#!/usr/bin/env python

from dendropy import TaxonNamespace, Tree
import sys
import time
from random import sample



treefile=sys.argv[1]

num_sub_taxa = 200//2
subsample_list =  [f"Tip{i}" for i in sample(list(range(200)),num_sub_taxa)]

start_time = time.time()
tns = TaxonNamespace()
tree = Tree.get(path=treefile, schema="newick", taxon_namespace=tns)
x = tree.extract_tree_with_taxa_labels(subsample_list)

print(f"Internal time: {(time.time()-start_time)}\n")
