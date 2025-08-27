#!/usr/bin/env python

from dendropy import TaxonNamespace, Tree
from dendropy.calculate import treecompare
import sys
import time
from random import sample



t1_file=sys.argv[1]
t2_file=sys.argv[2]

num_sub_taxa = 200//2
subsample_list =  [f"Tip{i}" for i in sample(list(range(200)),num_sub_taxa)]

tns = TaxonNamespace()
tree1 = Tree.get(
        path=t1_file,
        schema="newick",
        taxon_namespace=tns)
tree2 = Tree.get(
        path=t2_file,
        schema="newick",
        taxon_namespace=tns)
start_time = time.time()

tree1.encode_bipartitions()
tree2.encode_bipartitions()
x = treecompare.symmetric_difference(tree1, tree2)

print(x)
print(f"Internal time: {(time.time()-start_time)}\n")
