
import treeswift
from random import sample
from tqdm import tqdm
import time
import tracemalloc

num_iter = 100

with open("sim_trees", "r") as f:
    trees = [i.split("\t") for i in f.readlines()]

taxa_sizes = [len(tree) for tree in trees]

# print("tree sizes")
# with open("treeswift-sizes.csv", "w") as f:
#     tracemalloc.start()
#     for (t1_str, t2_str) in tqdm(trees):
#         all_diffs = 0
#         num_diffs = 0
#         for i in range(num_iter):
#             snapshot1 = tracemalloc.take_snapshot()
#             tree1 = treeswift.read_tree_newick(t1_str)
#             snapshot2 = tracemalloc.take_snapshot()
#             top_stats = snapshot2.compare_to(snapshot1, 'traceback')

#             if top_stats[0].size_diff>0:
#                 all_diffs += top_stats[0].size_diff
#                 num_diffs+=1
#         print(f"{len(list(tree1.labels()))},{(all_diffs/num_diffs)/1000}\n")
#         f.write(f"{len(list(tree1.labels()))},{(all_diffs/num_diffs)/1000}\n")
#     tracemalloc.stop()

print("Contract")
with open("treeswift-contract-times.csv", "w") as f:
    for (t1_str, t2_str) in tqdm(trees):
        num_sub_taxa = 200//2
        subsample_list =  [[str(i) for i in sample(list(range(200)),num_sub_taxa)] for _ in range(num_iter)]
        start_time = time.time()
        for i in range(num_iter):
            tree1 = treeswift.read_tree_newick(t1_str)

            x = tree1.extract_tree_with(subsample_list[i])
        f.write(f"{len(list(tree1.labels()))},{((time.time()-start_time)*1000)/num_iter}\n")

print("Postord")
with open("treeswift-postord-times.csv", "w") as f:
    for (t1_str, t2_str) in tqdm(trees):
        start_time = time.time()
        for i in range(num_iter):
            tree1 = treeswift.read_tree_newick(t1_str)

            x = list([node for node in tree1.traverse_postorder()])
        f.write(f"{len(list(tree1.labels()))},{((time.time()-start_time)*1000)/num_iter}\n")

print("MRCA")
with open("treeswift-mrca-times.csv", "w") as f:
    for (t1_str, t2_str) in tqdm(trees):
        subsample_list =  [[str(i) for i in sample(list(range(200)),2)] for _ in range(num_iter)]
        start_time = time.time()
        for i in range(num_iter):
            tree1 = treeswift.read_tree_newick(t1_str)

            x = tree1.mrca(subsample_list[i])
        print(f"{len(list(tree1.labels()))},{((time.time()-start_time)*1000)/num_iter}")
        f.write(f"{len(list(tree1.labels()))},{((time.time()-start_time)*1000)/num_iter}\n")