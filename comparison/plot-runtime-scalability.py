#!/usr/bin/env python

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib as mpl
import pandas as pd
import seaborn as sns
import string

sns.set(style='white', rc={'figure.figsize':(14,10)})


mpl.rcParams["figure.figsize"] = (10, 8)
mpl.rcParams["figure.dpi"] = 300

file = "./runtimes.csv"
data = dict()
with open(file, 'r') as f:
    lines = [i.split(",") for i in f.readlines()]
    entry_len = len(lines[0])
    entries = [i[:entry_len] for i in lines]
    data = {i[0]:[float(j) for j in i[1:]] for i in entries}
df = pd.DataFrame.from_dict(data)

methods = {"traverse":"Tree Traversal",
 "lca":"Least Common Ancestor Retrieval",
 "nni": "Nearest Neighbor Interchange",
 "yts":"Yule Tree Simulation",
 "contract": "Tree Contraction",
 "rfs":"Robinson Foulds metric computation"}
m_keys = list(methods.keys())
linestyles = ['solid', 'dashdot', 'dashed']
markers = ["x", ".","^","s"]
libs = ["phylo-rs", "phylotree", "CompactTree", "genesis", "gotree", "dendropy", "treeswift", "ape"]
palette = sns.color_palette("colorblind", len(libs))
colors = {i: palette[n] for n,i in enumerate(libs)}

fig, ax = plt.subplots(3, 2, sharex=True, sharey=True, figsize=(15, 10))
df_cols = list(df.columns[1:])
handles = [mlines.Line2D([],[],color=colors[c], label=c, linestyle=linestyles[0] if c=="phylo-rs" else linestyles[2], linewidth=1, marker=markers[0] if c=="phylo-rs" else markers[1]) for c in colors]
labels = list(colors.keys())
for m_idx in range(len(methods)):
    subplot = ax[m_idx%3,m_idx//3]
    method = methods[m_keys[m_idx]]
    m_cols = [i for i in df_cols if m_keys[m_idx] in i]
    sub_df = df[['algorithms', *m_cols]]
    for line in m_cols:
        l_key = "-".join(line.split("-")[:-1])
        val = sub_df[line].values * 1000
        if l_key=="phylo-rs":
            subplot.plot(range(200, 10001, 200), val, label=f"{l_key}", linestyle=linestyles[0], marker=markers[0], markevery=2, color=colors[l_key])
        else:
            subplot.plot(range(200, 10001, 200), val, label=f"{l_key}", linestyle=linestyles[2], marker=markers[1], markevery=2, color=colors[l_key])
    subplot.set_title(f"({string.ascii_uppercase[m_idx]}) {method}")

    subplot.set_yscale('log')
    subplot.grid(color='gray', linewidth=0.5, linestyle="--")
fig.legend(handles, labels, loc="center right", ncols=1)
fig.supxlabel("Taxa Size")
fig.supylabel("Time (ms)")
fig.suptitle("Runtime scalability analysis")
fig.tight_layout()
plt.subplots_adjust(right=0.88)
# plt.show()
plt.savefig("runtime-scalability.png")