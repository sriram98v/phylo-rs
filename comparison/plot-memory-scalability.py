#!/usr/bin/env python

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib as mpl
import pandas as pd
import seaborn as sns

sns.set(style='white', rc={'figure.figsize':(14,10)})


mpl.rcParams["figure.figsize"] = (10, 8)
mpl.rcParams["figure.dpi"] = 300

file = "./mem-util.csv"
data = dict()
with open(file, 'r') as f:
    lines = [i[:-1].split(",") for i in f.readlines()[1:]]
    entry_len = len(lines[0])
    entries = [i[:entry_len] for i in lines]
    data = {i[0]:[float(j) for j in i[1:]] for i in entries}
df = pd.DataFrame.from_dict(data)


linestyles = ['solid', 'dashdot', 'dashed']
markers = ["x", ".","^","s"]
libs = ["phylo-rs", "phylotree", "CompactTree", "genesis", "gotree", "dendropy", "treeswift", "ape"]
palette = sns.color_palette("colorblind", len(libs))
colors = {i: palette[n] for n,i in enumerate(libs)}

fig, ax = plt.subplots(1,1)
df_cols = list(df.columns[1:])
handles = [mlines.Line2D([],[],color=colors[c], label=c, linestyle=linestyles[0] if c=="phylo-rs" else linestyles[2], linewidth=1, marker=markers[0] if c=="phylo-rs" else markers[1]) for c in colors]
labels = list(colors.keys())
for n,(key,val) in enumerate(data.items()):
    # ax.set_title(f"Read Newick (Memory)")
    if key =="dendropy":
        ax.plot([1000, 2000, 5000, 10000, 20000, 50000], val[:6], label=f"{key}", linestyle=linestyles[2], marker=markers[1], markersize=5, color=colors[key])
    # elif key=="treeswift":
    #     ax.plot([1000, 2000, 5000, 10000, 20000, 50000, 100000], val[:7], label=f"{key}", linestyle=linestyles[2], marker=markers[1], markersize=5, color=colors[key])
    else:
        if key=="phylo-rs":
            ax.plot([1000, 2000, 5000, 10000, 20000, 50000,100000,200000,500000,1000000], val, label=f"{key}", linestyle=linestyles[0], marker=markers[0], markersize=5, color=colors[key])
        else:
            ax.plot([1000, 2000, 5000, 10000, 20000, 50000,100000,200000,500000,1000000], val, label=f"{key}", linestyle=linestyles[2], marker=markers[1], markersize=5, color=colors[key])
    ax.grid(color='gray', linewidth=0.5, linestyle="--")
plt.xlabel("Taxa Size")
plt.ylabel("Memory (Mb)")
plt.xticks([1000, 100000, 200000, 500000, 1000000], ["1K", "100K", "200K", "500K", "1M"])
ax.set_yscale('log')
fig.legend(handles, labels, loc="upper right", ncols=1)
fig.tight_layout()
plt.subplots_adjust(right=0.83)


# plt.show()
plt.savefig("memory-scalability.png")