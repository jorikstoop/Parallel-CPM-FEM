import pandas as pd
import matplotlib.pyplot as plt
import os
import glob
import re

# define sort key to sort the csv files by iteration and size int values
def sort_key(path):
    name = path.split("/")[-1]  # get filename (last part of path)
    match = re.search(r"t(\d+)_s(\d+)", name) #get value after "t" and "s"
    iteration = int(match.group(1)) # save as int
    domainsize = int(match.group(2))
    return (iteration, domainsize)

# ===============================
# Specify what to profile:
# whole --> shows simloop breakdown as four parts: trac force, fem, cpm, output
# fem_split --> splits the FEM into PCG time and non-PCG time
MODE = "fem_split"

# Specify input and output file names:
# read in all timing csv files from the output folder 
# they are saved with naming convention timing_t<iteration>_s<domainsize>.csv

# sort input files based on iteration and the file size 

INPUT_FILES = sorted(glob.glob("../output/timing_*.csv"), key=sort_key)

OUTPUT_FIG1 = "runtime_breakdown.png"
OUTPUT_FIG2 = "runtime_stacked.png"
# ===============================

def compute_breakdown(csv_file):
    # read csv file from output dir
    df = pd.read_csv(csv_file)

    # create dictionary
    data = dict(zip(df["Process"], df["Time"]))

    # Timing breakdown of simloop
    if MODE == "whole":
        breakdown = {
            "Traction Force": data["tf"],
            "FEM": data["fem"],
            "CPM": data["cpm"],
            "MPI Cell Comm": data["cell_comm"],
            "Output": data["out"]
        }
    elif MODE == "fem_split":
        fem_other = data["fem"] - data["pcg"]
        pcg_nocomm = data["pcg"] - data["node_comm"]
        breakdown = {
            "Traction Force": data["tf"],
            "FEM (non-PCG)": fem_other,
            "PCG (excluding comm)": data["pcg"],
            "MPI Node Comm": data["node_comm"],
            "CPM": data["cpm"],
            "MPI Cell Comm": data["cell_comm"],
            "Output": data["out"]
        }

    # Convert to percentages
    total = sum(breakdown.values())
    percent = {k: 100*v/total for k, v in breakdown.items()}
    return percent, breakdown, total

# ===============================
# Plot breakdown for each file
# ===============================

percent_results = []
raw_results = []
totals = []
labels = []

for file in INPUT_FILES:
    percent, raw, total = compute_breakdown(file)

    percent_results.append(percent)
    raw_results.append(raw)
    totals.append(total)

    # label from filename
    name = os.path.basename(file)
    match = re.search(r"t(\d+)_s(\d+)", name)
    iteration = match.group(1)
    domainsize = match.group(2)
    labels.append(f"Iter {iteration}, Size {domainsize}")

# 1. Plot runtime breakdown as pct out of 100%
fig, ax = plt.subplots()
categories = list(percent_results[0].keys())
colors = plt.get_cmap("Set3")
for i, percent in enumerate(percent_results):
    left = 0
    for j, cat in enumerate(categories):
        value = percent[cat]
        ax.barh(labels[i], value, left=left, color=colors(j),
                label=cat if i == 0 else "")
        left += value
ax.set_xlim(0, 100)
ax.set_xlabel("Percent of SimLoop Runtime")
ax.set_title("Runtime Breakdown")
ax.legend(ncol = 5, loc="lower center", bbox_to_anchor=(0.5, -0.5))
plt.tight_layout()
plt.savefig(OUTPUT_FIG1, dpi=300, bbox_inches="tight")
plt.close()

# 2. Plot raw runtime values
fig2, ax2 = plt.subplots()
categories = list(raw_results[0].keys())
colors = plt.get_cmap("Set3")
x = range(len(labels))
for j, cat in enumerate(categories):
    bottoms = [sum(r[c] for c in categories[:j]) for r in raw_results]
    values = [r[cat] for r in raw_results]
    ax2.bar(x, values, bottom=bottoms, color=colors(j),
            label=cat)
ax2.set_xticks(x)
ax2.set_xticklabels(labels)
#ax2.set_xticklabels(labels, rotation=45, ha="right") #rotated labels
ax2.set_ylabel("Runtime (seconds)")
ax2.set_title("SimLoop Runtime")
ax2.legend(ncol=5, loc="lower center", bbox_to_anchor=(0.5, -0.5))
plt.tight_layout()
plt.savefig(OUTPUT_FIG2, dpi=300, bbox_inches="tight")
plt.close()