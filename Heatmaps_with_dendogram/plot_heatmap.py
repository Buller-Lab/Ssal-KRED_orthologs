import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator
from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib.colors import LinearSegmentedColormap
import tempfile
import subprocess
import numpy as np
import shutil

if len(sys.argv) != 2:
    print("Usage: python plot_heatmap.py <input.csv>")
    sys.exit(1)

CSV_FILE = sys.argv[1]
base_name = os.path.splitext(os.path.basename(CSV_FILE))[0]
OUTPUT_FIG = f"{base_name}_multi_heatmap_phylogenetic_order.svg"

CLUSTALO_PATH = "clustalo"

if shutil.which(CLUSTALO_PATH) is None:
    print(f"Error: '{CLUSTALO_PATH}' not found in PATH.")
    sys.exit(1)

HEATMAP_COL_1 = "1b"
PAIRED_COLUMNS = {
    "2b": ["2b (R)", "2b (S)"],
    "3b": ["3b (R)", "3b (S)"],
    "4b": ["4b (R)", "4b (S)"],
    "5b": ["5b (R)", "5b (S)"]
}
PAIR_COLORS = ['#FFFFFF', '#3375BE']
PAIR_VMAX = 7.0
CMAP_1 = LinearSegmentedColormap.from_list('green_scale', ['#FFFFFF', '#5EAA74'])
VMAX_1 = 3.0
CMAP_PAIR = LinearSegmentedColormap.from_list('blue_scale', PAIR_COLORS)

CELL_HEIGHT = 0.18
DENDRO_WIDTH = 0.75
LEFT_MARGIN  = 1.6
RIGHT_MARGIN = 0.5
TOP_MARGIN   = 0.7
BOTTOM_MARGIN = 1.0
WSPACE = 0.08

df = pd.read_csv(CSV_FILE)
df[HEATMAP_COL_1] = pd.to_numeric(df[HEATMAP_COL_1], errors="coerce")
all_paired_cols = [col for pair in PAIRED_COLUMNS.values() for col in pair]
for col in all_paired_cols:
    df[col] = pd.to_numeric(df[col], errors="coerce")

required_cols = ["ID", "sequence", HEATMAP_COL_1] + all_paired_cols
df = df.dropna(subset=required_cols).reset_index(drop=True)
print(f"Loaded {len(df)} sequences")

n_rows = len(df)

with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
    fasta_file = f.name
    for _, row in df.iterrows():
        rec = SeqRecord(Seq(row["sequence"]), id=str(row["ID"]), description="")
        SeqIO.write(rec, f, "fasta")

aligned_file = fasta_file + ".aln"
subprocess.run([CLUSTALO_PATH, "-i", fasta_file, "-o", aligned_file,
                "--force", "--auto", "--outfmt=fasta", "--verbose"], check=True)

alignment = MultipleSeqAlignment(list(SeqIO.parse(aligned_file, "fasta")))
calculator = DistanceCalculator("blosum62")
dm = calculator.get_distance(alignment)
ids = list(dm.names)

n = len(ids)
condensed_dist = [dm[ids[i], ids[j]] for i in range(n) for j in range(i+1, n)]
Z = linkage(condensed_dist, method='average')
dendro = dendrogram(Z, labels=ids, no_plot=True)
order = dendro['ivl']

os.unlink(fasta_file)
os.unlink(aligned_file)

df = df.set_index("ID")
df_ordered = df.loc[order].reset_index().copy()

height_total = TOP_MARGIN + BOTTOM_MARGIN + n_rows * CELL_HEIGHT
width_total = LEFT_MARGIN + 0.55 + 4 * 0.70 + DENDRO_WIDTH + RIGHT_MARGIN

fig = plt.figure(figsize=(width_total, height_total))

gs = fig.add_gridspec(1, 6,
                      width_ratios=[0.55, 0.70, 0.70, 0.70, 0.70, DENDRO_WIDTH],
                      wspace=WSPACE,
                      left  = LEFT_MARGIN  / width_total,
                      right = 1 - RIGHT_MARGIN / width_total,
                      top   = 1 - TOP_MARGIN   / height_total,
                      bottom= BOTTOM_MARGIN    / height_total)

ax1     = fig.add_subplot(gs[0, 0])
ax2     = fig.add_subplot(gs[0, 1])
ax3     = fig.add_subplot(gs[0, 2])
ax4     = fig.add_subplot(gs[0, 3])
ax5     = fig.add_subplot(gs[0, 4])
ax_dend = fig.add_subplot(gs[0, 5])

axes_paired = [ax2, ax3, ax4, ax5]
pair_labels = list(PAIRED_COLUMNS.keys())

# ─── Heatmaps ─────────────────────────────────────────────
sns.heatmap(df_ordered[[HEATMAP_COL_1]].values,
            cmap=CMAP_1, vmin=0, vmax=VMAX_1, cbar=False,
            yticklabels=df_ordered["ID"], xticklabels=[HEATMAP_COL_1],
            linewidths=1.2, linecolor="white",
            ax=ax1)

ax1.tick_params(axis='both', labelsize=10)
ax1.tick_params(axis='x', rotation=90)

for spine in ax1.spines.values():
    spine.set_visible(True)
    spine.set_edgecolor('black')
    spine.set_linewidth(1.0)

for ax, label in zip(axes_paired, pair_labels):
    cols = PAIRED_COLUMNS[label]
    data = df_ordered[cols].values
    sns.heatmap(data,
                cmap=CMAP_PAIR, vmin=0, vmax=PAIR_VMAX + 0.05, cbar=False,
                yticklabels=False, xticklabels=cols,
                linewidths=1.2, linecolor="white",
                ax=ax)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_edgecolor('black')
        spine.set_linewidth(1.0)
    ax.tick_params(axis='both', labelsize=10)
    ax.tick_params(axis='x', rotation=90)

dendrogram(Z, labels=ids, ax=ax_dend, orientation='right',
           no_labels=False, color_threshold=0,
           above_threshold_color='black')

for line in ax_dend.get_lines():
    line.set_color('black')
    line.set_linewidth(1.0)

ax_dend.invert_yaxis()
ax_dend.axis('off')

plt.savefig(OUTPUT_FIG, bbox_inches="tight", dpi=400)
print(f"Figure saved: {OUTPUT_FIG}")