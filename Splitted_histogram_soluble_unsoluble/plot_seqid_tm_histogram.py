import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde

CSV_FILE = "histogram_data_solubility_tm_seqid.csv"

df = pd.read_csv(CSV_FILE)
df['Solubility'] = df['Solubility'].astype(bool)
df['seq_id'] = df['sequence_identity'] * 100

all_data = df
soluble_only = df[df['Solubility']]

n_all = len(all_data)
n_sol = len(soluble_only)

color_seq = "#5BB8E8"
color_tm = "#885A95"
legend_gray = "#3C3C3B"

fig, ax1 = plt.subplots(figsize=(9.5, 8))
ax2 = ax1.twinx()

def standard_kde(data):
    data = np.asarray(data)
    return gaussian_kde(data)

def draw_split_violin(ax, data_all, data_sol, pos, color, is_percent=False):
    width = 0.38
    box_width = 0.08

    def draw_half(data, side, hatch=None):
        if len(data) == 0:
            return
        if len(data) == 1:
            y = data.iloc[0]
            ax.hlines(y, pos - box_width/2, pos + box_width/2, color='black', lw=3)
            return

        vmin, vmax = data.min(), data.max()
        kde = standard_kde(data)
        bw = kde.factor * np.std(data)
        extension = 2.5 * bw
        y_min = max(0 if is_percent else -np.inf, vmin - extension)
        y_max = min(100 if is_percent else np.inf, vmax + extension)
        y_grid = np.linspace(y_min, y_max, 1000)
        density = kde(y_grid)
        density = density / density.max() * width

        x = pos + side * density
        xs = np.concatenate([[pos], x, [pos]])
        ys = np.concatenate([y_grid[[0]], y_grid, y_grid[[-1]]])

        ax.fill(xs, ys, color=color, alpha=0.7, hatch=hatch,
                edgecolor='white' if hatch else None, lw=0, zorder=2)
        ax.plot(x, y_grid, color='black', lw=2, zorder=4)

        q1, q3 = np.percentile(data, [25, 75])
        med = np.median(data)
        left = pos + side * box_width
        xL = min(left, pos)
        xR = max(left, pos)

        ax.add_patch(plt.Rectangle((xL, q1), xR-xL, q3-q1,
                     facecolor='white', edgecolor='black', lw=2.5, zorder=6))
        ax.plot([xL, xR], [med, med], color='black', lw=3, zorder=7)

        fmt = ".1f" if is_percent else ".3f"
        ha = 'right' if side < 0 else 'left'
        ax.text(xL - 0.03 if side < 0 else xR + 0.03, med,
                f'{med:{fmt}}', ha=ha, va='center',
                fontweight='bold', fontsize=11, color='black', zorder=8)

    draw_half(data_all, -1, hatch=None)
    draw_half(data_sol, +1, hatch='//////')

ax1.set_ylim(-2, 102)
ax1.set_ylabel('Sequence Identity (%)', color=color_seq, fontsize=14, fontweight='bold')
ax1.tick_params(axis='y', colors=color_seq, labelsize=12, width=2)
ax1.spines['left'].set_color(color_seq)
ax1.spines['left'].set_linewidth(2.5)

draw_split_violin(ax1, all_data['seq_id'], soluble_only['seq_id'],
                  pos=0.85, color=color_seq, is_percent=True)

ax2.set_ylim(-0.02, 1.02)
ax2.set_ylabel('TM-Score', color=color_tm, fontsize=14, fontweight='bold')
ax2.tick_params(axis='y', colors=color_tm, labelsize=12, width=2)
ax2.spines['right'].set_color(color_tm)
ax2.spines['right'].set_linewidth(2.5)
ax2.spines['right'].set_visible(True)

draw_split_violin(ax2, all_data['tm_score'], soluble_only['tm_score'],
                  pos=1.15, color=color_tm, is_percent=False)

ax1.set_xlim(0.45, 1.55)
ax1.set_xticks([])

for ax in (ax1, ax2):
    for s in ['top', 'bottom', 'left' if ax == ax2 else 'right']:
        ax.spines[s].set_visible(False)

legend_elements = [
    plt.Rectangle((0,0),1,1, facecolor=legend_gray, edgecolor='black',
                  label=f'All orthologs (n = {n_all})'),
    plt.Rectangle((0,0),1,1, facecolor=legend_gray, edgecolor='white',
                  hatch='//////', label=f'Soluble only (n = {n_sol})')
]

ax1.legend(handles=legend_elements,
           loc='upper center',
           bbox_to_anchor=(0.5, -0.02),
           ncol=2,
           fontsize=13,
           frameon=False)

plt.subplots_adjust(left=0.16, right=0.84, top=0.96, bottom=0.09)
plt.savefig("split_violin_solubility.pdf", format="pdf", bbox_inches="tight", transparent=False)
plt.show()