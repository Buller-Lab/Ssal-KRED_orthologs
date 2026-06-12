import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', required=True, help='Input CSV file')
args = parser.parse_args()

CSV_FILE = args.input

COL_DATASET = "dataset"

COL_VOL = "cavity_volume"
COL_NEG = "cavity_frequency_negatively_charged"
COL_ARO = "cavity_frequency_aromatic"

LEFT_GROUP  = "tested_orthologs"
RIGHT_GROUP = "SwissProt_homologs"

COLOR_VOL = "#5BB8E8"
COLOR_NEG = "#885A95"
COLOR_ARO = "#76b689"
LEGEND_GRAY = "#3C3C3B"

LEGEND_LEFT_TEXT  = "Far-distant tested orthologs"
LEGEND_RIGHT_TEXT = "Far-distant SwissProt homologs"

df = pd.read_csv(CSV_FILE)

data_left  = df[df[COL_DATASET] == LEFT_GROUP]
data_right = df[df[COL_DATASET] == RIGHT_GROUP]

n_left  = len(data_left)
n_right = len(data_right)

fig, (ax_vol, ax_neg, ax_aro) = plt.subplots(
    1, 3,
    figsize=(15, 8),
    sharey=False,
    gridspec_kw={'wspace': 0.12}
)

def standard_kde(data):
    return gaussian_kde(np.asarray(data))

def draw_split_violin(ax, data_all, data_sol, pos, color, fmt=".2f"):
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

        y_min = vmin - extension
        y_max = vmax + extension

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

        ha = 'right' if side < 0 else 'left'
        ax.text(xL - 0.03 if side < 0 else xR + 0.03, med,
                f'{med:{fmt}}', ha=ha, va='center',
                fontweight='bold', fontsize=15, color='black', zorder=8)

    draw_half(data_all, -1, hatch=None)
    draw_half(data_sol, +1, hatch='//////')

def set_dynamic_ylim(ax, data_all, data_sol, buffer_factor=0.05):
    all_data = pd.concat([data_all, data_sol])
    data_min = all_data.min()
    data_max = all_data.max()
    data_range = data_max - data_min if data_max > data_min else 1.0
    buffer = data_range * buffer_factor
    ax.set_ylim(data_min - buffer, data_max + buffer)

draw_split_violin(ax_vol, data_left[COL_VOL], data_right[COL_VOL],
                  pos=1, color=COLOR_VOL, fmt=".1f")
set_dynamic_ylim(ax_vol, data_left[COL_VOL], data_right[COL_VOL])
ax_vol.set_ylabel('Cavity Volume (Å³)', color=COLOR_VOL, fontsize=18, fontweight='bold')
ax_vol.tick_params(axis='y', colors=COLOR_VOL, labelsize=18, width=2)
ax_vol.spines['left'].set_color(COLOR_VOL)
ax_vol.spines['left'].set_linewidth(2.5)

draw_split_violin(ax_neg, data_left[COL_NEG], data_right[COL_NEG],
                  pos=1, color=COLOR_NEG, fmt=".3f")
set_dynamic_ylim(ax_neg, data_left[COL_NEG], data_right[COL_NEG])
ax_neg.set_ylabel('Negatively Charged Frequency', color=COLOR_NEG, fontsize=18, fontweight='bold')
ax_neg.tick_params(axis='y', colors=COLOR_NEG, labelsize=18, width=2)
ax_neg.spines['left'].set_color(COLOR_NEG)
ax_neg.spines['left'].set_linewidth(2.5)

draw_split_violin(ax_aro, data_left[COL_ARO], data_right[COL_ARO],
                  pos=1, color=COLOR_ARO, fmt=".3f")
set_dynamic_ylim(ax_aro, data_left[COL_ARO], data_right[COL_ARO])
ax_aro.set_ylabel('Aromatic Frequency', color=COLOR_ARO, fontsize=18, fontweight='bold')
ax_aro.tick_params(axis='y', colors=COLOR_ARO, labelsize=18, width=2)
ax_aro.spines['left'].set_color(COLOR_ARO)
ax_aro.spines['left'].set_linewidth(2.5)

for ax in [ax_vol, ax_neg, ax_aro]:
    ax.set_xlim(0.45, 1.55)
    ax.set_xticks([])
    for s in ['top', 'right', 'bottom']:
        ax.spines[s].set_visible(False)

legend_elements = [
    plt.Rectangle((0,0),1,1,
                  facecolor=LEGEND_GRAY, edgecolor='black',
                  label=f'{LEGEND_LEFT_TEXT} (n = {n_left})'),
    plt.Rectangle((0,0),1,1,
                  facecolor=LEGEND_GRAY, edgecolor='white',
                  hatch='//////',
                  label=f'{LEGEND_RIGHT_TEXT} (n = {n_right})')
]

ax_neg.legend(handles=legend_elements,
              loc='upper center',
              bbox_to_anchor=(0.5, -0.04),
              ncol=2,
              fontsize=18,
              frameon=False)

plt.subplots_adjust(left=0.11, right=0.96, top=0.96, bottom=0.12)

plt.savefig("split_violins_orthologs_vs_swissprot_homologs_cavity_properties.svg",
            format="svg", bbox_inches="tight", dpi=300)
plt.show()