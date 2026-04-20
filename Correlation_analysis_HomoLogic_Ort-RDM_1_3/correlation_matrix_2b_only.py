import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import scipy.stats
from matplotlib.colors import LinearSegmentedColormap

if len(sys.argv) != 2:
    print("Usage: python script.py input.csv")
    sys.exit(1)

csv_file = sys.argv[1]
base_name = os.path.splitext(os.path.basename(csv_file))[0]
output_svg = f"{base_name}_correlation_2b.svg"
output_png = f"{base_name}_correlation_2b.png"

FEATURES = [
    "n_bindingsite_residues", "cavity_volume", "cavity_area", "cavity_avg_depth", "cavity_avg_hydropathy",
    "cavity_frequency_Alipathic_apolar", "cavity_frequency_Aromatic", "cavity_frequency_Polar_uncharged",
    "cavity_frequency_Negatively_charged", "cavity_frequency_Positively_charged", "cavity_n_points",
    "sequence_identity", "rmsd", "tm_score", "rmsd_bindingsite", "bindingsite_similarity_score"
]

TARGETS = [ "2b (R)", "2b (S)", "ee 2b (%)"]

df = pd.read_csv(csv_file)

missing_markers = ['n.d.', 'nd', 'N.D.', 'NA', 'N/A', 'n/a', '-', 'none', 'None', '']

columns_of_interest = [c for c in FEATURES + TARGETS if c in df.columns]
df[columns_of_interest] = df[columns_of_interest].replace(missing_markers, np.nan)

for col in columns_of_interest:
    df[col] = pd.to_numeric(df[col], errors='coerce')

X = df[columns_of_interest].select_dtypes(include=np.number)

corr = X.corr()
sub = corr.loc[[f for f in FEATURES if f in corr.index],
               [t for t in TARGETS if t in corr.columns]]

plot_matrix = sub.copy()

target_sample_sizes = df[TARGETS].notna().sum().to_dict()

p_matrix = pd.DataFrame(index=plot_matrix.index, columns=plot_matrix.columns, dtype=float)

for fi, f in enumerate(plot_matrix.index):
    for ti, t in enumerate(plot_matrix.columns):
        x = df[f].values
        y = df[t].values
        mask = ~np.isnan(x) & ~np.isnan(y)
        if mask.sum() >= 3:
            _, p = scipy.stats.pearsonr(x[mask], y[mask])
            p_matrix.iloc[fi, ti] = p
        else:
            p_matrix.iloc[fi, ti] = np.nan

def format_annot(r, p):
    if np.isnan(r) or np.isnan(p):
        return ""
    r_str = f"{r:.2f}"
    if p < 0.001:
        sig = "***"
    elif p < 0.01:
        sig = "**"
    elif p < 0.05:
        sig = "*"
    else:
        sig = ""
    if sig:
        return f"{r_str}\n({sig})"
    else:
        return f"{r_str}\n "

annot_matrix = plot_matrix.copy().astype(str)
for i in range(len(plot_matrix.index)):
    for j in range(len(plot_matrix.columns)):
        annot_matrix.iloc[i,j] = format_annot(plot_matrix.iloc[i,j], p_matrix.iloc[i,j])

div_cmap = LinearSegmentedColormap.from_list(
    "custom_div", ["#885A95", "#ffffff", "#76b689"]
)

n_features = len(plot_matrix.index)
n_targets  = len(plot_matrix.columns)

cell_height = 1.00
cell_width  = 2.20

fig_w = 7 + n_targets * cell_width
fig_h = 2.6 + n_features * cell_height + 1.8

plt.figure(figsize=(fig_w, fig_h))

ax = sns.heatmap(
    plot_matrix,
    cmap=div_cmap,
    vmin=-1, vmax=1,
    center=0,
    annot=annot_matrix,
    fmt="",
    annot_kws={
        "size": 20,
        "weight": "medium",
        "linespacing": 0.95,
        "verticalalignment": "top",
    },
    linewidths=1.0,
    linecolor="white",
    cbar_kws={
        'label': 'Pearson r',
        'shrink': 1.0,
        'aspect': 30,
        'pad': 0.03
    },
    xticklabels=True,
    yticklabels=True,
    square=False,
)

for text in ax.texts:
    content = text.get_text()
    if content.startswith('(') and content.endswith(')'):
        text.set_fontsize(16)
        text.set_fontweight('normal')

for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_edgecolor('black')
    spine.set_linewidth(1.2)

ax.tick_params(axis='y', rotation=0, labelsize=20)
ax.tick_params(axis='x', rotation=0, labelsize=18, pad=8)

ax.set_xticklabels(plot_matrix.columns, ha='center', fontsize=18, fontweight='bold')

ax.set_ylabel("Features", fontsize=24, fontweight='bold', labelpad=16)

for j, t in enumerate(plot_matrix.columns):
    n = target_sample_sizes.get(t, 0)
    x_pos = j + 0.5
    y_pos_n = -0.6

    ax.text(
        x_pos, y_pos_n,
        f"(N={n})",
        ha='center', va='top',
        fontsize=16,
        fontweight='normal',
        transform=ax.transData
    )

cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=20)
cbar.set_label('Pearson r', fontsize=24, fontweight='bold', labelpad=12)

plt.tight_layout(rect=[0.02, 0.22, 0.94, 0.98])

plt.savefig(output_svg, format="svg", bbox_inches="tight", dpi=400)
plt.savefig(output_png, format="png", bbox_inches="tight", dpi=400)
print(f"Saved to: {output_svg}")

plt.show()