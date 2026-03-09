import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
from matplotlib.colors import LinearSegmentedColormap

if len(sys.argv) != 2:
    print("Usage: python script.py input.csv")
    sys.exit(1)

csv_file = sys.argv[1]
base_name = os.path.splitext(os.path.basename(csv_file))[0]
output_svg = f"{base_name}_correlation.svg"

FEATURES = [
    "sequence_identity", "rmsd", "tm_score", "rmsd_bindingsite",
    "n_bindingsite_residues", "n_aligned_residues", "aligned_fraction", "n_res_plddt", "n_res_plddt_bs",
    "median_plddt", "median_plddt_bs", "avg_plddt", "avg_plddt_bs", "std_plddt", "std_plddt_bs",
    "iqr_plddt", "iqr_plddt_bs", "p10_plddt", "p10_plddt_bs", "bindingsite_similarity_score",
    "non_equivalent_residues", "cavity_volume", "cavity_area", "cavity_avg_depth", "cavity_avg_hydropathy",
    "cavity_frequency_Alipathic_apolar", "cavity_frequency_Aromatic", "cavity_frequency_Polar_uncharged",
    "cavity_frequency_Negatively_charged", "cavity_frequency_Positively_charged",
    "cavity_frequency_Non-standard", "cavity_rmsd", "cavity_n_points"
]

TARGETS = ["1b", "2b (R)", "2b (S)", "ee 2b (%)", "3b (R)", "3b (S)", "ee 3b (%)",
           "4b (R)", "4b (S)", "ee 4b (%)", "5b (R)", "5b (S)", "ee 5b (%)", ]

df = pd.read_csv(csv_file)
X = df[FEATURES + TARGETS].select_dtypes(include=np.number)
corr = X.corr()
sub = corr.loc[FEATURES, TARGETS]
mean_abs_corr = sub.abs().mean(axis=1).sort_values(ascending=False)
ordered_features = mean_abs_corr.index.tolist()
plot_matrix = sub.loc[ordered_features]

div_cmap = LinearSegmentedColormap.from_list(
    "custom_div",
    ["#2e5da6", "#ffffff", "#c0392b"]
)

n_features = len(ordered_features)
cell_h = 0.22
fig_h = 1.2 + n_features * cell_h + 1.0
fig_w = 12 + len(TARGETS) * 0.65

plt.figure(figsize=(fig_w, fig_h))

ax = sns.heatmap(
    plot_matrix,
    cmap=div_cmap,
    vmin=-1, vmax=1,
    center=0,
    annot=True,
    fmt=".2f",
    annot_kws={"size": 9, "weight": "medium"},
    linewidths=1.2,
    linecolor="white",
    cbar_kws={'label': 'Pearson r', 'shrink': 0.6, 'aspect': 30},
    xticklabels=True,
    yticklabels=True
)

for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_edgecolor('black')
    spine.set_linewidth(1.2)

ax.tick_params(axis='both', which='major', labelsize=10)
ax.tick_params(axis='x', rotation=55)
ax.set_xticklabels(ax.get_xticklabels(), ha='right')
ax.tick_params(axis='y', rotation=0, labelsize=10)

ax.set_xlabel("Targets", fontsize=12, labelpad=12)
ax.set_ylabel("Features (sorted by mean |r|)", fontsize=12, labelpad=12)
ax.set_title("Feature–Target Pearson Correlations", fontsize=14, pad=18, weight='medium')

plt.tight_layout(rect=[0, 0, 1, 0.98])

plt.savefig(output_svg, format="svg", bbox_inches="tight", dpi=400)
print(f"Saved to: {output_svg}")

plt.show()