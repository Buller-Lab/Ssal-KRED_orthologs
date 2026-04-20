#!/usr/bin/env python3
import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer

parser = argparse.ArgumentParser(description="PCA on selected features colored by target column")
parser.add_argument("csv_file", type=str, help="Path to input CSV")
parser.add_argument("--features", nargs="+", required=True,
                    help="List of feature column names (space-separated)")
parser.add_argument("--color-col", "-c", type=str, required=True,
                    help="Column to color points by (NaN → faded grey)")
parser.add_argument("--random-state", type=int, default=42)

args = parser.parse_args()

OUT_FOLDER = "PCA_top3_features"
os.makedirs(OUT_FOLDER, exist_ok=True)

df = pd.read_csv(args.csv_file)

missing_feats = [f for f in args.features if f not in df.columns]
if missing_feats:
    raise ValueError(f"Features not found in CSV: {missing_feats}")

X_raw = df[args.features].select_dtypes(include=[np.number]).copy()

imputer = SimpleImputer(strategy='median')
X_imputed = pd.DataFrame(imputer.fit_transform(X_raw), columns=X_raw.columns, index=X_raw.index)

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_imputed)

pca = PCA(n_components=2, random_state=args.random_state)
X_pca = pca.fit_transform(X_scaled)

evr = pca.explained_variance_ratio_
cum_evr = np.cumsum(evr)

n_samples = X_scaled.shape[0]
print("\n" + "═" * 60)
print(f"PCA on {len(args.features)} features → colored by: {args.color_col}")
print(f"Number of samples: {n_samples}")
print(f"Explained variance:")
print(f"  PC1: {evr[0]:.3f}  ({evr[0]*100:.1f} %)")
print(f"  PC2: {evr[1]:.3f}  ({evr[1]*100:.1f} %)")
print(f"  Cumulative (PC1+PC2): {cum_evr[1]:.3f}  ({cum_evr[1]*100:.1f} %)")
print("═" * 60)

loadings = pd.DataFrame(
    pca.components_.T,
    columns=['PC1', 'PC2'],
    index=args.features
)

print("\nFeature loadings:")
print(loadings.round(3))

print("\nPC1 sorted by |loading| (strongest contributors):")
pc1_sorted = loadings['PC1'].abs().sort_values(ascending=False)
for feat, abs_load in pc1_sorted.items():
    sign = "+" if loadings.loc[feat, 'PC1'] >= 0 else "−"
    print(f"  {feat:35} {sign}{abs_load:.3f}")

print("\nQuick interpretation:")
if cum_evr[1] >= 0.80:
    strength = "very strong"
elif cum_evr[1] >= 0.65:
    strength = "good"
else:
    strength = "moderate"
print(f"→ The top-3 features capture {strength} amount of variance ({cum_evr[1]*100:.1f} % in 2D).")
if evr[0] > 0.50:
    print("  Separation mainly occurs along PC1.")
print("  Strong gradient along PC1 suggests these descriptors carry most of the signal for the colored property.")
print("═" * 60 + "\n")

color_values = df[args.color_col]
mask_missing = color_values.isna()
mask_valid   = ~mask_missing

valid_values = color_values[mask_valid]
if len(valid_values) > 0:
    vmin, vmax = valid_values.min(), valid_values.max()
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
else:
    norm = mcolors.Normalize(vmin=0, vmax=1)

cmap = mcolors.LinearSegmentedColormap.from_list(
    "purple_white_green",
    ["#885A95", "#ffffff", "#76b689"]
)

fig, ax = plt.subplots(figsize=(7, 7), dpi=100)

ax.scatter(
    X_pca[mask_missing, 0], X_pca[mask_missing, 1],
    c='lightgrey', s=60, edgecolor='none', alpha=0.5, zorder=1
)

sc = ax.scatter(
    X_pca[mask_valid, 0], X_pca[mask_valid, 1],
    c=valid_values, cmap=cmap, norm=norm,
    s=80, edgecolor='black', linewidth=0.6, alpha=1.0, zorder=2
)

ax.set_xlabel(f"PC1 ({evr[0]*100:.1f} % variance)", fontsize=14)
ax.set_ylabel(f"PC2 ({evr[1]*100:.1f} % variance)", fontsize=14)

LIM = 5.0
ax.set_xlim(-LIM, LIM)
ax.set_ylim(-LIM, LIM)

ax.grid(True, linestyle='--', alpha=0.3)
ax.set_aspect('equal')

cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label(f"{args.color_col}", fontsize=14)

plt.tight_layout()

safe_col = args.color_col.replace(" ", "_").replace("%", "pct").replace("(", "").replace(")", "")
fname_base = f"pca_top3_{safe_col}"
plt.savefig(os.path.join(OUT_FOLDER, f"{fname_base}.png"), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(OUT_FOLDER, f"{fname_base}.svg"), bbox_inches='tight')
plt.close(fig)

print(f"Plot saved to: {OUT_FOLDER}/{fname_base}.png  (and .svg)")