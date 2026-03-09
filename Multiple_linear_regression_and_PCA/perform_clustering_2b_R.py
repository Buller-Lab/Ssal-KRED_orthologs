#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.impute import SimpleImputer

parser = argparse.ArgumentParser(description="PCA vs t-SNE colored by selected column – grey & faded for missing")
parser.add_argument("csv_file", type=str, help="Path to the input CSV file")
parser.add_argument("--color-col", "-c", type=str, default="ee 2b (%)", help="Column to color by (NaN → light grey, faded)")
parser.add_argument("--perplexity", type=int, default=10)
parser.add_argument("--learning-rate", type=str, default="auto")
parser.add_argument("--random-state", type=int, default=42)

args = parser.parse_args()

FEATURES = [
    'tm_score', 'cavity_frequency_Aromatic','cavity_frequency_Polar_uncharged'
]

df = pd.read_csv(args.csv_file)
available_features = [f for f in FEATURES if f in df.columns]

if not available_features:
    raise ValueError("No known numeric feature columns found in the CSV")

X = df[available_features].select_dtypes(include=[np.number]).copy()

imputer = SimpleImputer(strategy='median')
X_imputed = pd.DataFrame(imputer.fit_transform(X), columns=X.columns, index=X.index)

if 'e_value' in X_imputed.columns:
    X_imputed['e_value'] = np.log1p(X_imputed['e_value'].clip(lower=0))

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_imputed)

pca = PCA(n_components=2, random_state=args.random_state)
X_pca = pca.fit_transform(X_scaled)

tsne = TSNE(
    n_components=2,
    perplexity=args.perplexity,
    learning_rate=args.learning_rate,
    max_iter=1000,
    random_state=args.random_state,
    init='pca'
)
X_tsne = tsne.fit_transform(X_scaled)

color_col = args.color_col

if color_col not in df.columns:
    raise ValueError(f"Column '{color_col}' not found. Available: {list(df.columns)}")

values = df[color_col]

mask_missing = values.isna()
mask_valid   = ~mask_missing

valid_values = values[mask_valid]
if len(valid_values) > 0:
    vmin = valid_values.min()
    vmax = valid_values.max()
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
else:
    norm = mcolors.Normalize(vmin=0, vmax=1)

from matplotlib.colors import LinearSegmentedColormap

cmap = LinearSegmentedColormap.from_list(
    "white_blue",
    ["#FFFFFF", "#2A5BA4"]
)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6.5))

ax1.scatter(
    X_pca[mask_missing, 0], X_pca[mask_missing, 1],
    c='lightgrey',         
    s=70,
    edgecolor='none',       
    linewidth=0.4,
    alpha=0.20,            
    zorder=1
)

ax1.scatter(
    X_pca[mask_valid, 0], X_pca[mask_valid, 1],
    c=valid_values,
    cmap=cmap,
    norm=norm,
    s=70,
    edgecolor='black',
    linewidth=0.6,
    alpha=1.0,
    zorder=2              
)

ax1.set_title(f"PCA\nExplained variance: {sum(pca.explained_variance_ratio_):.1%}")
ax1.set_xlabel("PC1")
ax1.set_ylabel("PC2")
ax1.grid(True, alpha=0.3)

ax2.scatter(
    X_tsne[mask_missing, 0], X_tsne[mask_missing, 1],
    c='lightgrey',
    s=70,
    edgecolor='none',
    linewidth=0.4,
    alpha=0.20,
    zorder=1
)

ax2.scatter(
    X_tsne[mask_valid, 0], X_tsne[mask_valid, 1],
    c=valid_values,
    cmap=cmap,
    norm=norm,
    s=70,
    edgecolor='black',
    linewidth=0.6,
    alpha=1.0,
    zorder=2
)

ax2.set_title(f"t-SNE\nperplexity = {args.perplexity}")
ax2.set_xlabel("t-SNE 1")
ax2.set_ylabel("t-SNE 2")
ax2.grid(True, alpha=0.3)

cbar_ax = fig.add_axes([0.88, 0.15, 0.015, 0.7])
cbar = fig.colorbar(
    plt.cm.ScalarMappable(norm=norm, cmap=cmap),
    cax=cbar_ax, orientation='vertical'
)
cbar.set_label(f"{color_col} value", fontsize=11)

plt.subplots_adjust(right=0.82)

plt.savefig("clustering.png", dpi=300, bbox_inches='tight')
plt.savefig("clustering.svg", bbox_inches='tight')
plt.close(fig)