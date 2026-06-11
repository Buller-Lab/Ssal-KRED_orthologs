import argparse
import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import LeaveOneOut
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import os
import subprocess
import tempfile
from pathlib import Path

FIXED_FEATURES = [
    "cavity_frequency_Polar_uncharged",
    "cavity_frequency_Aromatic",
    "tm_score"
]

TARGET = "2b (R)"
SEQ_COL = "sequence"

def run_mmseqs_clustering(sequences, threshold=0.9, tmp_dir=None):
    if tmp_dir is None:
        tmp_dir = tempfile.mkdtemp(prefix="mmseqs_")
    fasta_path = Path(tmp_dir) / "sequences.fasta"
    cluster_out = Path(tmp_dir) / "cluster"
    with open(fasta_path, "w") as f:
        for i, seq in enumerate(sequences):
            f.write(f">seq_{i}\n{seq}\n")
    cmd = [
        "mmseqs", "easy-cluster",
        str(fasta_path), str(cluster_out), str(Path(tmp_dir) / "tmp"),
        "--min-seq-id", str(threshold),
        "--cov-mode", "0",
        "--cluster-mode", "0",
        "--threads", "4"
    ]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    rep_path = cluster_out.with_name(cluster_out.name + "_rep_seq.fasta")
    centroids = []
    with open(rep_path) as f:
        for line in f:
            if line.startswith(">"):
                idx = int(line[1:].split("_")[1])
                centroids.append(idx)
    return sorted(centroids)

parser = argparse.ArgumentParser(description="Fixed features MLR with mmseqs2 clustering")
parser.add_argument("-i", "--input", required=True, help="Input CSV file")
args = parser.parse_args()

safe_target = TARGET.replace(" ", "_").replace("(", "").replace(")", "").replace("%", "pct").replace("/", "_")
output_folder = "generalization_capacity"
os.makedirs(output_folder, exist_ok=True)

df = pd.read_csv(args.input)
df = df[FIXED_FEATURES + [TARGET, SEQ_COL]].dropna().reset_index(drop=True)
print(f"Original samples: {len(df)}")

thresholds = [0.90, 0.80, 0.70, 0.60, 0.50]
results = []

for thresh in thresholds:
    sequences = df[SEQ_COL].tolist()
    centroid_indices = run_mmseqs_clustering(sequences, threshold=thresh)
    df_clust = df.iloc[centroid_indices].reset_index(drop=True)
    n_rem = len(df_clust)
    print(f"Clustering at {int(thresh*100)}%: {n_rem} representatives")
    if n_rem < 10:
        continue
    X_raw = df_clust[FIXED_FEATURES]
    y = df_clust[TARGET].values
    scaler = StandardScaler()
    X_std = scaler.fit_transform(X_raw)
    model_full = LinearRegression().fit(X_std, y)
    y_pred_full = model_full.predict(X_std)
    r2_full = r2_score(y, y_pred_full)
    rmse_full = np.sqrt(mean_squared_error(y, y_pred_full))
    cv = LeaveOneOut()
    y_pred_cv = np.zeros_like(y, dtype=float)
    ols = LinearRegression()
    for train_idx, test_idx in cv.split(X_std):
        ols.fit(X_std[train_idx], y[train_idx])
        y_pred_cv[test_idx] = ols.predict(X_std[test_idx])
    press = np.sum((y - y_pred_cv)**2)
    tss = np.sum((y - np.mean(y))**2)
    q2 = 1 - press / tss if tss > 1e-10 else np.nan
    rmse_loo = np.sqrt(mean_squared_error(y, y_pred_cv))
    mae_loo = np.mean(np.abs(y - y_pred_cv))
    pcc_loo, _ = pearsonr(y, y_pred_cv)
    delta = r2_full - q2 if not np.isnan(q2) else np.nan
    results.append({
        'threshold': f'{int(thresh*100)}%',
        'n_remaining': n_rem,
        'r2_full': round(r2_full, 3),
        'adj_r2': round(1 - (1 - r2_full) * (n_rem - 1) / (n_rem - 4), 3) if n_rem > 4 else round(r2_full, 3),
        'rmse_full': round(rmse_full, 3),
        'q2_loo': round(q2, 3) if not np.isnan(q2) else None,
        'rmse_loo': round(rmse_loo, 3),
        'mae_loo': round(mae_loo, 3),
        'pcc_loo': round(pcc_loo, 3),
        'delta_r2_q2': round(delta, 3) if not np.isnan(delta) else None
    })

summary_df = pd.DataFrame(results)
print("\n" + "═"*110)
print("MMSEQS2 CLUSTERING + FIXED FEATURES PERFORMANCE SUMMARY")
print("═"*110)
print(summary_df.to_string(index=False))
summary_path = os.path.join(output_folder, f"clustering_summary_{safe_target}.csv")
summary_df.to_csv(summary_path, index=False)
print(f"\nSummary table saved to: {summary_path}")

COLORS = ['#1f77b4', '#76b689', '#d62728']
fig, ax1 = plt.subplots(figsize=(9, 6))
x = np.arange(len(summary_df))
width = 0.35
bars = ax1.bar(x - width/2, summary_df['n_remaining'], width, color=COLORS[0], alpha=0.85, label='Remaining Datapoints')
ax1.set_xlabel('Sequence Identity Threshold', fontsize=12)
ax1.set_ylabel('Number of Remaining Datapoints', fontsize=12, color=COLORS[0])
ax1.tick_params(axis='y', labelcolor=COLORS[0])
ax1.set_xticks(x)
ax1.set_xticklabels(summary_df['threshold'])
ax2 = ax1.twinx()
ax2.plot(x, summary_df['delta_r2_q2'], color=COLORS[2], marker='o', linestyle='-', linewidth=2.5, markersize=8, label='R² - Q² Gap')
ax2.set_ylabel('R² - Q² Gap (Overfitting)', fontsize=12, color=COLORS[2])
ax2.tick_params(axis='y', labelcolor=COLORS[2])
lines, labels = ax2.get_legend_handles_labels()
ax1.legend(bars, ['Remaining Datapoints'], loc='upper right')
ax2.legend(lines, ['R² - Q² Gap'], loc='upper left')
plt.title(f'Generalization Capacity – {TARGET}\nEffect of Homology Reduction', fontsize=14, pad=20)
plt.grid(True, linestyle='--', alpha=0.3)
plt.tight_layout()
plot_path = os.path.join(output_folder, f"generalization_barplot_{safe_target}.png")
plt.savefig(plot_path, dpi=300, bbox_inches="tight")
plt.savefig(plot_path.replace(".png", ".svg"), format="svg", bbox_inches="tight")
plt.close()
print(f"Barplot saved to: {plot_path}")