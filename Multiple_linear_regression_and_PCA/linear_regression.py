import argparse
import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression, LassoCV
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score, mean_squared_error
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import seaborn as sns
import os

FEATURES = [
    "sequence_identity", "rmsd", "tm_score", "rmsd_bindingsite",
    "n_bindingsite_residues", "n_aligned_residues", "aligned_fraction", "bindingsite_similarity_score",
    "non_equivalent_residues", "cavity_volume", "cavity_area", "cavity_avg_depth", "cavity_avg_hydropathy",
    "cavity_frequency_Alipathic_apolar", "cavity_frequency_Aromatic", "cavity_frequency_Polar_uncharged",
    "cavity_frequency_Negatively_charged", "cavity_frequency_Positively_charged",
    "cavity_frequency_Non-standard", "cavity_rmsd", "cavity_n_points"
]

#TARGETS = ["2b (R)"]
TARGETS = ["ee 2b (%)"]

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", required=True)
args = parser.parse_args()

plot_folder = "mlr_lasso_top3"
os.makedirs(plot_folder, exist_ok=True)

df = pd.read_csv(args.input)
available_features = [f for f in FEATURES if f in df.columns]
df = df[available_features + TARGETS].dropna()

print(f"Final modeling shape: {df.shape}")

X_raw = df[available_features]
y_dict = {t: df[t] for t in TARGETS if t in df.columns}

scaler = StandardScaler()
X_std = pd.DataFrame(
    scaler.fit_transform(X_raw),
    columns=X_raw.columns,
    index=X_raw.index
)

results = []

distinct_colors = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
    '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'
]

for target_name, y in y_dict.items():
    lasso = LassoCV(cv=5, random_state=42, max_iter=30000, n_alphas=200)
    lasso.fit(X_std, y)

    coef_lasso = pd.Series(lasso.coef_, index=X_std.columns)

    top_features = coef_lasso.abs().sort_values(ascending=False).head(3).index.tolist()
    n_selected = len(top_features)

    if n_selected == 0:
        continue

    X_selected = X_std[top_features]
    model = LinearRegression().fit(X_selected, y)

    coef_series = pd.Series(model.coef_, index=top_features, name="β_std")
    coef_series = coef_series.sort_values(key=abs, ascending=False)

    y_pred = model.predict(X_selected)

    r2_val = r2_score(y, y_pred)
    rmse_val = np.sqrt(mean_squared_error(y, y_pred))
    pcc, _ = pearsonr(y, y_pred)
    n = len(y)
    p = n_selected
    adj_r2 = 1 - (1 - r2_val) * (n - 1) / (n - p - 1) if n > p + 1 else r2_val

    sizes = np.abs(coef_series)
    labels = [f"{feat}\n{coef:+.3f}" for feat, coef in coef_series.items()]

    pie_colors = distinct_colors[:n_selected]

    fig, ax = plt.subplots(figsize=(9, 9))
    wedges, texts, autotexts = ax.pie(
        sizes,
        labels=labels,
        colors=pie_colors,
        autopct='%1.1f%%',
        startangle=90,
        pctdistance=0.82,
        labeldistance=1.05,
        textprops={'fontsize': 11, 'fontweight': 'medium'},
        wedgeprops={'edgecolor': 'white', 'linewidth': 1.2}
    )

    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontweight('bold')

    ax.set_title(f"{target_name} — Feature Importance (Top {n_selected})\nLasso → OLS   |   adj R² = {adj_r2:.3f}   |   n = {n}", fontsize=14, pad=20)

    centre_circle = plt.Circle((0,0), 0.65, fc='white')
    fig.gca().add_artist(centre_circle)

    plt.axis('equal')
    plt.tight_layout()
    safe_name = target_name.replace(' ', '_').replace('%', 'pct')
    plt.savefig(f"{plot_folder}/importance_pie_{safe_name}.png", dpi=220, bbox_inches='tight')
    plt.close()

    fig, ax = plt.subplots(figsize=(6, 6))
    sns.regplot(
        x=y, y=y_pred,
        scatter_kws={'alpha':0.7, 's':70, 'edgecolor':'none'},
        line_kws={'color':'red', 'lw':1.5, 'linestyle':'--'},
        ci=None,
        ax=ax
    )
    ax.set_xlabel(f"Observed {target_name}", fontsize=12)
    ax.set_ylabel(f"Predicted {target_name}", fontsize=12)
    ax.set_title(f"{target_name} — Observed vs Predicted\nn = {n}", fontsize=13, pad=15)

    metrics_text = f"R²   = {r2_val:.3f}\nadj R² = {adj_r2:.3f}\nRMSE = {rmse_val:.3f}\nPCC  = {pcc:.3f}"
    ax.text(
        0.05, 0.95, metrics_text,
        transform=ax.transAxes,
        fontsize=11,
        verticalalignment='top',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.85, edgecolor='gray')
    )

    ax.grid(True, alpha=0.3, linestyle='--')
    ax.set_aspect('equal')
    plt.tight_layout()
    plt.savefig(f"{plot_folder}/obs_vs_pred_with_fit_{safe_name}.png", dpi=220, bbox_inches='tight')
    plt.close()

    res_row = {
        "target": target_name,
        "n": n,
        "n_features": n_selected,
        "adj_R2": adj_r2,
        "R2": r2_val,
        "RMSE": rmse_val,
        "PCC": pcc,
        **{f"coef_{feat}": val for feat, val in coef_series.items()}
    }
    results.append(res_row)

pd.DataFrame(results).to_excel(f"{plot_folder}/summary_lasso_top3.xlsx", index=False)
print(f"\nResults & plots saved to: {plot_folder}")