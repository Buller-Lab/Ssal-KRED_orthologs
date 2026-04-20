import argparse
import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression, LassoCV
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import LeaveOneOut
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import os
import pickle
import json

FEATURES = [
    "sequence_identity", "rmsd", "tm_score", "rmsd_bindingsite",
    "n_bindingsite_residues", "n_aligned_residues", "aligned_fraction", "bindingsite_similarity_score",
    "non_equivalent_residues", "cavity_volume", "cavity_area", "cavity_avg_depth", "cavity_avg_hydropathy",
    "cavity_frequency_Alipathic_apolar", "cavity_frequency_Aromatic", "cavity_frequency_Polar_uncharged",
    "cavity_frequency_Negatively_charged", "cavity_frequency_Positively_charged",
    "cavity_frequency_Non-standard", "cavity_rmsd", "cavity_n_points"
]

parser = argparse.ArgumentParser(description="LASSO + MLR with top-3 features, LOO-CV, and model saving")
parser.add_argument("-i", "--input", required=True, help="Input CSV file")
parser.add_argument("-t", "--target", required=True, help="Name of the target column to predict (e.g. '2b (R)')")
args = parser.parse_args()

plot_folder = "mlr_lasso_top3_fixed_features_LOO"
model_folder = "saved_models"
os.makedirs(plot_folder, exist_ok=True)
os.makedirs(model_folder, exist_ok=True)

df = pd.read_csv(args.input)

if args.target not in df.columns:
    print(f"Error: Target column '{args.target}' not found in the input file.")
    exit(1)

available_features = [f for f in FEATURES if f in df.columns]
df = df[available_features + [args.target]].dropna().reset_index(drop=True)

X_raw = df[available_features]
y = df[args.target].values
y_name = args.target

scaler = StandardScaler()
X_std = scaler.fit_transform(X_raw)
feature_names = available_features

COLORS = ['#1f77b4', '#76b689', '#d62728']

lasso = LassoCV(cv=5, random_state=42, max_iter=30000, alphas=200, n_jobs=-1)
lasso.fit(X_std, y)

coef = pd.Series(lasso.coef_, index=feature_names)
top_features = coef.abs().sort_values(ascending=False).head(3).index.tolist()
n_selected = len(top_features)

if n_selected == 0:
    print(f"No features selected for target: {y_name}")
    exit(0)

sel_idx = [feature_names.index(f) for f in top_features]
X_sel = X_std[:, sel_idx]

model_full = LinearRegression().fit(X_sel, y)
y_pred_full = model_full.predict(X_sel)

r2_full   = r2_score(y, y_pred_full)
rmse_full = np.sqrt(mean_squared_error(y, y_pred_full))
adj_r2    = 1 - (1 - r2_full) * (len(y) - 1) / (len(y) - n_selected - 1) if len(y) > n_selected + 1 else r2_full

cv = LeaveOneOut()
y_pred_cv = np.zeros_like(y, dtype=float)
ols = LinearRegression()

for train_idx, test_idx in cv.split(X_sel):
    ols.fit(X_sel[train_idx], y[train_idx])
    y_pred_cv[test_idx] = ols.predict(X_sel[test_idx])

press     = np.sum((y - y_pred_cv)**2)
tss       = np.sum((y - np.mean(y))**2)
q2        = 1 - press / tss if tss > 1e-10 else np.nan
rmse_loo  = np.sqrt(mean_squared_error(y, y_pred_cv))
mae_loo   = np.mean(np.abs(y - y_pred_cv))
pcc_loo, _ = pearsonr(y, y_pred_cv)

delta = r2_full - q2 if not np.isnan(q2) else np.nan

print(f"{'═' * 70}")
print(f"Target: {y_name}")
print(f"Number of samples: {len(y)}")
print(f"{'─' * 70}")
print("Top 3 features (LASSO):")
for i, feat in enumerate(top_features, 1):
    beta = coef[feat]
    print(f"  {i}. {feat} (β = {beta:>6.3f})")
print(f"\nFull-data fit:")
print(f"  R²     = {r2_full:>6.3f}")
print(f"  adj R² = {adj_r2:>6.3f}")
print(f"  RMSE   = {rmse_full:>6.3f}")
print(f"\nLOO-CV performance:")
print(f"  Q²     = {q2:>6.3f}")
print(f"  RMSE   = {rmse_loo:>6.3f}")
print(f"  MAE    = {mae_loo:>6.3f}")
print(f"  PCC    = {pcc_loo:>6.3f}")
print(f"  Δ(R² - Q²) = {delta:>6.3f}")
print("\nInterpretation:")
if np.isnan(q2):
    print("  → Q² could not be computed (very low variance in y?)")
elif delta < 0.05:
    print("  → Excellent internal predictivity — very little overfitting")
elif delta < 0.15:
    print("  → Good internal predictivity — mild overfitting acceptable")
elif delta < 0.30:
    print("  → Moderate overfitting detected — still potentially useful")
else:
    print("  → Strong overfitting or limited generalizability")
print(f"{'═' * 70}\n")

safe_name = y_name.replace(" ", "_").replace("(", "").replace(")", "").replace("%", "pct").replace("/", "_")

model_data = {
    "target": y_name,
    "top_features": top_features,
    "scaler": scaler,
    "model": model_full,
    "feature_names": feature_names,
    "selected_indices": sel_idx,
    "coef": model_full.coef_.tolist(),
    "intercept": float(model_full.intercept_),
    "training_samples": len(y),
    "r2_full": float(r2_full),
    "q2_loo": float(q2) if not np.isnan(q2) else None,
}

model_path = os.path.join(model_folder, f"model_{safe_name}.pkl")
with open(model_path, "wb") as f:
    pickle.dump(model_data, f)

json_path = os.path.join(model_folder, f"model_{safe_name}_info.json")
with open(json_path, "w") as f:
    json.dump({
        "target": y_name,
        "top_features": top_features,
        "coefficients": dict(zip(top_features, model_full.coef_.tolist())),
        "intercept": float(model_full.intercept_),
        "r2": float(r2_full),
        "q2_loo": float(q2) if not np.isnan(q2) else None,
        "n_samples": len(y),
        "model_file": f"model_{safe_name}.pkl"
    }, f, indent=4)

print(f"Model successfully saved as: {model_path}")
print(f"Model info summary saved as: {json_path}")

abs_coef = coef[top_features].abs().values
abs_coef_norm = abs_coef / abs_coef.sum() if abs_coef.sum() > 0 else abs_coef

fig_pie, ax_pie = plt.subplots(figsize=(7, 5.5))
wedges, _, autotexts = ax_pie.pie(
    abs_coef_norm,
    labels=None,
    colors=COLORS[:n_selected],
    autopct='%1.1f%%' if abs_coef.sum() > 0 else None,
    startangle=90,
    textprops={'fontsize': 11},
    pctdistance=0.78,
    wedgeprops={'edgecolor': 'white', 'linewidth': 1.2}
)

ax_pie.legend(
    wedges,
    [f.replace("_", " ").title() for f in top_features],
    title="Top features (LASSO |β|)",
    loc="center left",
    bbox_to_anchor=(1.0, 0.5, 0.35, 0),
    fontsize=10.5,
    title_fontsize=11.5,
    frameon=True,
    edgecolor='black',
    facecolor='white',
    labelspacing=0.8
)

centre_circle = plt.Circle((0,0), 0.55, fc='white', linewidth=1.2, edgecolor='lightgray')
ax_pie.add_artist(centre_circle)

ax_pie.set_title(f"Feature importance – {y_name}\n(LASSO absolute coefficients)",
                 fontsize=13, pad=18)

plt.tight_layout()
plt.savefig(os.path.join(plot_folder, f"pie_importance_{safe_name}.png"), dpi=300, bbox_inches="tight")
plt.savefig(os.path.join(plot_folder, f"pie_importance_{safe_name}.svg"), format="svg", bbox_inches="tight")
plt.close(fig_pie)

fig_scatter, ax_scatter = plt.subplots(figsize=(6.5, 6))
ax_scatter.scatter(y, y_pred_full, color='#1f77b4', alpha=0.7, edgecolor='white', s=60, linewidth=0.8)

min_val = min(np.min(y), np.min(y_pred_full))
max_val = max(np.max(y), np.max(y_pred_full))
ax_scatter.plot([min_val, max_val], [min_val, max_val], color='gray', linestyle='--', linewidth=1.5, alpha=0.7)

ax_scatter.set_xlabel(f"Observed {y_name}", fontsize=12)
ax_scatter.set_ylabel(f"Predicted {y_name}", fontsize=12)
ax_scatter.set_title(f"Linear Regression Fit – {y_name}\n(R² = {r2_full:.3f})", fontsize=13, pad=15)

ax_scatter.grid(True, linestyle='--', alpha=0.4)
ax_scatter.set_aspect('equal')

plt.tight_layout()
plt.savefig(os.path.join(plot_folder, f"scatter_fit_{safe_name}.png"), dpi=300, bbox_inches="tight")
plt.savefig(os.path.join(plot_folder, f"scatter_fit_{safe_name}.svg"), format="svg", bbox_inches="tight")
plt.close(fig_scatter)

print(f"Plots (pie + scatter) saved to: {plot_folder}")