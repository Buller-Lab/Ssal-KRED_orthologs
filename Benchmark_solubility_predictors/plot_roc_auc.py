import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score, matthews_corrcoef
import argparse

parser = argparse.ArgumentParser(description="Evaluate solubility predictors with ROC curves + cutoff metrics")
parser.add_argument("-i", "--input", required=True,
                    help="Path to the input CSV file (e.g., soluble_expression.csv)")
args = parser.parse_args()

CSV_FILE = args.input
TARGET_COL = "solubly_expressed"

PREDICTOR_COLS = [
    "Protein-sol",
    "pI",
    "SoluProt",
    "NetSolP",
    "NetSolP (ESM1b distilled)",
    "NetSolP (ESM1b)",     
    "NetSolP (ESM12)",        
    "TM-score",
    "SeqID"
]


COLORS = ['#4a84c5', '#8898a9', '#FFA500','#76b689', '#badac4','#49895c','#a382af', '#90CAEE']

CUT_OFFS = [0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8]

df = pd.read_csv(CSV_FILE)

if TARGET_COL not in df.columns:
    raise ValueError(f"Target column '{TARGET_COL}' not found in the CSV")

missing_cols = [col for col in PREDICTOR_COLS if col not in df.columns]
if missing_cols:
    print(f"Warning: The following predictor columns are missing → {missing_cols}")
    PREDICTOR_COLS = [col for col in PREDICTOR_COLS if col in df.columns]

df = df.dropna(subset=[TARGET_COL] + PREDICTOR_COLS).copy()

y_true = df[TARGET_COL].values.astype(int)

print(f"Dataset shape after cleaning: {df.shape}")
print(f"Positive rate (solubly_expressed = 1): {y_true.mean():.3f}\n")

results = []
curves = []

plt.figure(figsize=(9, 7))
plt.plot([0, 1], [0, 1], color='gray', lw=1.0, linestyle='--', alpha=0.7,
         label='Random (AUC = 0.50)')

for i, col in enumerate(PREDICTOR_COLS):
    y_score = df[col].values
    try:
        auc = roc_auc_score(y_true, y_score)
    except ValueError as e:
        print(f"Skipping {col}: {e}")
        continue
    
    fpr, tpr, _ = roc_curve(y_true, y_score)
    curves.append((auc, col, fpr, tpr, COLORS[i % len(COLORS)]))
    results.append({"Predictor": col, "AUC": auc})

curves.sort(reverse=True, key=lambda x: x[0])

for auc, col, fpr, tpr, color in curves:
    plt.plot(fpr, tpr, lw=2, color=color,
             label=f'{col}  ({auc:.2f})')

plt.xlim([-0.02, 1.02])
plt.ylim([-0.02, 1.02])
plt.xlabel('False Positive Rate', fontsize=14, fontweight='bold')
plt.ylabel('True Positive Rate', fontsize=14, fontweight='bold')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc='lower right', fontsize=13, framealpha=0.95)
plt.grid(True, alpha=0.25, linestyle=':')
plt.tight_layout()
plt.savefig("roc_curves_solubility.svg", bbox_inches='tight', dpi=180)

cutoff_rows = []

for predictor in PREDICTOR_COLS:
    if predictor not in df.columns:
        print(f"Cannot compute cutoffs — column missing: {predictor}")
        continue
    
    y_score = df[predictor].values
    
    for cutoff in CUT_OFFS:
        y_pred = (y_score >= cutoff).astype(int)
        
        tp = np.sum((y_pred == 1) & (y_true == 1))
        tn = np.sum((y_pred == 0) & (y_true == 0))
        fp = np.sum((y_pred == 1) & (y_true == 0))
        fn = np.sum((y_pred == 0) & (y_true == 1))
        
        mcc = matthews_corrcoef(y_true, y_pred) if (tp + fp > 0 and tp + fn > 0) else np.nan
        
        row = {
            "Predictor": predictor,
            "Cutoff": cutoff,
            "MCC": round(mcc, 4) if not np.isnan(mcc) else np.nan,
            "TP": int(tp),
            "TN": int(tn),
            "FP": int(fp),
            "FN": int(fn),
            "Total": int(tp + tn + fp + fn),
            "Sensitivity": round(tp / (tp + fn), 4) if (tp + fn) > 0 else 0.0,
            "Specificity": round(tn / (tn + fp), 4) if (tn + fp) > 0 else 0.0,
        }
        cutoff_rows.append(row)

cutoff_df = pd.DataFrame(cutoff_rows)
cutoff_df = cutoff_df.sort_values(["Predictor", "Cutoff"])

print("\n" + "="*80)
print("Performance at different solubility score cutoffs")
print("="*80)
print(cutoff_df.to_string(index=False))
print("="*80 + "\n")

cutoff_df.to_csv("solubility_cutoffs_metrics.csv", index=False)
print("→ Saved detailed cutoff metrics to: solubility_cutoffs_metrics.csv")