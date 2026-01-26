import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score
import argparse

parser = argparse.ArgumentParser(description="Evaluate solubility predictors with ROC curves")
parser.add_argument("-i", "--input", required=True,
                    help="Path to the input CSV file (e.g., soluble_expression.csv)")
args = parser.parse_args()

CSV_FILE = args.input
TARGET_COL = "solubly_expressed"

PREDICTOR_COLS = [
    "Protein-sol",
    "pI",
    "SoluProt",
    "TM-score",
    "SeqID"
]

COLORS = ['#4a84c5', '#8898a9', '#76b689', '#a382af', '#90CAEE']

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

# Sort curves by AUC descending → best on top in legend
curves.sort(reverse=True, key=lambda x: x[0])

for auc, col, fpr, tpr, color in curves:
    plt.plot(fpr, tpr, lw=2.2, color=color,
             label=f'{col}  ({auc:.3f})')

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
plt.show()