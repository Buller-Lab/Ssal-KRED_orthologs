import sys, warnings, pandas as pd, numpy as np
from scipy.stats import pearsonr, ConstantInputWarning
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap

TARGET_HOMOLOG_ID = "Q9UUN9"
TARGET_LABEL = "Ssal-KRED"

def sig_stars(p):
    if p < 0.001: return "***"
    elif p < 0.01: return "**"
    elif p < 0.05: return "*"
    else: return "ns"

def compute_correlations(df, target_col="melting_temperature"):
    results = []
    feature_cols = [c for c in df.columns if c not in ["homolog_id", target_col]]
    df_target_row = df[df["homolog_id"] == TARGET_HOMOLOG_ID]

    for col in feature_cols:
        x = pd.to_numeric(df[col], errors="coerce")
        y = pd.to_numeric(df[target_col], errors="coerce")
        mask = x.notna() & y.notna()
        x, y = x[mask], y[mask]

        if len(x) < 3 or x.nunique() < 2 or np.std(x) == 0 or np.std(y) == 0: continue

        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", ConstantInputWarning)
                r, p = pearsonr(x, y)
        except: continue

        top_rank = bottom_rank = np.nan

        if not df_target_row.empty:
            v = pd.to_numeric(df_target_row[col], errors="coerce").values[0]
            if not np.isnan(v):
                tmp = df[["homolog_id", col]].copy()
                tmp[col] = pd.to_numeric(tmp[col], errors="coerce")
                tmp = tmp.dropna(subset=[col])

                if tmp[col].nunique() > 1:
                    t1 = tmp.sort_values(col, ascending=False).reset_index(drop=True)
                    t1["rank"] = np.arange(1, len(t1) + 1)
                    m1 = t1[t1["homolog_id"] == TARGET_HOMOLOG_ID]
                    if not m1.empty: top_rank = int(m1["rank"].values[0])

                    t2 = tmp.sort_values(col, ascending=True).reset_index(drop=True)
                    t2["rank"] = np.arange(1, len(t2) + 1)
                    m2 = t2[t2["homolog_id"] == TARGET_HOMOLOG_ID]
                    if not m2.empty: bottom_rank = int(m2["rank"].values[0])

        results.append({"feature": col, "pearson_r": r, "p_value": p,
                        "top_rank": top_rank, "bottom_rank": bottom_rank})

    return pd.DataFrame(results)

def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py <features_csv> <target_csv>")
        sys.exit(1)

    df = pd.merge(pd.read_csv(sys.argv[1]), pd.read_csv(sys.argv[2]),
                on="homolog_id", how="inner")
    df = df.dropna(subset=["melting_temperature"])

    corr = compute_correlations(df)
    corr = corr.reindex(corr["pearson_r"].abs().sort_values(ascending=False).index)

    print("\n" + "="*80)
    print(" SIGNIFICANT CORRELATION OBSERVATIONS (p < 0.05)")
    print("="*80)

    significant_found = False
    
    for _, row in corr.iterrows():
        p_val = row["p_value"]
        if p_val >= 0.05:
            continue
            
        significant_found = True
        feature = row["feature"]
        r_val = row["pearson_r"]
        stars = sig_stars(p_val)
        t_rank = row["top_rank"]
        b_rank = row["bottom_rank"]
        
        direction = "POSITIVE" if r_val > 0 else "NEGATIVE"
        
        print(f"\n• Feature: {feature}")
        print(f"  Correlation: {direction} (r = {r_val:.3f}, p = {p_val:.4e} {stars})")
        
        if direction == "POSITIVE" and not pd.isna(t_rank):
            print(f"  Enzyme Position: {TARGET_LABEL} ranks as Top {int(t_rank)} for this feature.")
            if int(t_rank) <= 5:
                print(f"    -> Observation: Highly optimized! High values link to thermal stability, and {TARGET_LABEL} sits at the top tier.")
        elif direction == "NEGATIVE" and not pd.isna(b_rank):
            print(f"  Enzyme Position: {TARGET_LABEL} ranks as Bottom {int(b_rank)} for this feature (where 1 is the lowest absolute value).")
            if int(b_rank) <= 5:
                print(f"    -> Observation: Favorable configuration! Lower values link to thermal stability, and {TARGET_LABEL} successfully minimizes this descriptor.")
        else:
            if not pd.isna(t_rank):
                print(f"  Enzyme Position: {TARGET_LABEL} ranks as Top {int(t_rank)} (Ranked high, opposing the target axis).")

    if not significant_found:
        print("\nNo statistically significant correlations (p < 0.05) were identified in this dataset.")
    print("="*80 + "\n")

    top = corr.head(10).iloc[::-1]

    if top.empty:
        print("No correlations found to plot.")
        return

    cmap = LinearSegmentedColormap.from_list("custom_div", ["#885A95", "#ffffff", "#76b689"])
    norm = plt.Normalize(vmin=-1, vmax=1)
    bar_colors = cmap(norm(top["pearson_r"].values))

    plt.figure(figsize=(13, 7))
    bars = plt.barh(top["feature"], top["pearson_r"], color=bar_colors)
    ax = plt.gca()

    used_stars = set()

    for bar, r, p, t, b in zip(bars, top["pearson_r"], top["p_value"], top["top_rank"], top["bottom_rank"]):
        label = None
        if r >= 0 and not pd.isna(t): label = f"{TARGET_LABEL}=Top {int(t)}"
        elif r < 0 and not pd.isna(b): label = f"{TARGET_LABEL}=Bottom {int(b)}"

        if label:
            ax.text(bar.get_width()/2,
                    bar.get_y()+bar.get_height()/2,
                    label,
                    va="center", ha="center",
                    fontsize=9, color="black", fontweight="bold")

        stars = sig_stars(p)
        used_stars.add(stars)

        ax.text(bar.get_width(),
                bar.get_y()+bar.get_height()/2,
                f"{r:.2f}\n({stars})",
                va="center",
                ha="left" if r >= 0 else "right",
                fontsize=9, color="black", fontweight="bold")

    plt.xlabel("Pearson correlation (r)")
    plt.title("Top 10 feature correlations with melting temperature")

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, orientation="vertical", shrink=0.75, pad=0.08)
    cbar.set_label("Correlation Scale (r)", fontsize=10, fontweight="bold")

    legend_map = {
        "***": "p < 0.001",
        "**": "p < 0.01",
        "*": "p < 0.05",
        "ns": "p ≥ 0.05"
    }

    handles = [
        Line2D([0], [0], marker='', color='w', label=f"{k} = {legend_map[k]}")
        for k in ["***", "**", "*", "ns"] if k in used_stars
    ]

    plt.legend(handles=handles, loc="lower right", frameon=True)

    plt.tight_layout()
    plt.savefig("correlations_top10.png", dpi=300)
    plt.close()

if __name__ == "__main__":
    main()