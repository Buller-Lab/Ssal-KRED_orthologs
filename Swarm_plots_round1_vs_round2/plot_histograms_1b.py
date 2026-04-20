import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

MAIN_COLOR = '#d2e0f3'
SECOND_COLOR = '#73b587'

if len(sys.argv) < 3 or len(sys.argv) > 4:
    print("Usage:")
    print(" python script.py <csv_file> <main_column>")
    print(" python script.py <csv_file> <main_column> <second_column>")
    sys.exit(1)

file_path = sys.argv[1]
main_col = sys.argv[2]
second_col = sys.argv[3] if len(sys.argv) == 4 else None

columns = [main_col, 'refined']
if second_col:
    columns.append(second_col)

df = pd.read_csv(file_path)[columns].dropna()
df['Round'] = df['refined'].map({0: 'Round 1', 1: 'Round 2'})

fig, (ax_left, ax_right) = plt.subplots(
    nrows=1,
    ncols=2,
    figsize=(10.5, 5.4),
    sharey=False,
    gridspec_kw={'wspace': 0.08}
)

sns.swarmplot(
    data=df,
    x='Round',
    y=main_col,
    hue='Round',
    palette=[MAIN_COLOR, MAIN_COLOR],
    size=15,
    alpha=0.90,
    edgecolor='black',
    linewidth=1.0,
    dodge=False,
    legend=False,
    ax=ax_left
)

ax_left.set_ylabel("Tm (°C)", fontsize=18, fontweight='bold', color='black')
ax_left.tick_params(axis='both', labelsize=18)

for label in ax_left.get_yticklabels():
    label.set_fontweight('bold')
for label in ax_left.get_xticklabels():
    label.set_fontweight('bold')

ax_left.spines['left'].set_linewidth(2.0)
ax_left.spines['bottom'].set_linewidth(2.0)
ax_left.spines['top'].set_visible(False)
ax_left.spines['right'].set_visible(False)

if second_col:
    sns.swarmplot(
        data=df,
        x='Round',
        y=second_col,
        hue='Round',
        palette=[SECOND_COLOR, SECOND_COLOR],
        size=16,
        alpha=0.90,
        edgecolor='black',
        linewidth=0.5,
        dodge=False,
        legend=False,
        ax=ax_right
    )

    ax_right.set_ylabel("FIO_Ssal-KRED 1b", fontsize=18, fontweight='bold', color='black')
    ax_right.tick_params(axis='both', labelsize=18)

    for label in ax_right.get_yticklabels():
        label.set_fontweight('bold')
    for label in ax_right.get_xticklabels():
        label.set_fontweight('bold')

    ax_right.spines['left'].set_visible(False)
    ax_right.spines['right'].set_linewidth(2.0)
    ax_right.spines['bottom'].set_linewidth(2.0)
    ax_right.spines['top'].set_visible(False)

    ax_right.yaxis.set_label_position("right")
    ax_right.yaxis.tick_right()

else:
    ax_right.axis('off')

sns.despine(fig=fig, top=True, right=False, left=False, bottom=False)

ax_left.set_xlabel('')
ax_right.set_xlabel('')

plt.tight_layout(rect=[0, 0, 1, 0.96])

safe_main = main_col.replace(' ', '_').replace('/', '_').replace('\\', '_')
output_name = f"Round_comparison_{safe_main}"
if second_col:
    safe_sec = second_col.replace(' ', '_').replace('/', '_').replace('\\', '_')
    output_name += f"_{safe_sec}"
output_name += ".svg"

plt.savefig(output_name, format='svg', bbox_inches='tight', dpi=300)
print(f"Figure saved as: {output_name}")
plt.show()