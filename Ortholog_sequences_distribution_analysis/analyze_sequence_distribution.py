import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import squareform
from scipy.stats import skew, kurtosis

def calculate_percent_identity(seq1, seq2):
    if len(seq1) == 0 or len(seq2) == 0:
        return 0
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 1.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -5.0
    aligner.extend_gap_score = -1.0
    alignments = aligner.align(seq1, seq2)
    alignment = alignments[0]
    aligned_seq1 = alignment[0]
    aligned_seq2 = alignment[1]
    matches = 0
    total = 0
    for a, b in zip(aligned_seq1, aligned_seq2):
        if a != '-' and b != '-':
            total += 1
            if a == b:
                matches += 1
    return (matches / total * 100) if total > 0 else 0

parser = argparse.ArgumentParser(description='Quality control for protein sequences in FASTA file.')
parser.add_argument('-i', required=True, help='Input FASTA file')
args = parser.parse_args()

ids = []
sequences = []
for record in SeqIO.parse(args.i, 'fasta'):
    ids.append(record.id)
    sequences.append(str(record.seq))

n = len(sequences)
if n == 0:
    raise ValueError("No sequences found in the input file.")

lengths = [len(seq) for seq in sequences]

identity_matrix = np.zeros((n, n))
for i in range(n):
    identity_matrix[i, i] = 100
    for j in range(i + 1, n):
        identity = calculate_percent_identity(sequences[i], sequences[j])
        identity_matrix[i, j] = identity
        identity_matrix[j, i] = identity

dist_matrix = 100 - identity_matrix
condensed_dist = squareform(dist_matrix)

Z = linkage(condensed_dist, method='average')
clusters = fcluster(Z, t=20, criterion='distance')

df = pd.DataFrame({
    'ID': ids,
    'Length': lengths,
    'Cluster': clusters
})

df.to_csv('sequence_distribution_table.csv', index=False)
print("CSV table saved as 'sequence_distribution_table.csv'")

pairwise_identities = identity_matrix[np.triu_indices(n, k=1)]
cluster_counts = pd.Series(clusters).value_counts().sort_index()

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans', 'Liberation Sans'],
    'axes.titlesize': 13,
    'axes.labelsize': 11,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'figure.titlesize': 14,
})

fig = plt.figure(figsize=(14, 11.2))

ax_dend = plt.subplot2grid((5, 2), (0, 0), colspan=2, rowspan=2)
dendrogram(Z, labels=ids, leaf_rotation=90, ax=ax_dend, above_threshold_color='grey')
ax_dend.set_title('a.  Hierarchical clustering dendrogram', loc='left', pad=12)
ax_dend.set_ylabel('Distance (100 − % identity)')
ax_dend.set_xlabel('Sequence IDs')

ax_cluster = plt.subplot2grid((5, 2), (2, 0), colspan=2)
ax_cluster.bar(cluster_counts.index.astype(str), cluster_counts.values, edgecolor='black', linewidth=0.6)
ax_cluster.set_title('b.  Cluster size distribution  (>80% identity clusters)', loc='left', pad=12)
ax_cluster.set_xlabel('Cluster ID')
ax_cluster.set_ylabel('Number of sequences')
ax_cluster.tick_params(axis='x', rotation=90, labelsize=9.5)

ax_length = plt.subplot2grid((5, 2), (3, 0), rowspan=2)
ax_length.hist(lengths, bins=min(30, n//5 + 1), edgecolor='black', linewidth=0.6)
ax_length.set_title('c.  Sequence length distribution', loc='left', pad=10)
ax_length.set_xlabel('Length (amino acids)')
ax_length.set_ylabel('Frequency')

ax_identity = plt.subplot2grid((5, 2), (3, 1), rowspan=2)
ax_identity.hist(pairwise_identities, bins=30, edgecolor='black', linewidth=0.6)
ax_identity.set_title('d.  Pairwise % identity distribution', loc='left', pad=10)
ax_identity.set_xlabel('% identity')
ax_identity.set_ylabel('Frequency')

plt.tight_layout(h_pad=1.4, rect=[0, 0, 1, 0.98])
plt.savefig('sequence_qc_summary.png', dpi=180, bbox_inches='tight')
print("Main summary figure saved as 'sequence_qc_summary.png'")