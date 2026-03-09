import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import squareform
from collections import Counter
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

cluster_counts = pd.Series(clusters).value_counts()

standard_aas = 'ACDEFGHIKLMNPQRSTVWY'
aa_freqs = []
for seq in sequences:
    count = Counter(seq)
    total = len(seq)
    freq = {aa: count.get(aa, 0) / total for aa in standard_aas}
    aa_freqs.append(freq)
df_aa = pd.DataFrame(aa_freqs)
mean_aa = df_aa.mean().sort_values(ascending=False)

fig, axs = plt.subplots(2, 2, figsize=(14, 10))

axs[0, 0].hist(lengths, bins=min(30, n//5 + 1), edgecolor='black')
axs[0, 0].set_title('Sequence Length Distribution')
axs[0, 0].set_xlabel('Length')
axs[0, 0].set_ylabel('Frequency')

axs[0, 1].hist(pairwise_identities, bins=30, edgecolor='black')
axs[0, 1].set_title('Pairwise % Identity Distribution')
axs[0, 1].set_xlabel('% Identity')
axs[0, 1].set_ylabel('Frequency')

axs[1, 0].bar(cluster_counts.index.astype(str), cluster_counts.values, edgecolor='black')
axs[1, 0].set_title('Cluster Size Distribution (>80% Identity Clusters)')
axs[1, 0].set_xlabel('Cluster ID')
axs[1, 0].set_ylabel('Number of Sequences')
axs[1, 0].tick_params(axis='x', rotation=90)

axs[1, 1].bar(mean_aa.index, mean_aa.values, edgecolor='black')
axs[1, 1].set_title('Average Amino Acid Composition')
axs[1, 1].set_xlabel('Amino Acid')
axs[1, 1].set_ylabel('Mean Frequency')

plt.tight_layout()
plt.savefig('sequence_distribution_figure.png')
print("Figure saved as 'sequence_distribution_figure.png'")

plt.figure(figsize=(10, 7))
dendrogram(Z, labels=ids, leaf_rotation=90)
plt.title('Hierarchical Clustering Dendrogram')
plt.xlabel('Sequence IDs')
plt.ylabel('Distance (100 - %Identity)')
plt.savefig('sequence_distribution_dendrogram.png')
print("Dendrogram saved as 'sequence_distribution_dendrogram.png'")

print("\nDataset Statistics for ML Suitability:")

print("\nSequence Length Statistics:")
print(f"Number of sequences: {n}")
print(f"Min length: {np.min(lengths)}")
print(f"Max length: {np.max(lengths)}")
print(f"Mean length: {np.mean(lengths):.2f}")
print(f"Median length: {np.median(lengths):.2f}")
print(f"Std length: {np.std(lengths):.2f}")
print(f"Skewness: {skew(lengths):.2f}")
print(f"Kurtosis: {kurtosis(lengths):.2f}")

print("\nPairwise Identity Statistics (non-self pairs):")
print(f"Number of pairs: {len(pairwise_identities)}")
print(f"Min %ID: {np.min(pairwise_identities):.2f}")
print(f"Max %ID: {np.max(pairwise_identities):.2f}")
print(f"Mean %ID: {np.mean(pairwise_identities):.2f}")
print(f"Median %ID: {np.median(pairwise_identities):.2f}")
print(f"Std %ID: {np.std(pairwise_identities):.2f}")

print("\nClustering Statistics (threshold >80% avg. identity):")
print(f"Number of clusters: {len(cluster_counts)}")
print(f"Number of singletons: {(cluster_counts == 1).sum()}")
print(f"Max cluster size: {cluster_counts.max()}")
print(f"Mean cluster size: {cluster_counts.mean():.2f}")
print(f"Percentage of sequences in singletons: {((cluster_counts == 1).sum() / n * 100):.2f}%")
print(f"Percentage of sequences in largest cluster: {(cluster_counts.max() / n * 100):.2f}%")

print("\nAverage Amino Acid Frequencies:")
for aa, freq in mean_aa.items():
    print(f"{aa}: {freq:.4f}")