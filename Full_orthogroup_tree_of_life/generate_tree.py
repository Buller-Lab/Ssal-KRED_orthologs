import sys
import os
import subprocess
import tempfile
from Bio import AlignIO, SeqIO

if sys.version_info < (3, 6):
    print("Error: Python 3.6 or higher is required.")
    sys.exit(1)

if len(sys.argv) != 2:
    print("Usage: python generate_tree.py your_sequences.fasta")
    sys.exit(1)

fasta_file = sys.argv[1]
output_base = fasta_file.replace('.fasta', '_pretty_cyclic_tree')

records = list(SeqIO.parse(fasta_file, "fasta"))
seq_lengths = [len(rec.seq) for rec in records]

if len(set(seq_lengths)) > 1:
    print("Aligning sequences with MAFFT...")
    with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as temp_aligned:
        aligned_fasta = temp_aligned.name
    subprocess.run(["mafft", "--auto", "--leavegappyregion", fasta_file], 
                   stdout=open(aligned_fasta, "w"), check=True)
else:
    print("Sequences are already aligned.")
    aligned_fasta = fasta_file

alignment = AlignIO.read(aligned_fasta, "fasta")
seq_dict = {rec.id: str(rec.seq) for rec in alignment}

ref_seq = list(seq_dict.values())[0]
groups = {}
for name, seq in seq_dict.items():
    matches = sum(a == b for a, b in zip(ref_seq, seq) if a != '-' and b != '-')
    total = sum(1 for a, b in zip(ref_seq, seq) if a != '-' and b != '-')
    sim = (matches / total * 100) if total > 0 else 0
    if sim >= 50:
        group = "High"
    elif sim >= 40:
        group = "Medium"
    elif sim >= 30:
        group = "Low"
    else:
        group = "Far distant"
    groups[name] = (group, sim)

print("Running FastTree (protein mode)...")
tree_file = output_base + ".nwk"

try:
    with open(tree_file, "w") as f_out:
        subprocess.run([
            "FastTree",
            "-lg",
            "-gamma",
            aligned_fasta
        ], stdout=f_out, check=True, text=True)
except Exception as e:
    print(f"FastTree failed: {e}")
    print("Make sure FastTree is installed and in your PATH (e.g. conda install -c bioconda fasttree)")
    sys.exit(1)

annot_file = output_base + "_itol_colors.txt"
with open(annot_file, "w") as f:
    f.write("DATASET_COLORSTRIP\n")
    f.write("SEPARATOR TAB\n")
    f.write("DATASET_LABEL\tSimilarity\n")
    f.write("COLOR\t#ff0000\n")
    f.write("DATA\n")
    for seq_id, (group, sim) in groups.items():
        safe_id = seq_id.replace(":", "_").replace("|", "_").replace("(", "_").replace(")", "_").replace(",", "_")
        if group == "High":
            color = "#59C959"
        elif group == "Medium":
            color = "#EAFF80"
        elif group == "Low":
            color = "#87CEEB"
        else:
            color = "#905799"
        f.write(f"{safe_id}\t{color}\t{group} ({sim:.1f}%)\n")

print(f"Tree and annotation files generated: {tree_file} and {annot_file}")
print("Upload to iTOL:")
print(f"1. Upload {tree_file}")
print(f"2. Add {annot_file} as Color strip dataset")
print("3. Set display mode to Circular")
print(f"Files are in: {os.getcwd()}")

if aligned_fasta != fasta_file and os.path.exists(aligned_fasta):
    os.remove(aligned_fasta)