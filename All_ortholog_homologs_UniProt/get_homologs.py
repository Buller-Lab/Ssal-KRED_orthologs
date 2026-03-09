import os
import pandas as pd
import subprocess
import sys
import argparse
from pathlib import Path

parser = argparse.ArgumentParser(description='MMseqs2 easy-search + linclust with per-seed %id and BacDive trait columns')
parser.add_argument('-i', '--input', required=True, help='Input seeds FASTA file')
parser.add_argument('-d', '--database', required=True, help='Path to MMseqs2 database')
parser.add_argument('-o', '--output-dir', required=True, help='Output directory / project name')
parser.add_argument('--bacdive-dir', help='Path to bacdive_taxid folder containing *_ncbi.csv files (optional)')
args = parser.parse_args()

input_fasta = args.input
database = args.database
output_dir = args.output_dir
bacdive_dir = args.bacdive_dir

if not os.path.exists(input_fasta):
    print(f"Error: Input file not found: {input_fasta}")
    sys.exit(1)

os.makedirs(output_dir, exist_ok=True)

trait_sets = {}
trait_column_order = []

if bacdive_dir and os.path.isdir(bacdive_dir):
    for file in Path(bacdive_dir).glob("*_ncbi.csv"):
        trait_name = file.stem.replace("_ncbi", "").strip()
        try:
            df_trait = pd.read_csv(file, dtype={'tax_id': int})
            val_col = [c for c in df_trait.columns if c != 'tax_id'][0]
            true_rows = df_trait[val_col].astype(str).str.upper().isin(['TRUE', 'T', 'YES', '1'])
            taxid_set = set(df_trait.loc[true_rows, 'tax_id'])
            trait_sets[trait_name] = taxid_set
            trait_column_order.append(trait_name)
        except:
            pass

format_fields = "query,target,tseq,pident,tlen,taxid,taxname,taxlineage"

subprocess.run(
    [
        "mmseqs", "easy-search",
        input_fasta, database,
        f"{output_dir}/alnRes.m8", "tmp",
        "--format-output", format_fields,
        "--min-seq-id", "0.25",
        "--max-seqs", "100000000000",
        "--split-memory-limit", "90G"
    ],
    check=True
)

hits = pd.read_csv(
    f'{output_dir}/alnRes.m8',
    sep='\t',
    header=None,
    names=['query', 'target', 'tseq', 'pident', 'tlen', 'taxid', 'taxname', 'taxlineage'],
    comment='#',
    dtype={'taxid': 'Int64'}
)

hits['seed_accession'] = hits['query'].str.split().str[0]
hits['col_name'] = hits['seed_accession'].apply(lambda x: f"SeqID ({x})")

for trait, taxset in trait_sets.items():
    hits[trait] = hits['taxid'].isin(taxset)

wide_identities = hits.pivot_table(
    index='target',
    columns='col_name',
    values='pident',
    aggfunc='first'
).reset_index()

hits_best = hits.sort_values(['target', 'pident'], ascending=[True, False]) \
                .drop_duplicates('target', keep='first')

df_final = hits_best.merge(wide_identities, on='target', how='left')

df_final = df_final.rename(columns={
    'tseq':    'Sequence',
    'tlen':    'Length',
    'taxname': 'Organism'
})

for trait in trait_column_order:
    if trait not in df_final.columns:
        df_final[trait] = df_final['taxid'].isin(trait_sets.get(trait, set()))

extremophile_traits = [
    'acidophilic',
    'alkaliphilic',
    'halophilic',
    'psychrophilic',
    'thermophilic',
]

def create_annotation(row):
    present = []
    for trait in extremophile_traits:
        if trait in row and row[trait] is True:
            present.append(trait)
    if 'thermophilic' in present and 'psychrophilic' in present:
        present = [t for t in present if t not in ('thermophilic', 'psychrophilic')]
    if 'alkaliphilic' in present and 'acidophilic' in present:
        present = [t for t in present if t not in ('alkaliphilic', 'acidophilic')]
    if len(present) == 0:
        return ""
    elif len(present) == 1:
        return present[0]
    else:
        return " & ".join(sorted(present))

if any(t in df_final.columns for t in extremophile_traits):
    df_final['Annotation'] = df_final.apply(create_annotation, axis=1)

id_cols = [c for c in df_final.columns if c.startswith("SeqID (")]
trait_cols = [t for t in trait_column_order if t in df_final.columns]
core_cols = ['target', 'Sequence', 'Length', 'taxid', 'Organism', 'taxlineage']
new_order = core_cols + trait_cols + ['Annotation'] + id_cols + ['query']
df_final = df_final[[c for c in new_order if c in df_final.columns]]

hits_renamed = hits.rename(columns={
    'tseq':    'Sequence',
    'tlen':    'Length',
    'taxname': 'Organism',
    'pident':  'SeqID (%)',
    'query':   'Seed',
    'target':  'target'
})

if any(t in hits_renamed.columns for t in extremophile_traits):
    hits_renamed['Annotation'] = hits_renamed.apply(create_annotation, axis=1)

long_cols = ['Seed', 'target', 'Sequence', 'SeqID (%)', 'Length', 'taxid', 'Organism', 'taxlineage']
trait_cols_long = [t for t in trait_column_order if t in hits_renamed.columns]
final_long_cols = long_cols + trait_cols_long + ['Annotation']
hits_renamed = hits_renamed[[c for c in final_long_cols if c in hits_renamed.columns]]

df_final.to_csv(
    f'{output_dir}/alnRes_wide.m8',
    sep='\t',
    index=False,
    header=True
)

hits_renamed.to_csv(
    f'{output_dir}/alnRes.m8',
    sep='\t',
    index=False,
    header=True
)

with open(f'{output_dir}/hits.fasta', 'w') as f:
    for _, row in hits.iterrows():
        f.write(f">{row['target']}\n{row['tseq']}\n")

subprocess.run(
    [
        "mmseqs", "easy-linclust",
        f"{output_dir}/hits.fasta",
        f"{output_dir}/clusterRes_60", "tmp",
        "--min-seq-id", "0.6", "-c", "0.8", "--cov-mode", "0"
    ],
    check=True
)

cluster_60 = pd.read_csv(
    f'{output_dir}/clusterRes_60_cluster.tsv',
    sep='\t', header=None, names=['representative', 'member'], comment='#'
)

cluster_60.to_csv(
    f'{output_dir}/clusterRes_60_cluster.tsv',
    sep='\t', index=False, header=['representative', 'member']
)

print("Done.")