#!/usr/bin/env python3
import sys
import pandas as pd

def main():
    if len(sys.argv) not in (5, 7):
        print(f"Usage: {sys.argv[0]} input.m8 min_len max_len output_prefix [--max_seqid FLOAT]")
        print("   --max_seqid : upper cutoff for any SeqID* column (e.g. 90.0 means ≥90% identity → filtered out)")
        print("                 if not given → no sequence identity filtering")
        sys.exit(1)

    input_file = sys.argv[1]
    min_len = int(sys.argv[2])
    max_len = int(sys.argv[3])
    prefix = sys.argv[4]

    max_seqid = None
    if len(sys.argv) == 7:
        if sys.argv[5] != "--max_seqid":
            print("Error: fifth argument must be --max_seqid when 6 arguments are given")
            sys.exit(1)
        try:
            max_seqid = float(sys.argv[6])
        except ValueError:
            print("Error: --max_seqid must be followed by a number (e.g. 90.0)")
            sys.exit(1)

    out_m8 = f"{prefix}.m8"
    out_fasta = f"{prefix}.fasta"

    df = pd.read_csv(
        input_file,
        sep='\t',
        comment='#',
        on_bad_lines='skip',
        low_memory=False
    )

    if "Length" in df.columns:
        df["Length"] = pd.to_numeric(df["Length"], errors='coerce')
    
    df = df.dropna(subset=["Length"]) if "Length" in df.columns else df

    length_mask = (df["Length"] >= min_len) & (df["Length"] <= max_len)

    seqid_mask = pd.Series(True, index=df.index)

    if max_seqid is not None:
        seqid_cols = [col for col in df.columns if col.startswith("SeqID")]
        if seqid_cols:
            for col in seqid_cols:
                df[col] = pd.to_numeric(df[col], errors='coerce')
            for col in seqid_cols:
                seqid_mask = seqid_mask & (df[col].isna() | (df[col] < max_seqid))

    df_filtered = df[length_mask & seqid_mask].copy()

    if "target" in df_filtered.columns:
        df_unique = df_filtered.drop_duplicates(subset=["target"], keep='first')
    else:
        df_unique = df_filtered

    print(f"Kept {len(df_unique)} non-redundant targets "
          f"in length range {min_len}–{max_len}", end="")
    if max_seqid is not None:
        print(f" and max seq. identity < {max_seqid}%")
    else:
        print()

    df_unique.to_csv(out_m8, sep='\t', header=False, index=False)

    with open(out_fasta, 'w') as f:
        for _, row in df_unique.iterrows():
            target_id = row["target"] if "target" in df.columns else row.iloc[0]

            if len(row) >= 2:
                sequence = str(row.iloc[1]).replace(" ", "").replace("\t", "")
            else:
                sequence = ""

            if not sequence:
                continue

            tlen = int(row["Length"]) if "Length" in row and pd.notna(row["Length"]) else len(sequence)

            header = f">{target_id}"
            print(header, file=f)

            for i in range(0, len(sequence), 60):
                print(sequence[i:i+60], file=f)

    print(f"Wrote: {out_m8}")
    print(f"Wrote: {out_fasta} (with sequences!)")

if __name__ == "__main__":
    main()