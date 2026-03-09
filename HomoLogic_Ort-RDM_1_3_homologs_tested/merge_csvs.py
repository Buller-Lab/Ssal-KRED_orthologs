#!/usr/bin/env python3

import sys
import pandas as pd

def main(csv1_path, csv2_path, output_path):
    # Read CSV files
    df1 = pd.read_csv(csv1_path)
    df2 = pd.read_csv(csv2_path)

    # Check required column
    if "homolog_id" not in df1.columns or "homolog_id" not in df2.columns:
        raise ValueError("Both CSV files must contain a 'homolog_id' column")

    # Merge: keep all rows from csv1, add matching columns from csv2
    merged_df = df1.merge(
        df2,
        on="homolog_id",
        how="left",
        suffixes=("", "_csv2")
    )

    # Write output CSV
    merged_df.to_csv(output_path, index=False)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python merge_csvs.py <csv1.csv> <csv2.csv> <output.csv>")
        sys.exit(1)

    csv1_path = sys.argv[1]
    csv2_path = sys.argv[2]
    output_path = sys.argv[3]

    main(csv1_path, csv2_path, output_path)
