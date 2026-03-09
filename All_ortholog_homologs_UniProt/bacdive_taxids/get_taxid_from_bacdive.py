#!/usr/bin/env python3
import sys
import pandas as pd
from pathlib import Path
import argparse
import csv
from concurrent.futures import ThreadPoolExecutor, as_completed
from bacdive import BacdiveClient

def extract_ncbi_tax(strain):
    if not isinstance(strain, dict):
        return None

    def find_tax_paths(d, path=""):
        if isinstance(d, dict):
            for k, v in d.items():
                new_path = f"{path}['{k}']" if path else f"['{k}']"
                if 'tax' in k.lower():
                    if isinstance(v, (int, str)) and str(v).strip().isdigit():
                        try:
                            return int(v)
                        except (ValueError, TypeError):
                            pass
                if isinstance(v, (dict, list)):
                    found = find_tax_paths(v, new_path)
                    if found is not None:
                        return found
        elif isinstance(d, list):
            for i, item in enumerate(d):
                new_path = f"{path}[{i}]"
                found = find_tax_paths(item, new_path)
                if found is not None:
                    return found
        return None

    tax_id = find_tax_paths(strain)
    if tax_id is not None:
        return tax_id

    return None

def process_one(bid, client):
    ncbi_tax = None
    status = "not found"

    try:
        client.search(id=bid)
        gen = client.retrieve()
        strain = next(gen)
        ncbi_tax = extract_ncbi_tax(strain)
        if ncbi_tax:
            status = "found"
    except StopIteration:
        status = "no result"
    except Exception as e:
        status = f"error: {type(e).__name__}"
        print(f" Error for {bid}: {type(e).__name__} {str(e)}", flush=True)

    print(f"{bid:>8} → {status} tax_id={ncbi_tax}", flush=True)
    return bid, ncbi_tax if ncbi_tax else ""

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_csv")
    args = parser.parse_args()

    p = Path(args.input_csv)
    if not p.is_file():
        sys.exit(f"File not found: {p}")

    stem = p.stem
    out = p.with_name(f"{stem}_ncbi.csv")

    df = pd.read_csv(p, usecols=["ID"], dtype=str)
    df["ID"] = df["ID"].astype(str).str.strip()
    search_ids = df["ID"].dropna()
    search_ids = search_ids[search_ids.str.strip() != ""]
    search_ids = search_ids.drop_duplicates()

    numeric_ids = []
    for raw in search_ids:
        if raw.isdigit():
            numeric_ids.append(int(raw))

    if not numeric_ids:
        print("No valid numeric BacDive IDs found in column 'ID'.")
        return

    total = len(numeric_ids)
    print(f"Processing {total} numeric BacDive IDs (100 workers, batches of 100) ...\n")

    client = BacdiveClient()

    file_existed = out.exists() and out.stat().st_size > 0

    if not file_existed:
        with open(out, 'w', newline='', encoding='utf-8') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["tax_id", stem])

    processed_count = 0

    for batch_start in range(0, total, 100):
        batch_end = min(batch_start + 100, total)
        batch_ids = numeric_ids[batch_start:batch_end]

        print(f"\nStarting batch {batch_start+1}–{batch_end} ({len(batch_ids)} items)")

        batch_results = []

        with ThreadPoolExecutor(max_workers=100) as executor:
            future_to_bid = {
                executor.submit(process_one, bid, client): bid
                for bid in batch_ids
            }

            for future in as_completed(future_to_bid):
                _, tax_id = future.result()
                batch_results.append([tax_id, "TRUE"])
                processed_count += 1
                print(f"Progress: {processed_count}/{total}", flush=True)

        with open(out, 'a', newline='', encoding='utf-8') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerows(batch_results)

        print(f"Batch {batch_start+1}–{batch_end} saved ({len(batch_results)} rows appended)")

    print("\nCleaning output file: removing empty tax_id rows and duplicates...")

    df_out = pd.read_csv(out, dtype=str)
    df_out = df_out[df_out["tax_id"].notna() & (df_out["tax_id"] != "")]
    df_out = df_out.drop_duplicates()
    df_out.to_csv(out, index=False)

    final_rows = len(df_out)
    print(f"Cleaning complete. Final rows: {final_rows}")
    print(f"File saved → {out}")

if __name__ == "__main__":
    main()