##Repository accompanying the publication:
#"Iterative and data-driven ortholog mining enables reliable discovery of stereoselective ketoreductases"

[![Python](https://img.shields.io/badge/python-3.14%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Ketoreductases (KREDs) enable the stereoselective synthesis of chiral alcohols and are widely used in pharmaceutical manufacturing. However, reliably identifying enzymes with high activity and selectivity toward non-nativesubstrates remains challenging and often requires extensive empirical screening.

This repository contains all code, scripts, data, and analysis results for the iterative ortholog mining workflow, phylogenetic analyses, screening visualizations, statistical evaluations, and data-driven modeling presented inthe study.

**HomoLogic** — the fully automated computational pipeline for functional homology quantification — is available in a separate repository:  
→ **[HomoLogic](https://github.com/Buller-Lab/HomoLogic)**

## Overview

The study introduces an iterative, evolutionarily balanced ortholog mining strategy combined with data-driven modeling using HomoLogic:

- Random sampling of 48 distant orthologs from the seed KRED (*Ssal*-KRED) identified five promising variants (de/ee between 42.8% and >99%).
- Refined sampling of 60 orthologs in close phylogenetic proximity yielded ten superior ketoreductases (de/ee between 98% and >99%).
- MMseqs2-based homology search expanded the sequence space, and HomoLogic increased the searchable functional space ~15-fold.
- Interpretable models trained on HomoLogic-derived descriptors reliably predicted activity and stereoselectivity.

## Repository Structure & Methods Mapping

# Ssal-KRED_orthologs
# Full Experimental & Computational Datasets
only experimental.csv
full data.csv

## Dendogram & heatmap Plot Round 1
python command plot

## Dendogram & heatmap Plot Round 2
python command plot

## Swarm Plots improvement Round 1 to Round 2
python codes & commands

## HomoLogic analysis of experimentally tested orthologs
python run_workflow_orthologs.py

## Identification of Ssal-KRED homologs in SwissProt database via MMseqs2
mmseqs.yml
python get_homologs.py -i reference.fasta -o SwissProt_homologs -d /mnt/bkup/tools/SwissProt
(SwissProt can be downloaded from MMSeqs via: XYZ – replace with actual download command when available)

## HomoLogic analysis of Ssal-KRED homologs from SwissProt
python run_workflow_swissprot.py

## t-test
statistics.yml
python perform_t_test.py --orthologs /mnt/bkup/Ssal-KRED_orthologs/Ssal-KRED_orthologs_tested/12_cavity_analysis_results/bindingsite_cavity_homology.csv --swissprot /mnt/bkup/Ssal-KRED_orthologs/Ssal-KRED_homologs_SwissProt2_cavity_analysis_results/bindingsite_cavity_homology.csv --max_seqid 0.3

## Plot splitted violine distributions for sequence, structure and bindingsite distributions
COMMAND MISSING

## Identification of Ort-EZM-20_4 homologs in UniProt
python get_homologs.py -i reference.fasta -o Ort-EZM-20_4_homologs -d /mnt/bkup/tools/UniProt --bacdive-dir /mnt/bkup/Ssal-KRED_orthologs/Ort-EZM-20_4_homologs_UniProt/bacdive_taxids

## Identification of homologs of all orthologs in UniProt
python get_homologs.py -i Ssal-KRED_ortholog_library_sequences.fasta -o UniProt_homologs -d /mnt/bkup/tools/UniProt

## Lasso and multiple linear regression
 python linear_regression.py -i full_data_table_experimental_and_homologic.csv -t "2b (R)"
 python linear_regression.py -i full_data_table_experimental_and_homologic.csv -t "ee 2b (%)"

## PCA on top 3 features
python perform_clustering_2b_R.py full_data_table_experimental_and_homologic.csv --features cavity_frequency_Polar_uncharged cavity_frequency_Aromatic tm_score --color-col "2b (R)"
python perform_clustering_2b_ee.py full_data_table_experimental_and_homologic.csv --features bindingsite_similarity_score n_bindingsite_residues cavity_n_points --color-col "ee 2b (%)"

**Conda environments** (for easy reproduction):
- `homologic.yml` — HomoLogic descriptor generation and modeling
- `mmseqs2.yml` — MMseqs2 homology searches
- `treegen.yml` — Orthogroup inference and phylogenetic tree construction
- `plots.yml` — Visualization scripts (heatmaps, dendrograms, swarm plots, histograms)
- `statistics.yml` — Statistical analyses (t-tests, correlations, LASSO regression)

## System Requirements

### Operating Systems
- Linux (Ubuntu 22.04 / 24.04 recommended)
- macOS (Ventura or newer)
- Windows 10/11 (WSL2 recommended)

### Software Dependencies
- Python ≥ 3.9 (tested up to 3.14)
- External tools: MMseqs2, OrthoFinder, MAFFT, Clustal Omega, FastTree, Boltz-2 (v2.2.0)

### Hardware Requirements
- Minimum: 16 GB RAM (32 GB+ strongly recommended for OrthoFinder and large MMseqs2 searches)
- OrthoFinder originally run on 32 CPU cores (2× Xeon Gold 6142)
- No GPU required
- Non-standard hardware: None

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/YOUR_USERNAME/YOUR_REPO.git
   cd YOUR_REPO
   ```

2. Create and activate a conda environment (recommended):
   ```bash
   conda env create -f homologic.yml     # or mmseqs2.yml, treegen.yml, plots.yml, statistics.yml
   conda activate homologic
   ```

**Typical install time** on a normal desktop computer: **10–20 minutes** (Python packages and one environment).

For **HomoLogic** installation instructions, visit [https://github.com/Buller-Lab/HomoLogic](https://github.com/Buller-Lab/HomoLogic).

## Demo

A small demo dataset is included. Run:
```bash
conda activate plots
python demo.py --input data/demo/ --output results/demo/
```

**Expected Run Time for Demo** on a normal desktop (16 GB RAM): **~6–15 minutes**.

## Instructions for Use

Use the concrete commands listed in the Repository Structure section above to reproduce each analysis.

After generating candidate sequences, export them and process with the **HomoLogic** repository for full descriptor calculation and predictive modeling.

## Reproducing Results

To reproduce the full set of results from the manuscript, run the commands in the order shown in the "Repository Structure & Methods Mapping" section using the appropriate conda environments.

**Expected reproduction time** on a normal desktop (32 GB RAM):
- OrthoFinder + MMseqs2: several hours
- Remaining steps (alignments, heatmaps, statistics, LASSO): **~60–180 minutes**

## Data Availability

All datasets, intermediate files, and final analysis outputs are organized in the folders listed above.  
Supplementary data and Zenodo archive link will be added upon publication.

## License

This project is licensed under the MIT License – see the [LICENSE](LICENSE) file for details.

## Citation

If you use this code, data, or workflows, please cite:

> "Iterative and data-driven ortholog mining enables reliable discovery of stereoselective ketoreductases"  
> *Authors, Journal, Year, DOI* (insert the full citation here)

For the HomoLogic pipeline and the `get_homologs.py` script, please cite the respective repositories and the associated ligase/polymerase paper when published.
