## Repository accompanying the publication:
# Iterative and data-driven ortholog mining enables reliable discovery of stereoselective ketoreductases

[![Python](https://img.shields.io/badge/python-3.14%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/1142494427.svg)](https://doi.org/10.5281/zenodo.20756197)

Ketoreductases (KREDs) enable the stereoselective synthesis of chiral alcohols and are widely used in pharmaceutical manufacturing. However, reliably identifying enzymes with high activity and selectivity toward non-native substrates remains challenging and often requires extensive empirical screening.

This repository contains description data and code to provide figures, tables and analysis results presented in the main article. Code & data for SI figure/table elements is also provided.

**HomoLogic** — the fully automated computational pipeline for functional homology quantification — is available in a separate repository:  
→ **[HomoLogic](https://github.com/Buller-Lab/HomoLogic)**

**Conda environments** (for easy reproduction):
- `homologic.yml` — HomoLogic descriptor generation and modeling
- `mmseqs2.yml` — MMseqs2 homology searches
- `treegen.yml` — Orthogroup inference and phylogenetic tree construction
- `plots.yml` — Visualization scripts (heatmaps, dendrograms, swarm plots, histograms)
- `statistics.yml` — Statistical analyses (t-tests, correlations, LASSO regression)

## System Requirements

### Operating Systems
- Linux (Ubuntu 24.04 recommended)

### Software Dependencies
- Python ≥ 3.9 (tested with 3.14)
- External tools: MMseqs2, OrthoFinder, MAFFT, Clustal Omega, FastTree, Boltz-2 (v2.2.0), please see conda yml files or methods in publication for specific versions

### Hardware Requirements
- Minimum: 32 GB RAM (64 GB+ strongly recommended for OrthoFinder and large MMseqs2 searches)
- OrthoFinder originally run on 32 CPU cores (e.g. 2× Xeon Gold 6142)
- GPU recommended for Boltz-2 modeling (tested on NVIDIA RTX5090)

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/Buller-Lab/Ssal-KRED_orthologs.git
   cd Ssal-KRED_orthologs
   ```

2. Create and activate a conda environment (recommended):
   ```bash
   conda env create -f homologic.yml     # or mmseqs2.yml, treegen.yml, plots.yml, statistics.yml
   conda activate homologic
   ```

**Typical install time** on a normal ubuntu machine : **1-3 minutes each** (Python packages and one environment).

For **HomoLogic** installation instructions, visit [https://github.com/Buller-Lab/HomoLogic](https://github.com/Buller-Lab/HomoLogic).


## Instructions for Use

Use the concrete commands listed in the Repository Structure section above to reproduce each analysis.

After generating candidate sequences, export them and process with the **HomoLogic** repository for full descriptor calculation and predictive modeling.

## Data Availability

All datasets, intermediate files, and final analysis outputs are organized in the folders listed below.  
Supplementary data and Zenodo archive link will be provided upon publication.

## Dendogram & heatmap Plot Round 1
Data: Ssal-KRED_orthologs/Heatmaps_with_dendogram/Orthologs_Experimental_Data_Round_1.csv
```bash
conda env create -f plots.yml      # if not already done
conda activate plots
cd Heatmaps_with_dendogram
python plot_heatmap.py Orthologs_Experimental_Data_Round_1.csv
cd ..
```
Expected output: see Orthologs_Experimental_Data_Round_1_multi_heatmap_phylogenetic_order.svg
## Dendogram & heatmap Plot Round 2
Data: Ssal-KRED_orthologs/Heatmaps_with_dendogram/Orthologs_Experimental_Data_Round_2.csv
```bash
conda env create -f plots.yml     # if not already done
conda activate plots
cd Heatmaps_with_dendogram
python plot_heatmap.py Orthologs_Experimental_Data_Round_2.csv
cd ..
```
Expected output: see Orthologs_Experimental_Data_Round_2_multi_heatmap_phylogenetic_order.svg
## Swarm Plots improvement Round 1 to Round 2
Data: Ssal-KRED_orthologs/Swarm_plots_round1_vs_round2/round_1_vs_round_2.csv

Columns to plot can be selected in command line argument.
```bash
conda env create -f plots.yml     # if not already done
conda activate plots
cd Swarm_plots_round1_vs_round2
python script.py <csv_file> <main_column> <second_column> # script and columns need to be customized based on product formation to visualize
cd ..
```

## HomoLogic analysis of experimentally tested orthologs
Input data: reference_fasta="reference.fasta", reference_pdb="reference.pdb", reference_ligands_pdb="reference_ligands.pdb", homologs_fasta="Ssal-KRED_ortholog_library_sequences.fasta"
```bash
conda env create -f homologic.yml     # if not already done
conda activate homologic
cd Ssal-KRED_orthologs
python run_workflow.py
cd ..
```
Expected output: see output folders and csv files

## Identification of Ssal-KRED homologs in SwissProt database via MMseqs2
mmseqs.yml
Input data: Ssal-KRED_homologs_SwissProt/reference.fasta

Databases such as SwissProt need to be installed according to MMseqs2 documentation.
```bash
conda env create -f mmseqs2.yml     # if not already done
conda activate mmseqs2
python get_homologs.py -i reference.fasta -o SwissProt_homologs -d /mnt/bkup/tools/SwissProt # adapt database path according to where it was downloaded
cd ..
```

Expected output: see output folders and m8 & fasta files

## HomoLogic analysis of Ssal-KRED homologs from SwissProt
Input data: reference_fasta="reference.fasta", reference_pdb="reference.pdb", reference_ligands_pdb="reference_ligands.pdb", homologs_fasta="../Ssal-KRED_homologs_SwissProt/SwissProt_homologs/hits.fasta"
```bash
conda env create -f homologic.yml     # if not already done
conda activate homologic
cd HomoLogic_Ssal-KRED_homologs_SwissProt 
python run_workflow_swissprot.py
cd ..
```
Expected output: see output folders and csv files

## t-test

Input data: reference_fasta="reference.fasta", reference_pdb="reference.pdb", reference_ligands_pdb="reference_ligands.pdb", homologs_fasta="../Ssal-KRED_homologs_SwissProt/SwissProt_homologs/hits.fasta"

Additionally, analysis depents on HomoLogic analysis of orthologs and swissprot homologs

```bash
conda env create -f statistics.yml     # if not already done
conda activate statistics
cd Statistical_analysis_orthologs_vs_homologs
python perform_t_test.py --orthologs ../Ssal-KRED_orthologs/Ssal-KRED_orthologs_tested/12_cavity_analysis_results/bindingsite_cavity_homology.csv --swissprot ../Ssal-KRED_homologs_SwissProt2_cavity_analysis_results/bindingsite_cavity_homology.csv --max_seqid 0.3
cd ..
```
Expected output: see t_test_results.xlsx
## Plot splitted violine distributions for soluble vs all orthologs
Data: Ssal-KRED_orthologs/Splitted_histogram_soluble_unsoluble/soluble_expression_orthologs.csv
```bash
conda env create -f plots.yml     # if not already done
conda activate plots
cd Splitted_histogram_soluble_unsoluble
python seqid_tm_histogram.py -i soluble_expression_orthologs.csv
cd ..
```

Expected output: see split_violin_solubility.svg/pdf

## Plot splitted violine distributions for sequence, structure and bindingsite distribution
Data: Splitted_histogram_Orthologs_vs_SwissProt/far_distant_orthologs_vs_swissprot_homologs.csv
```bash
conda env create -f plots.yml     # if not already done
conda activate plots
cd Splitted_histogram_Orthologs_vs_SwissProt
python plot_seqid_tm_bindingsite_histogram.py -i far_distant_orthologs_vs_swissprot_homologs.csv
cd ..
```

Expected output: split_violins_orthologs_vs_swissprot_homologs.svg/pdf

## Identification of Ssal-KRED homologs in SwissProt database via MMseqs2
Input data: Ssal-KRED_orthologs/All_ortholog_homologs_UniProt/Ssal-KRED_ortholog_library_sequences.fasta

Databases such as UniProt need to be installed according to MMseqs2 documentation.
```bash
conda env create -f mmseqs2.yml    # if not already done
conda activate mmseqs2
python get_homologs.py -i Ssal-KRED_ortholog_library_sequences.fasta -o UniProt_homologs -d /mnt/bkup/tools/UniProt # adapt database path according to where it was downloaded
cd ..
```
Expected output: see output folders and csv files

## Lasso and multiple linear regression & PCA on top 3 features
Data: Multiple_linear_regression_and_PCA/full_data_table_experimental_and_homologic.csv
```bash
conda env create -f statistics.yml     # if not already done
conda activate statistics
cd Multiple_linear_regression_and_PCA
python linear_regression.py -i full_data_table_experimental_and_homologic.csv -t "2b (R)"
python linear_regression.py -i full_data_table_experimental_and_homologic.csv -t "ee 2b (%)"
python perform_clustering_2b_R.py full_data_table_experimental_and_homologic.csv --features cavity_frequency_Polar_uncharged cavity_frequency_Aromatic tm_score --color-col "2b (R)"
python perform_clustering_2b_ee.py full_data_table_experimental_and_homologic.csv --features bindingsite_similarity_score n_bindingsite_residues cavity_n_points --color-col "ee 2b (%)"
cd ..
```
Expected output: see output folders

## Analysis of generalization capabilities (exploiting mmseqs2 easy-cluster)
Data: Regressor_generalization_capacity/full_data_table_experimental_and_homologic.csv
```bash
conda env create -f homologic.yml     # if not already done
conda activate homologic
cd Regressor_generalization_capacity
python analyze_generalization_2b_R.py -i full_data_table_experimental_and_homologic.csv 
python analyze_generalization_2b_ee.py -i full_data_table_experimental_and_homologic.csv 
cd ..
```
Expected output: see output folders

## Investigation of factors contributing to stability
Data: Multiple_linear_regressionCorrelation_analysis_stability_factors/full_data_table_experimental_and_homologic.csv
```bash
conda env create -f homologic.yml     # if not already done
conda activate homologic
cd Regressor_generalization_capacity
python calculate_stability_descriptors.py ../HomoLogic_Ssal-KRED_orthologs_tested/04_structure_models stability_features.csv 
python analyze_correlations_with_melting_temperatures.py full_data_table_experimental_and_homologic.csv stability_features.csv
cd ..
```
Expected output: see output folders

## License

This project is licensed under the MIT License – see file for details.

## Citation

If you use this code, data, or workflows, please cite:
> "Iterative and data-driven ortholog mining enables reliable discovery of stereoselective ketoreductases"  
(DOI specified after publication)

For the `get_homologs.py` script, please cite:
https://pubs.acs.org/doi/full/10.1021/acscatal.4c04474
