# Ssal-KRED_orthologs

# Full Experimental & Computational Datasets
only experimental.csv
full data.csv

## Dendogram & heatmap Plot Round 1

python command plot

## Dendogram % heatmap Plot Round 2
python command plot

## Swarm Plots improvement Round 1 to Round 2

python codes & commands

## HomoLogic analysis of experimentally tested orthologs

python run_workflow_orthologs.py

## Identification of Ssal-KRED homologs in SwissProt database via MMseqs2

mmseqs.yml

python get_homologs.py -i reference.fasta -o SwissProt_homologs -d /mnt/bkup/tools/SwissProt

whereby SwissProt can be downloaded from MMSeqs via: XYZ


## HomoLogic analysis of Ssal-KRED homologs from SwissProt

python run_workflow_swissprot.py

## t-test

statistics.yml

python perform_t_test.py --orthologs /mnt/bkup/Ssal-KRED_orthologs/Ssal-KRED_orthologs_tested/12_cavity_analysis_results/bindingsite_cavity_homology.csv --swissprot /mnt/bkup/Ssal-KRED_orthologs/Ssal-KRED_homologs_SwissProt/12_cavity_analysis_results/bindingsite_cavity_homology.csv --max_seqid 0.3

## Plot splitted violine distributions for sequence, structure and bindingsite distributions

## Identification of Ort-EZM-20_4 homologs in UniProt
python get_homologs.py -i reference.fasta -o Ort-EZM-20_4_homologs -d /mnt/bkup/tools/UniProt --bacdive-dir /mnt/bkup/Ssal-KRED_orthologs/Ort-EZM-20_4_homologs_UniProt/bacdive_taxids

##Identification of homologs of all orthologs in UniProt
python get_homologs.py -i Ssal-KRED_ortholog_library_sequences.fasta -o UniProt_homologs -d /mnt/bkup/tools/UniProt
