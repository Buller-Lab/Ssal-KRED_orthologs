from homologic import HomoLogic

# Define reference enzyme sequence, structure, and ligands, as well as homologous sequences (search space)
hl = HomoLogic(reference_fasta="reference.fasta", reference_pdb="reference.pdb", reference_ligands_pdb="reference_ligands.pdb", homologs_fasta="../Ssal-KRED_homologs_SwissProt/SwissProt_homologs/hits.fasta")

# Calculate sequence homology between reference and homologs via MMseqs2
hl.calculate_sequence_homology(output_folder='01_sequence_homology_results')

# Generate input files for Boltz-2 cofolding
hl.generate_boltz_input(smiles_code=["C1C=CN(C=C1C(=O)N)[C@H]2[C@@H]([C@@H]([C@H](O2)COP(=O)(O)OP(=O)(O)OC[C@@H]3[C@H]([C@H]([C@@H](O3)N4C=NC5=C(N=CN=C54)N)OP(=O)(O)O)O)O)O"], output_folder='02_boltz_input')

# Perform Boltz-2 cofolding and generate protein only structure models
hl.perform_boltz_modeling(input_folder = "02_boltz_input", boltz_results_folder='03_boltz_results', protein_only_folder='04_structure_models')

# Calculate structure homology between reference and homologs via TM-align
hl.calculate_structure_homology(input_csv_file='01_sequence_homology_results/sequence_homology.csv', input_folder='04_structure_models', output_folder='05_structure_homology_results')

# Superpose all structure models onto the reference structure via PyMol for consistent spatial orientation
hl.superpose_structures(input_folder='04_structure_models', output_folder='06_superposed_structures')

# Based on given distance from reference ligand, extract binding sites from reference structure
hl.extract_reference_binding_site(distance_cutoff=6.0)


# Based on given distance from reference ligand, extract binding sites from superposed homolog structures
hl.extract_homolog_binding_sites(distance_cutoff=6.0, input_folder='06_superposed_structures', output_folder='07_homolog_bindingsites')

# Superpose all homolog binding sites onto the reference binding site via PyMol
hl.superpose_binding_sites(csv_file='05_structure_homology_results/structure_homology.csv', input_folder='07_homolog_bindingsites', output_folder='08_homolog_bindingsites_superposed')

# Analyze local "metasequences" of binding sites 
hl.analyze_binding_sites(csv_file='08_homolog_bindingsites_superposed/bindingsite_geometry_homology.csv', input_folder='08_homolog_bindingsites_superposed', bindingsite_metasequence_folder='09_bindingsite_metasequences', bindingsite_similarity_results='10_bindingsite_similarity_results', tolerated_misalignment=1.0)

# Analyze cavity properties of reference binding site
hl.analyze_reference_cavity_properties(input='reference_bindingsite.pdb', output_results_folder='reference_cavity_analysis', probe_out=4.0, volume_cutoff=100.0)

# And finally, detect cavities analyze their properties
hl.analyze_cavity_properties(input_folder='08_homolog_bindingsites_superposed', output_cavities_folder='11_detected_cavities', output_results_folder='12_cavity_analysis_results', probe_out=4.0, volume_cutoff=100.0, csv_input='10_bindingsite_similarity_results/bindingsite_homology.csv')