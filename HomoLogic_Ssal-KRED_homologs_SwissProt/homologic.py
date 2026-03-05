import os
import sys
import shutil
import subprocess
import glob
import csv
import re
import pandas as pd
import shutil
from pathlib import Path
import pandas as pd
import numpy as np
import mdtraj as md
from datetime import datetime
from Bio import SeqIO
from Bio.PDB import PDBParser, NeighborSearch, PDBIO, is_aa
from skbio import Protein
from Bio.Align import substitution_matrices
from pymol import cmd
import pyKVFinder
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
import open3d as o3d
from skbio.alignment import global_pairwise_align_protein  
from skbio.alignment import local_pairwise_align_protein

class HomoLogic:
    def __init__(self, reference_fasta="reference.fasta", reference_pdb="reference.pdb", reference_ligands_pdb="reference_ligands.pdb", homologs_fasta="homologs.fasta"):
        self.reference_fasta = reference_fasta
        self.reference_pdb = reference_pdb
        self.reference_ligands_pdb = reference_ligands_pdb
        self.homologs_fasta = homologs_fasta
        
    def print_progress(self, message, start_time=None):
        """Consistent progress print with timestamp."""
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        if start_time:
            elapsed = (datetime.now() - start_time).total_seconds()
            print(f"[{timestamp}] PROGRESS: {message} (Elapsed: {elapsed:.1f}s)")
        else:
            print(f"[{timestamp}] PROGRESS: {message}")

    def print_success(self, message, start_time=None):
        """Consistent success print with timestamp."""
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        if start_time:
            elapsed = (datetime.now() - start_time).total_seconds()
            print(f"[{timestamp}] ✅ SUCCESS: {message} (Total time: {elapsed:.1f}s)")
        else:
            print(f"[{timestamp}] ✅ SUCCESS: {message}")

    def print_error(self, message):
        """Consistent error print with timestamp."""
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[{timestamp}] ❌ ERROR: {message}", file=sys.stderr)

    def print_warning(self, message):
        """Consistent warning print with timestamp."""
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[{timestamp}] ⚠️  WARNING: {message}")

    def safe_create_output_folder(self, folder_path):
        """Safely create output folder, deleting existing one with warning."""
        if os.path.exists(folder_path):
            if os.path.isdir(folder_path):
                self.print_warning(f"Output folder '{folder_path}' already exists and will be overwritten")
                shutil.rmtree(folder_path)
            else:
                self.print_warning(f"File '{folder_path}' exists and will be removed")
                os.remove(folder_path)

        os.makedirs(folder_path, exist_ok=True)
        self.print_progress(f"Output folder prepared: {folder_path}")

    def calculate_sequence_homology(self, output_folder='sequence_homology_results'):
        """
        Run MMseqs2 easy-search with seed FASTA (query) and database FASTA (target).
        Results saved in output_folder.
        """
        start_time = datetime.now()
        self.print_progress(f"Starting MMseqs2 sequence homology calculation")
        self.print_progress(f"Query: {self.reference_fasta}, Database: {self.homologs_fasta}")

        # Safely create output folder
        self.safe_create_output_folder(output_folder)

        # Create temporary directory for MMseqs2
        tmp_dir = os.path.join(output_folder, 'tmp_mmseqs')
        self.safe_create_output_folder(tmp_dir)

        result_m8 = os.path.join(output_folder, 'result.m8')
        output_csv = os.path.join(output_folder, 'sequence_homology.csv')

        try:
            # Run MMseqs2 easy-search
            self.print_progress("Running MMseqs2 easy-search...")
            cmd = ['mmseqs', 'easy-search', self.reference_fasta, self.homologs_fasta, result_m8, tmp_dir]
            subprocess.run(cmd, check=True)
            self.print_progress("MMseqs2 search completed")

            def parse_fasta_ids(fasta_file):
                ids = set()
                with open(fasta_file, 'r') as f:
                    for line in f:
                        if line.startswith('>'):
                            seq_id = line[1:].split()[0].strip()
                            ids.add(seq_id)
                return ids

            all_targets = parse_fasta_ids(self.homologs_fasta)
            self.print_progress(f"Found {len(all_targets)} target sequences in database")

            # Process results
            hit_targets = set()
            rows = []
    
            # Define column headers based on MMseqs2 documentation for m8 format
            headers = [
                'reference_id', 'homolog_id', 'sequence_identity', 'alignment_length',
                'mismatches', 'gap_openings', 'query_start', 'query_end',
                'target_start', 'target_end', 'e_value', 'bit_score'
            ]

            def parse_fasta_ids(fasta_file):
                """Parse a FASTA file to extract all sequence IDs."""
                ids = set()
                with open(fasta_file, 'r') as f:
                    for line in f:
                        if line.startswith('>'):
                            seq_id = line[1:].split()[0].strip()
                            ids.add(seq_id)
                return ids

            # Parse all target IDs from database FASTA
            all_targets = parse_fasta_ids(self.homologs_fasta)

            # Collect hit targets and rows from .m8
            hit_targets = set()
            rows = []
            with open(result_m8, 'r') as f:
                reader = csv.reader(f, delimiter='\t')
                for row in reader:
                    if len(row) == 12:  # Ensure it's a valid m8 row
                        rows.append(row)
                        hit_targets.add(row[1])  # target_id is in column 2 (0-indexed)

            # Add rows for unmatched targets
            unmatched = all_targets - hit_targets
            for target_id in unmatched:
                # Create a row with 'n.d.' for all fields except target_id
                row = ['n.d.', target_id] + ['n.d.'] * 10
                rows.append(row)

            # Write to CSV with headers
            with open(output_csv, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(headers)
                writer.writerows(rows)

            # Cleanup temporary files
            try:
                os.remove(result_m8)
                shutil.rmtree(tmp_dir)
            except OSError:
                print("Warning: Could not clean up temporary files")

            print(f"Results saved to {output_csv}")
            return output_csv

        except subprocess.CalledProcessError as e:
            self.print_error(f"MMseqs2 failed: {e}")
            raise
        except Exception as e:
            self.print_error(f"Sequence homology calculation failed: {str(e)}")
            raise

    def generate_boltz_input(self, smiles_code=False, output_folder='02_boltz_input'):
        """
        Generates separate FASTA files for each protein sequence from the input multi-FASTA file.
        - smiles_code=False           → no ligands added
        - smiles_code="SMILES"        → use default ligands (from self.default_smiles or similar)
        - smiles_code=list[str]       → use exactly these SMILES strings (max 4)
        """
        from datetime import datetime
        import os

        start_time = datetime.now()
        self.print_progress(f"Starting FASTA generation from {self.homologs_fasta}")

        self.safe_create_output_folder(output_folder)

        # Parse multi-FASTA
        records = []
        current_id = None
        current_seq = []

        try:
            with open(self.homologs_fasta, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        if current_id:
                            records.append((current_id, ''.join(current_seq)))
                        current_id = line[1:].split()[0]
                        current_seq = []
                    else:
                        current_seq.append(line)

            if current_id:
                records.append((current_id, ''.join(current_seq)))

            self.print_progress(f"Parsed {len(records)} protein sequences")

            # Decide whether and which ligands to add
            ligands = []
            if smiles_code == "SMILES":
                ligands = getattr(self, 'default_smiles', []) 
            elif smiles_code is not False and smiles_code is not None:
                if isinstance(smiles_code, str):
                    ligands = [smiles_code.strip()]
                elif isinstance(smiles_code, (list, tuple)):
                    ligands = [str(s).strip() for s in smiles_code if str(s).strip()]
                else:
                    ligands = [str(smiles_code).strip()]

            # Limit to 4 ligands max
            ligands = ligands[:4]

            # Generate files
            for i, (prot_id, sequence) in enumerate(records, 1):
                filename = os.path.join(output_folder, f"{prot_id}.fasta")

                with open(filename, 'w') as fasta_file:
                    fasta_file.write(">A|protein\n")
                    fasta_file.write(f"{sequence}\n")

                    if ligands:
                        chain_letters = "BCDE"
                        for j, smiles in enumerate(ligands):
                            if j >= 4:
                                break
                            chain = chain_letters[j]
                            fasta_file.write(f">{chain}|smiles\n")
                            fasta_file.write(f"{smiles}\n")

                if i % 10 == 0 or i == len(records):
                    self.print_progress(f"Generated {i}/{len(records)} FASTA files")

            self.print_success(
                f"Generated {len(records)} FASTA files in {output_folder} "
                f"(ligands added: {bool(ligands)})",
                start_time
            )

        except Exception as e:
            self.print_error(f"Failed to generate FASTA files: {str(e)}")
            raise

    def perform_boltz_modeling(self, input_folder = "02_boltz_input", boltz_results_folder='03_boltz_results', protein_only_folder='04_structure_models'):
        """
        Performs Boltz modeling on FASTA files in input_folder and generates protein-only PDB files.
        """
        start_time = datetime.now()
        input_folder = os.path.abspath(input_folder)

        self.print_progress(f"Starting Boltz modeling from {input_folder}")

        if not os.path.exists(input_folder) or not os.path.isdir(input_folder):
            self.print_error(f"Input folder '{input_folder}' does not exist or is not a directory")
            return

        fasta_files = [f for f in os.listdir(input_folder) if f.endswith('.fasta')]
        if not fasta_files:
            self.print_progress("No FASTA files found - skipping")
            return

        self.print_progress(f"Found {len(fasta_files)} FASTA files to process")

        # Safely create results folder
        self.safe_create_output_folder(boltz_results_folder)

        # Safely create protein-only folder
        self.safe_create_output_folder(protein_only_folder)

        original_dir = os.getcwd()
        successful_modeling = 0

        try:
            for i, filename in enumerate(fasta_files, 1):
                self.print_progress(f"Processing {i}/{len(fasta_files)}: {filename}")

                try:
                    # Copy FASTA to results folder
                    shutil.copy2(os.path.join(input_folder, filename), boltz_results_folder)

                    os.chdir(boltz_results_folder)
                    command = f"boltz predict {filename} --use_msa_server"

                    subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
                    self.print_success(f"Completed modeling for {filename}", None)
                    successful_modeling += 1

                except subprocess.CalledProcessError as e:
                    self.print_error(f"Boltz modeling failed for {filename}: {e.stderr}")

                # Cleanup
                try:
                    os.remove(filename)
                except OSError:
                    pass
                
                os.chdir(original_dir)

        finally:
            os.chdir(original_dir)

        self.print_progress(f"Boltz modeling complete: {successful_modeling}/{len(fasta_files)} successful")

        # Find CIF files
        patterns = [
            os.path.join(boltz_results_folder, "**/*_model_0.cif"),
            os.path.join(boltz_results_folder, "**/*.cif"),
        ]

        cif_files = []
        for pattern in patterns:
            cif_files = glob.glob(pattern, recursive=True)
            if cif_files:
                break
            
        if not cif_files:
            self.print_error("No CIF files found after modeling")
            return

        cif_files = list(dict.fromkeys(cif_files))
        self.print_progress(f"Processing {len(cif_files)} CIF files to PDB")

        # Process CIF files to PDB (protein_only_folder already created above)
        processed = 0

        for i, cif_file in enumerate(cif_files, 1):
            if i % 5 == 0 or i == len(cif_files):
                self.print_progress(f"Converting CIF to PDB: {i}/{len(cif_files)}")       
            protein_id = os.path.basename(cif_file).split("_model_0")[0]
            output_pdb = os.path.join(protein_only_folder, f"{protein_id}.pdb")
            traj = md.load(cif_file)
            traj.save_pdb(output_pdb)
            traj = md.load(output_pdb)
            protein_atoms = traj.topology.select("protein")
            protein_traj = traj.atom_slice(protein_atoms)
            protein_traj.save_pdb(output_pdb)
            processed += 1
        self.print_success(f"Generated {processed}/{len(cif_files)} protein-only PDB files", start_time)

    def calculate_structure_homology(self, input_csv_file='01_sequence_homology_results/sequence_homology.csv', input_folder='04_structure_models', output_folder='05_structure_homology_results'):
        """
        Add RMSD and TM-score columns to MMseqs2 CSV using TM-align.
        Results saved in output_folder.
        """
        start_time = datetime.now()
        self.print_progress(f"Starting structural homology calculation")
        self.print_progress(f"Reference PDB: {self.reference_pdb}")
        self.print_progress(f"Boltz protein folder: {input_folder}")

        # Safely create output folder
        self.safe_create_output_folder(output_folder)

        if not os.path.exists(self.reference_pdb):
            self.print_error(f"Reference PDB not found: {self.reference_pdb}")
            return pd.DataFrame()

        try:
            df = pd.read_csv(input_csv_file)
            df['rmsd'] = np.nan
            df['tm_score'] = np.nan

            if not os.path.exists(self.reference_pdb):
                return df

            pdb_files = [f for f in os.listdir(input_folder) if f.endswith('.pdb')]
            pdb_dict = {Path(f).stem: os.path.join(input_folder, f) for f in pdb_files}

            for idx, row in df.iterrows():
                target_id = row.get('homolog_id', row.get('reference_id', 'n.d.'))
                if target_id == 'n.d.':
                    continue
                
                target_accession = Path(str(target_id)).stem
                pdb_path = pdb_dict.get(target_accession)
                if not pdb_path or not os.path.exists(pdb_path):
                    continue
                
                try:
                    cmd = ['TMalign', self.reference_pdb, pdb_path]
                    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)

                    # Parse RMSD and TM-score from TM-align output
                    rmsd = None
                    tm_score = None

                    for line in result.stdout.splitlines():
                        if 'RMSD=' in line and 'A' in line and '=' in line:
                            match = re.search(r'RMSD[=:]\s*([\d.]+)', line)
                            if match:
                                rmsd = float(match.group(1))

                        if "TM-score=" in line:
                            match = (re.search(r'TM-score[=:]\s*([\d.]+)', line) or 
                                    re.search(r'TM-score=\s*([\d.]+)', line))
                            if match:
                                tm_score = float(match.group(1))

                    if rmsd is not None and tm_score is not None:
                        df.at[idx, 'rmsd'] = rmsd
                        df.at[idx, 'tm_score'] = tm_score

                except Exception as e:
                    self.print_error(f"TM align failed: {str(e)}")

            output_csv = os.path.join(output_folder, 'structure_homology.csv')
            df.to_csv(output_csv, index=False)

            self.print_success(f"Structural homology evaluation saved to {output_csv}")

            return df

        except Exception as e:
            self.print_error(f"Structural homology calculation failed: {str(e)}")
            raise
    
    def superpose_structures(self, input_folder='04_structure_models', output_folder='06_superposed_structures'):
        """
        Superpose all protein PDBs in input_folder to the reference PDB using MDTraj,
        saving superposed structures to output_folder.
        """
        start_time = datetime.now()
        self.print_progress(f"Starting superposition to {self.reference_pdb}")

        self.safe_create_output_folder(output_folder)  

        try:
        
            pdb_files = [f for f in os.listdir(input_folder) if f.endswith('.pdb')]
            self.print_progress(f"Found {len(pdb_files)} PDB files to superpose")

            successful = 0
            for i, filename in enumerate(pdb_files, 1):
                input_path = os.path.join(input_folder, filename)
                output_path = os.path.join(output_folder, f"{os.path.splitext(filename)[0]}_superposed.pdb")

                try:
                    cmd.reinitialize()
                    cmd.load(self.reference_pdb, "reference")
                    cmd.load(input_path, "target")
                    cmd.align("target and name CA", "reference and name CA")
                    cmd.save(output_path, "target")
                    cmd.delete("all")
                    successful += 1
                    if i % 10 == 0 or i == len(pdb_files):
                        self.print_progress(f"Superposed {i}/{len(pdb_files)} proteins")

                except Exception as e:
                    self.print_error(f"Failed to superpose {filename}: {str(e)}")

            self.print_success(f"Superposed {successful}/{len(pdb_files)} proteins in {output_folder}", start_time)

        except Exception as e:
            self.print_error(f"Superposition failed: {str(e)}")
            raise

    class SelectBindingSiteResidues:
        """Simple selector for binding site residues (from your example)"""
        def __init__(self, binding_residues):
            self.binding_residues = binding_residues

        def accept_residue(self, residue):
            return residue in self.binding_residues

        def accept_chain(self, chain): return True
        def accept_model(self, model): return True
        def accept_atom(self, atom): return True

    def extract_reference_binding_site(self, distance_cutoff=6.0):

        start_time = datetime.now()
        self.print_progress(f"Starting binding site extraction with cutoff {distance_cutoff}Å using {self.reference_ligands_pdb}")

        parser = PDBParser(QUIET=True)

        try:
            # Load ligand structure once
            ligand_structure = parser.get_structure("ligand", self.reference_ligands_pdb)
            ligand_atoms = [atom for atom in ligand_structure.get_atoms()]
            output_path = os.path.join(f"{os.path.splitext(self.reference_pdb)[0]}_bindingsite.pdb")

            # Load protein structure
            protein_structure = parser.get_structure("protein", self.reference_pdb)
            protein_atoms = [atom for atom in protein_structure.get_atoms() if is_aa(atom.parent)]

            # Find neighbors within cutoff
            ns = NeighborSearch(protein_atoms)
            binding_residues = set()

            for ligand_atom in ligand_atoms:
                nearby_atoms = ns.search(ligand_atom.coord, distance_cutoff)
                for atom in nearby_atoms:
                    binding_residues.add(atom.get_parent())

            io = PDBIO()
            io.set_structure(protein_structure)
            io.save(output_path, self.SelectBindingSiteResidues(binding_residues))

            self.print_success(f"Extracted binding sites for reference {self.reference_pdb}", start_time)
 
        except Exception as e:
            self.print_error(f"Failed to extract binding site for reference {self.reference_pdb}: {str(e)}")

    def extract_homolog_binding_sites(self, distance_cutoff=6.0, input_folder='06_superposed_structures', output_folder='07_homolog_bindingsites'):

        start_time = datetime.now()
        self.print_progress(f"Starting binding site extraction with cutoff {distance_cutoff}Å using {self.reference_ligands_pdb}")

        self.safe_create_output_folder(output_folder)

        parser = PDBParser(QUIET=True)

        try:
            # Load ligand structure once
            ligand_structure = parser.get_structure("ligand", self.reference_ligands_pdb)
            ligand_atoms = [atom for atom in ligand_structure.get_atoms()]

            pdb_files = [f for f in os.listdir(input_folder) if f.endswith('.pdb')]
            self.print_progress(f"Found {len(pdb_files)} PDB files for binding site extraction")

            successful = 0
            for i, filename in enumerate(pdb_files, 1):
                input_path = os.path.join(input_folder, filename)
                output_path = os.path.join(output_folder, f"{os.path.splitext(filename)[0]}_bindingsite.pdb")

                try:
                    # Load protein structure
                    protein_structure = parser.get_structure("protein", input_path)
                    protein_atoms = [atom for atom in protein_structure.get_atoms() if is_aa(atom.parent)]

                    # Find neighbors within cutoff
                    ns = NeighborSearch(protein_atoms)
                    binding_residues = set()

                    for ligand_atom in ligand_atoms:
                        nearby_atoms = ns.search(ligand_atom.coord, distance_cutoff)
                        for atom in nearby_atoms:
                            binding_residues.add(atom.get_parent())

                    if not binding_residues:
                        self.print_warning(f"No binding residues found for {filename}")
                        continue

                    io = PDBIO()
                    io.set_structure(protein_structure)
                    io.save(output_path, self.SelectBindingSiteResidues(binding_residues))

                    successful += 1
                    if i % 10 == 0 or i == len(pdb_files):
                        self.print_progress(f"Extracted binding sites for {i}/{len(pdb_files)} proteins")

                except Exception as e:
                    self.print_error(f"Failed to extract binding site for {filename}: {str(e)}")

            self.print_success(f"Extracted binding sites for {successful}/{len(pdb_files)} homologs in {output_folder}", start_time)

        except Exception as e:
            self.print_error(f"Homologs binding site extraction failed: {str(e)}")
            raise

    def superpose_binding_sites(self, csv_file='05_structure_homology_results/structure_homology.csv', input_folder='07_homolog_bindingsites', output_folder='08_homolog_bindingsites_superposed', min_residues_total=5, min_aligned_atoms=5):
        """
        Refine superposition of binding sites from homologs onto the reference binding site.
        Saves each superposed target binding site, records RMSD, number of residues,
        number of aligned atoms and alignment fraction.

        Parameters:
            min_residues_total: Minimum number of CA atoms required in target binding site (default: 10)
            min_aligned_atoms:  Minimum number of atoms that must be aligned to trust the RMSD (default: 10)
        """
        start_time = datetime.now()
        self.print_progress("Starting binding site refinement superposition", start_time)

        self.safe_create_output_folder(output_folder)

        try:
            df = pd.read_csv(csv_file)
        except Exception as e:
            self.print_error(f"Cannot read input CSV {csv_file}: {str(e)}")
            return

        # Initialize necessary columns
        for col in ['rmsd_bindingsite', 'n_bindingsite_residues', 'n_aligned_residues', 'aligned_fraction']:
            if col not in df.columns:
                df[col] = np.nan

        reference_bindingsite_pdb = f"{os.path.splitext(self.reference_pdb)[0]}_bindingsite.pdb"
        ref_stem = Path(self.reference_pdb).stem

        if not os.path.exists(reference_bindingsite_pdb):
            self.print_error(f"Reference binding site not found: {reference_bindingsite_pdb}")
            return

        processed = 0
        valid_rmsd = 0

        for idx, row in df.iterrows():
            target_id = row.get('homolog_id', row.get('reference_id', 'n.d.'))
            if target_id == 'n.d.':
                continue

            target_accession = Path(str(target_id)).stem
            pdb_path = os.path.join(input_folder, f"{target_accession}_superposed_bindingsite.pdb")

            if not os.path.exists(pdb_path):
                continue

            processed += 1

            if target_accession == ref_stem:
                df.at[idx, 'rmsd_bindingsite'] = 0.0
                df.at[idx, 'n_bindingsite_residues'] = np.nan
                df.at[idx, 'n_aligned_residues'] = np.nan
                df.at[idx, 'aligned_fraction'] = np.nan
                continue

            try:
                cmd.reinitialize()
                cmd.load(reference_bindingsite_pdb, "ref_bs")
                cmd.load(pdb_path, "tgt_bs")

                cmd.select("ref_ca", "name CA and ref_bs")
                cmd.select("tgt_ca", "name CA and tgt_bs")

                n_ref = cmd.count_atoms("ref_ca")
                n_tgt = cmd.count_atoms("tgt_ca")

                df.at[idx, 'n_bindingsite_residues'] = n_tgt

                if n_tgt < min_residues_total:
                    df.at[idx, 'rmsd_bindingsite'] = np.nan
                    df.at[idx, 'n_aligned_residues'] = 0
                    df.at[idx, 'aligned_fraction'] = 0.0
                    cmd.delete("all")
                    continue

                rmsd_result = cmd.super("tgt_ca", "ref_ca")

                if isinstance(rmsd_result, tuple) and len(rmsd_result) >= 2:
                    rmsd = rmsd_result[0] if len(rmsd_result) >= 1 else np.nan
                    n_aligned = rmsd_result[1] if len(rmsd_result) >= 2 else 0
                else:
                    rmsd = np.nan
                    n_aligned = 0

                df.at[idx, 'n_aligned_residues'] = n_aligned

                if n_tgt > 0:
                    df.at[idx, 'aligned_fraction'] = n_aligned / n_tgt
                else:
                    df.at[idx, 'aligned_fraction'] = np.nan

                # Quality checks and warnings
                if n_aligned < min_aligned_atoms:
                    rmsd = np.nan
                elif rmsd is not np.nan and rmsd < 0.05 and n_aligned >= min_aligned_atoms:
                    self.print_warning(
                        f"{target_accession}: Very low RMSD ({rmsd:.3f} Å) over {n_aligned} atoms – "
                        "check if inputs are identical or consider increasing the distance cutoff "
                        "during binding site definition"
                    )

                if df.at[idx, 'aligned_fraction'] < 0.3 and n_tgt >= min_residues_total:
                    self.print_warning(
                        f"{target_accession}: Low alignment fraction ({df.at[idx, 'aligned_fraction']:.2f}) "
                        f"despite {n_tgt} residues – consider increasing the distance cutoff "
                        "during binding site definition"
                    )

                output_filename = f"{target_accession}_superposed_bindingsite_refined.pdb"
                output_path = os.path.join(output_folder, output_filename)
                cmd.save(output_path, "tgt_bs")

                df.at[idx, 'rmsd_bindingsite'] = rmsd

                if not np.isnan(rmsd):
                    valid_rmsd += 1

                cmd.delete("all")

            except Exception as e:
                self.print_error(f"Superposition failed for {target_accession}: {str(e)}")
                df.at[idx, 'rmsd_bindingsite'] = np.nan
                df.at[idx, 'n_aligned_residues'] = np.nan
                df.at[idx, 'aligned_fraction'] = np.nan
        
        # Final save
        df = df.fillna('n.d.')
        output_csv = os.path.join(output_folder, 'bindingsite_geometry_homology.csv')
        df.to_csv(output_csv, index=False)

        self.print_success(
            f"Binding site refinement completed → {processed} processed, "
            f"{valid_rmsd} Bindingsite RMSDs saved to: {output_csv} - increase distance cutoff if many 'n.d.' RMSDs",
            start_time
        )

    def analyze_binding_sites(self,
                              csv_file='08_homolog_bindingsites_superposed/structure_homology.csv',
                              input_folder='08_homolog_bindingsites_superposed',
                              bindingsite_metasequence_folder='09_bindingsite_metasequences',
                              bindingsite_similarity_results='10_bindingsite_similarity_results',
                              boltz_results_base='03_boltz_results',          # added for pLDDT
                              tolerated_misalignment=1,
                              substitution_matrix_name='BLOSUM62'):
        """
        Analyze binding sites by calculating metasequence identity.
        """
        reference_bindingsite_pdb = f"{os.path.splitext(self.reference_pdb)[0]}_bindingsite.pdb"
        if not os.path.exists(reference_bindingsite_pdb):
            self.print_error(f"Reference binding site PDB not found: {reference_bindingsite_pdb}")
            return

        aa_mapping = {
            'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
            'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
            'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
            'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
        }

        def convert_to_one_letter(input_dict):
            sequence = ""
            for key, value in input_dict.items():
                three_letter_code = str(value)[:3]
                if three_letter_code in aa_mapping:
                    sequence += aa_mapping[three_letter_code]
            return sequence

        def extract_metasequence(reference_pdb, target_pdb, tolerated_misalignment):
            reference = md.load(reference_pdb)
            target = md.load(target_pdb)
            target_residues = {i: res for i, res in enumerate(target.topology.residues)}
            target_seq = convert_to_one_letter(target_residues)
            original_length = len(target_seq)

            reference_ca_indices = [atom.index for atom in reference.topology.atoms if atom.name == 'CA']
            target_ca_indices = [atom.index for atom in target.topology.atoms if atom.name == 'CA']

            reference_ca_positions = reference.xyz[0][reference_ca_indices, :].reshape(-1, 3)
            target_ca_positions = target.xyz[0][target_ca_indices, :].reshape(-1, 3)

            distances = np.linalg.norm(reference_ca_positions[:, np.newaxis] - target_ca_positions[np.newaxis, :], axis=2)
            closest_target_indices = np.argmin(distances, axis=1)
            min_distances = distances[np.arange(len(reference_ca_indices)), closest_target_indices]

            reference_residues_with_partner = set()
            target_residues_with_partner = set()
            for ref_idx, min_distance, tgt_idx in zip(range(len(reference_ca_indices)), min_distances, closest_target_indices):
                if min_distance < tolerated_misalignment:
                    ref_residue = reference.topology.atom(reference_ca_indices[ref_idx]).residue
                    tgt_residue = target.topology.atom(target_ca_indices[tgt_idx]).residue
                    reference_residues_with_partner.add(ref_residue)
                    target_residues_with_partner.add(tgt_residue)

            target_residues = [residue for residue in target.topology.residues]
            target_residues_without_partner = [residue for residue in target_residues if residue not in target_residues_with_partner]
            atoms_to_keep_target = [atom.index for atom in target.topology.atoms if atom.residue not in target_residues_without_partner]
            target = target.atom_slice(atoms_to_keep_target)

            reference_residues = {i: res for i, res in enumerate(reference.topology.residues)}
            target_residues = {i: res for i, res in enumerate(target.topology.residues)}

            reference_seq = convert_to_one_letter(reference_residues)
            target_seq = convert_to_one_letter(target_residues)

            non_equivalent_residues = original_length - len(target_seq)
            return reference_seq, target_seq, non_equivalent_residues

        def write_fasta_file(sequence, output_file, seq_id):
            with open(output_file, 'w') as f:
                f.write(f'>{seq_id}\n{sequence}\n')

        def calculate_sequence_identity(seq1, seq2, matrix_name):
            try:
                seq1_obj = Protein(seq1)
                seq2_obj = Protein(seq2)
                matrix = substitution_matrices.load(matrix_name)
                alignment = global_pairwise_align_protein(seq1_obj, seq2_obj, substitution_matrix=matrix)
                score = alignment[1]
                return score
            except Exception as e:
                self.print_error(f"Alignment score calculation failed ({matrix_name}): {str(e)}")
                return np.nan

        def get_plddt_statistics(plddt_values):
            if len(plddt_values) == 0:
                return {k: np.nan for k in [
                    'n_res_plddt', 'median_plddt', 'avg_plddt', 'std_plddt',
                    'iqr_plddt', 'p10_plddt'
                ]}

            p = np.asarray(plddt_values, dtype=float)

            return {
                'n_res_plddt': len(p),
                'median_plddt': np.median(p),
                'avg_plddt': np.mean(p),
                'std_plddt': np.std(p),
                'iqr_plddt': np.percentile(p, 75) - np.percentile(p, 25),
                'p10_plddt': np.percentile(p, 10),
            }

        def get_plddt_stats(target_id, mode='full'):
            try:
                protein_id = os.path.splitext(target_id)[0] if '.' in str(target_id) else str(target_id)
                plddt_path = os.path.join(
                    boltz_results_base,
                    f'boltz_results_{protein_id}',
                    'predictions',
                    protein_id,
                    f'plddt_{protein_id}_model_0.npz'
                )
                if not os.path.exists(plddt_path):
                    print(f"pLDDT file not found for {protein_id}")
                    return None

                plddt_data = np.load(plddt_path)
                plddt_full = plddt_data["plddt"]

                if mode == 'full':
                    return get_plddt_statistics(plddt_full)

                # bindingsite mode
                bs_pdb = os.path.join(input_folder, f"{protein_id}_superposed_bindingsite_refined.pdb")
                if not os.path.exists(bs_pdb):
                    print(f"Binding site PDB not found for {protein_id}")
                    return None

                parser = PDBParser(QUIET=True)
                structure = parser.get_structure("bs", bs_pdb)

                bs_plddt = []
                bs_res_nums = []
                out_of_range_res = []
                for model in structure:
                    for chain in model:
                        for residue in chain:
                            if not is_aa(residue):
                                continue
                            res_num = residue.id[1]
                            idx = res_num - 1
                            bs_res_nums.append(res_num)
                            if 0 <= idx < len(plddt_full):
                                bs_plddt.append(plddt_full[idx])
                            else:
                                out_of_range_res.append(res_num)
                if out_of_range_res:
                    print(f"Out-of-range residues for {protein_id} (idx out of {len(plddt_full)}): {out_of_range_res}")

                if not bs_plddt:
                    print(f"No valid pLDDT values for binding site of {protein_id}")
                    return None
                return get_plddt_statistics(bs_plddt)

            except Exception as e:
                print(f"Exception in get_plddt_stats for {protein_id} ({mode}): {str(e)}")
                return None
        # ────────────────────────────────────────────────────────────────

        start_time = datetime.now()
        self.print_progress("Starting binding site analysis")
        self.print_progress(f"Reference binding site: {reference_bindingsite_pdb}")
        self.safe_create_output_folder(bindingsite_metasequence_folder)
        self.safe_create_output_folder(bindingsite_similarity_results)

        if not os.path.exists(csv_file):
            self.print_error(f"Input CSV file not found: {csv_file}")
            return

        try:
            df = pd.read_csv(csv_file)
            self.print_progress(f"Loaded CSV with {len(df)} entries")

            # Add pLDDT columns
            plddt_cols = [
                'n_res_plddt', 'median_plddt', 'avg_plddt', 'std_plddt', 'iqr_plddt',
                'p10_plddt'
            ]
            for col in plddt_cols:
                df[col] = np.nan                # full protein
                df[f'{col}_bs'] = np.nan        # binding site

            df['bindingsite_similarity_score'] = np.nan
            df['non_equivalent_residues'] = np.nan

            processed = 0
            for idx, row in df.iterrows():
                target_id = row.get('homolog_id', row.get('reference_id', 'n.d.'))
                if target_id == 'n.d.':
                    continue

                target = f"{input_folder}/{target_id}_superposed_bindingsite_refined.pdb"
                try:
                    common_ref_seq, common_target_seq, non_equivalent_residues = extract_metasequence(
                        reference_bindingsite_pdb, target, tolerated_misalignment
                    )

                    fasta_output = os.path.join(bindingsite_metasequence_folder, f"{target_id}.fasta")
                    write_fasta_file(common_target_seq, fasta_output, target_id)

                    alignment_score = calculate_sequence_identity(common_ref_seq, common_target_seq, substitution_matrix_name)
                    maximal_alignment_score = calculate_sequence_identity(common_ref_seq, common_ref_seq, substitution_matrix_name)

                    norm_score = alignment_score / maximal_alignment_score if maximal_alignment_score != 0 else np.nan
                    df.at[idx, 'bindingsite_similarity_score'] = norm_score
                    df.at[idx, 'non_equivalent_residues'] = non_equivalent_residues

                    # ── pLDDT for full protein ──
                    full_stats = get_plddt_stats(target_id, mode='full')
                    if full_stats:
                        for k, v in full_stats.items():
                            df.at[idx, k] = round(v, 3) if pd.notna(v) else np.nan

                    # ── pLDDT for binding site ──
                    bs_stats = get_plddt_stats(target_id, mode='bindingsite')
                    if bs_stats:
                        for k, v in bs_stats.items():
                            df.at[idx, f'{k}_bs'] = round(v, 3) if pd.notna(v) else np.nan

                    processed += 1
                    if processed % 10 == 0:
                        self.print_progress(f"Processed {processed} binding sites")

                except Exception as e:
                    self.print_error(f"Failed to analyze binding site for {target_id}: {str(e)}")

            output_csv = os.path.join(bindingsite_similarity_results, 'bindingsite_homology.csv')
            df.to_csv(output_csv, index=False)

            self.print_success(f"Generated {processed} FASTA files in {bindingsite_metasequence_folder}")
            self.print_success(f"Processed {processed} binding sites → {output_csv}", start_time)
            return df

        except Exception as e:
            self.print_error(f"Binding site analysis failed: {str(e)}")
            raise

    def analyze_reference_cavity_properties(self, input='reference_bindingsite.pdb', output_results_folder='reference_cavity_analysis', probe_out=4.0, volume_cutoff=100.0):

        start_time = datetime.now()
        self.print_progress("Starting reference cavity analysis with pyKVFinder", start_time)
    
        self.safe_create_output_folder(output_results_folder)
    
        input_path = os.path.join(input)
    
        cavity_pdb_out    = os.path.join(output_results_folder, f"reference_cavities.pdb")
        cavity_toml_out   = os.path.join(output_results_folder, f"reference_results.toml")
        cavity_barplots_out = os.path.join(output_results_folder, f"reference_barplots.pdf")
        
        try:
            results = pyKVFinder.run_workflow(
                input_path,
                probe_out=probe_out,
                volume_cutoff=volume_cutoff,
                include_depth=True,
                include_hydropathy=True,
                ignore_backbone=True
            )
    
            results.export_all(
                fn=cavity_toml_out,
                output=cavity_pdb_out,
                include_frequencies_pdf=True,
                pdf=cavity_barplots_out
            )
           
        except Exception as e:
            self.print_error(f"pyKVFinder failed for reference bindingsite")  
    
        self.print_success(
            f"Reference cavity files (PDB/TOML/PDF): {output_results_folder}\n",
            start_time
        )
    
    def analyze_cavity_properties(self,
                                  input_folder='08_homolog_bindingsites_superposed',
                                  output_cavities_folder='11_detected_cavities',
                                  output_results_folder='12_cavity_analysis_results',
                                  probe_out=4.0,
                                  volume_cutoff=100.0,
                                  csv_input='10_bindingsite_similarity_results/bindingsite_homology.csv',
                                  reference_pdb='reference_bindingsite.pdb',
                                  reference_cavity="reference_cavity_analysis/reference_cavities.pdb"):

        import os
        import numpy as np
        import pandas as pd
        import open3d as o3d
        from pathlib import Path
        from datetime import datetime

        def read_point_cloud_from_pdb(pdb_file):
            coords = []
            with open(pdb_file) as f:
                for line in f:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        coords.append([x, y, z])
            return np.array(coords)

        def run_icp_on_cavity(target_pdb, ref_pdb):
            try:
                ref_pts = read_point_cloud_from_pdb(ref_pdb)
                tgt_pts = read_point_cloud_from_pdb(target_pdb)
                if len(ref_pts) < 3 or len(tgt_pts) < 3:
                    return np.nan, 0
                pcd_ref = o3d.geometry.PointCloud()
                pcd_ref.points = o3d.utility.Vector3dVector(ref_pts)
                pcd_tgt = o3d.geometry.PointCloud()
                pcd_tgt.points = o3d.utility.Vector3dVector(tgt_pts)
                result = o3d.pipelines.registration.registration_icp(
                    pcd_tgt,
                    pcd_ref,
                    max_correspondence_distance=2.0,
                    init=np.eye(4),
                    estimation_method=o3d.pipelines.registration.TransformationEstimationPointToPoint()
                )
                rmse = result.inlier_rmse
                n_corr = len(result.correspondence_set)
                return round(float(rmse), 3), int(n_corr)
            except Exception as e:
                self.print_error(f"ICP failed for {target_pdb}: {str(e)}")
                return np.nan, 0

        output_csv_name = 'bindingsite_cavity_homology.csv'
        start_time = datetime.now()
        self.print_progress("Starting cavity analysis + ICP superposition", start_time)

        self.safe_create_output_folder(output_cavities_folder)
        self.safe_create_output_folder(output_results_folder)

        output_csv = os.path.join(output_results_folder, output_csv_name)

        if not os.path.isfile(reference_pdb):
            self.print_error(f"Reference PDB not found: {reference_pdb}")
            return pd.DataFrame()

        try:
            df = pd.read_csv(csv_input)
            self.print_progress(f"Loaded input CSV with {len(df)} rows")
        except Exception:
            self.print_error(f"Cannot read input CSV: {csv_input} → starting with empty DataFrame")
            df = pd.DataFrame()

        new_cols = [
            'cavity_volume', 'cavity_area', 'cavity_avg_depth', 'cavity_avg_hydropathy',
            'cavity_frequency_Alipathic_apolar', 'cavity_frequency_Aromatic',
            'cavity_frequency_Polar_uncharged', 'cavity_frequency_Negatively_charged',
            'cavity_frequency_Positively_charged', 'cavity_frequency_Non-standard',
            'cavity_rmsd', 'cavity_n_points'
        ]

        for col in new_cols:
            if col not in df.columns:
                df[col] = np.nan

        pdb_files = [f for f in os.listdir(input_folder) if f.endswith('.pdb')]
        if not pdb_files:
            self.print_error(f"No PDB files found in {input_folder}")
            df.fillna('n.d.').to_csv(output_csv, index=False)
            return df

        self.print_progress(f"Found {len(pdb_files)} PDB files to process")

        processed = detected = 0

        for i, filename in enumerate(pdb_files, 1):
            target_id = Path(filename).stem.split('_')[0]
            input_path = os.path.join(input_folder, filename)

            cavity_pdb_out = os.path.join(output_cavities_folder, f"{target_id}_cavities.pdb")
            cavity_toml_out = os.path.join(output_cavities_folder, f"{target_id}_results.toml")
            cavity_barplots_out = os.path.join(output_cavities_folder, f"{target_id}_barplots.pdf")

            self.print_progress(f"[{i}/{len(pdb_files)}] Analyzing → {filename}")

            try:
                results = pyKVFinder.run_workflow(
                    input_path,
                    probe_out=probe_out,
                    volume_cutoff=volume_cutoff,
                    include_depth=True,
                    include_hydropathy=True,
                    ignore_backbone=True
                )

                results.export_all(
                    fn=cavity_toml_out,
                    output=cavity_pdb_out,
                    include_frequencies_pdf=True,
                    pdf=cavity_barplots_out
                )

                volume_dict = results.volume
                if not volume_dict:
                    self.print_warning(f"No cavities detected in {filename}")
                    continue

                detected += 1

                largest_label = max(volume_dict, key=volume_dict.get)
                largest_volume = volume_dict[largest_label]

                rmsd, n_atoms = run_icp_on_cavity(
                    target_pdb=cavity_pdb_out,
                    ref_pdb=reference_cavity
                )

                mask = df['homolog_id'].isin([target_id, f"{target_id}_superposed_bindingsite_refined"])
                if not mask.any():
                    self.print_warning(f"No matching row for '{target_id}' → skipping")
                    continue

                row_idx = df.index[mask].tolist()[0]

                df.loc[row_idx, 'cavity_volume'] = largest_volume
                df.loc[row_idx, 'cavity_area'] = results.area.get(largest_label, np.nan)
                df.loc[row_idx, 'cavity_avg_depth'] = results.avg_depth.get(largest_label, np.nan)
                df.loc[row_idx, 'cavity_avg_hydropathy'] = results.avg_hydropathy.get(largest_label, np.nan)

                freq = results.frequencies.get(largest_label, {}).get("CLASS", {})
                df.loc[row_idx, 'cavity_frequency_Alipathic_apolar'] = freq.get("R1", np.nan)
                df.loc[row_idx, 'cavity_frequency_Aromatic'] = freq.get("R2", np.nan)
                df.loc[row_idx, 'cavity_frequency_Polar_uncharged'] = freq.get("R3", np.nan)
                df.loc[row_idx, 'cavity_frequency_Negatively_charged'] = freq.get("R4", np.nan)
                df.loc[row_idx, 'cavity_frequency_Positively_charged'] = freq.get("R5", np.nan)
                df.loc[row_idx, 'cavity_frequency_Non-standard'] = freq.get("RX", np.nan)

                df.loc[row_idx, 'cavity_rmsd'] = rmsd
                df.loc[row_idx, 'cavity_n_points'] = n_atoms

                processed += 1

            except Exception as e:
                self.print_error(f"Processing failed for {filename}: {str(e)}")
                continue

        df.fillna('n.d.').to_csv(output_csv, index=False)

        self.print_success(
            f"Cavity + ICP analysis finished → {processed}/{len(pdb_files)} proteins processed | {detected} cavities detected\n"
            f"• Cavity files:       {output_cavities_folder}\n"
            f"• Results table:      {output_csv}\n"
            f"• Reference used:     {reference_pdb}",
            start_time
        )

        return df

