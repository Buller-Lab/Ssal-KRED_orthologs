import os
import sys
import shutil
import subprocess
import glob
import csv
import re
import pandas as pd
import numpy as np
import mdtraj as md
from datetime import datetime
from Bio.PDB import PDBParser, NeighborSearch, PDBIO, is_aa
from Bio.Align import substitution_matrices
from pathlib import Path
import pyKVFinder
import open3d as o3d
from skbio.alignment import global_pairwise_align_protein
from skbio import Protein, io

class HomoLogic:
    def __init__(self, reference_fasta="reference.fasta", reference_pdb="reference.pdb",
                 reference_ligands_pdb="reference_ligands.pdb", homologs_fasta="homologs.fasta"):
        self.reference_fasta = reference_fasta
        self.reference_pdb = reference_pdb
        self.reference_ligands_pdb = reference_ligands_pdb
        self.homologs_fasta = homologs_fasta

    def print_progress(self, message, start_time=None):
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        if start_time:
            elapsed = (datetime.now() - start_time).total_seconds()
            print(f"[{timestamp}] {message} ({elapsed:.1f}s)")
        else:
            print(f"[{timestamp}] {message}")

    def print_success(self, message, start_time=None):
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        if start_time:
            elapsed = (datetime.now() - start_time).total_seconds()
            print(f"[{timestamp}] SUCCESS: {message} ({elapsed:.1f}s)")
        else:
            print(f"[{timestamp}] SUCCESS: {message}")

    def print_error(self, message):
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[{timestamp}] ERROR: {message}", file=sys.stderr)

    def print_warning(self, message):
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[{timestamp}] WARNING: {message}")

    def safe_create_output_folder(self, folder_path):
        if os.path.exists(folder_path):
            if os.path.isdir(folder_path):
                self.print_warning(f"Output folder '{folder_path}' exists → overwriting")
                shutil.rmtree(folder_path)
            else:
                self.print_warning(f"File '{folder_path}' exists → removing")
                os.remove(folder_path)
        os.makedirs(folder_path, exist_ok=True)
        self.print_progress(f"Output folder prepared: {folder_path}")

    @staticmethod
    def compute_identity_gaps_mismatches_openings(ref_aln_str, tgt_aln_str):
        matches = 0
        aln_len = 0
        for a, b in zip(ref_aln_str, tgt_aln_str):
            if a != '-' and b != '-':
                aln_len += 1
                if a == b:
                    matches += 1
        identity = matches / aln_len if aln_len > 0 else 0.0
        mismatches = aln_len - matches
        opens = 0
        in_gap = False
        for a, b in zip(ref_aln_str, tgt_aln_str):
            is_gap = (a == '-' or b == '-')
            if is_gap and not in_gap:
                opens += 1
                in_gap = True
            elif not is_gap:
                in_gap = False
        return identity, mismatches, opens, aln_len

    def align_and_compute_identity_stats(self, ref_seq, tgt_seq, matrix_name='BLOSUM62',
                                         gap_open=11, gap_extend=1, return_alignment=False):
        alignment, _, _ = global_pairwise_align_protein(
            ref_seq, tgt_seq,
            gap_open_penalty=gap_open,
            gap_extend_penalty=gap_extend,
            substitution_matrix=substitution_matrices.load(matrix_name),
            penalize_terminal_gaps=False
        )
        if len(alignment) != 2:
            return None
        ref_aln = str(alignment[0])
        tgt_aln = str(alignment[1])
        identity, mismatches, gap_openings, aln_len = self.compute_identity_gaps_mismatches_openings(ref_aln, tgt_aln)
        result = {
            'identity': identity,
            'identity_pct_str': f"{identity:.3f}",
            'mismatches': mismatches,
            'gap_openings': gap_openings,
            'alignment_length': aln_len,
            'full_alignment_length': len(ref_aln)
        }
        if return_alignment:
            result['ref_aln'] = ref_aln
            result['tgt_aln'] = tgt_aln
        return result

    def calculate_sequence_homology(
        self,
        output_folder='sequence_homology_results',
        substitution_matrix_name='BLOSUM62',
        seqid_beyond_mmseqs=True
    ):
        start_time = datetime.now()
        self.print_progress("Starting sequence homology calculation", start_time)

        self.safe_create_output_folder(output_folder)
        tmp_dir = os.path.join(output_folder, 'tmp_mmseqs')
        self.safe_create_output_folder(tmp_dir)

        result_m8 = os.path.join(output_folder, 'result.m8')
        output_csv = os.path.join(output_folder, 'sequence_homology.csv')

        cmd = ['mmseqs', 'easy-search', '--max-seqs', '100000',
               self.reference_fasta, self.homologs_fasta, result_m8, tmp_dir]
        subprocess.run(cmd, check=True)

        homolog_seqs = {
            seq.metadata['id'].strip(): seq
            for seq in io.read(self.homologs_fasta, format='fasta', constructor=Protein)
        }
        self.print_progress(f"Loaded {len(homolog_seqs)} sequences")

        ref_seq_record = next(io.read(self.reference_fasta, format='fasta', constructor=Protein))
        ref_id = ref_seq_record.metadata['id'].strip()
        ref_seq = ref_seq_record

        headers = [
            'reference_id', 'homolog_id', 'sequence_identity', 'alignment_length',
            'mismatches', 'gap_openings', 'query_start', 'query_end',
            'target_start', 'target_end', 'e_value', 'bit_score'
        ]
        rows = []
        good_hits = set()
        fallback_targets = set()

        with open(result_m8, 'r') as f:
            for line in f:
                row = line.strip().split('\t')
                if len(row) != 12:
                    continue
                target_id = row[1].strip()
                try:
                    pident = float(row[2])
                    if 0 < pident <= 1:
                        good_hits.add(target_id)
                        rows.append(row)
                        continue
                except ValueError:
                    pass

                fallback_targets.add(target_id)
                rows.append(row)

        all_targets = set(homolog_seqs)
        unmatched = all_targets - good_hits - fallback_targets

        self.print_progress(
            f"MMseqs2: {len(good_hits)} good | {len(fallback_targets)} fallback | {len(unmatched)} unmatched"
        )

        if fallback_targets:
            self.print_progress(f"Global fallback for {len(fallback_targets)} sequences")
            for target_id in sorted(fallback_targets):
                if target_id not in homolog_seqs:
                    continue
                tgt_seq = homolog_seqs[target_id]
                stats = self.align_and_compute_identity_stats(ref_seq, tgt_seq, substitution_matrix_name)
                if stats:
                    for i, r in enumerate(rows):
                        if r[1].strip() == target_id:
                            rows[i] = [
                                ref_id, target_id,
                                stats['identity_pct_str'],
                                str(stats['full_alignment_length']),
                                str(stats['mismatches']),
                                str(stats['gap_openings']),
                                '1', str(len(ref_seq)),
                                '1', str(len(tgt_seq)),
                                r[10], r[11]          
                            ]
                            break

        if unmatched and seqid_beyond_mmseqs:
            self.print_progress(f"Global alignment for {len(unmatched)} unmatched sequences")
            for target_id in sorted(unmatched):
                if target_id not in homolog_seqs:
                    rows.append([ref_id, target_id] + ['n.d.'] * 10)
                    continue

                tgt_seq = homolog_seqs[target_id]
                stats = self.align_and_compute_identity_stats(ref_seq, tgt_seq, substitution_matrix_name)

                if stats:
                    rows.append([
                        ref_id, target_id,
                        stats['identity_pct_str'],
                        str(stats['full_alignment_length']),
                        str(stats['mismatches']),
                        str(stats['gap_openings']),
                        '1', str(len(ref_seq)),
                        '1', str(len(tgt_seq)),
                        'n.d.', 'n.d.'
                    ])
                else:
                    rows.append([ref_id, target_id] + ['n.d.'] * 10)

        elif unmatched and not seqid_beyond_mmseqs:
            self.print_progress(
                f"Skipping global alignment for {len(unmatched)} unmatched sequences "
                f"(seqid_beyond_mmseqs=False)"
            )
            for target_id in sorted(unmatched):
                rows.append([ref_id, target_id] + ['n.d.'] * 10)

        with open(output_csv, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(headers)
            writer.writerows(rows)

        self.print_success(f"Results saved → {output_csv} ({len(rows)} entries)", start_time)

        try:
            os.remove(result_m8)
            shutil.rmtree(tmp_dir)
        except OSError:
            pass

        return output_csv

    def generate_boltz_input(self, smiles_code=False, output_folder='02_boltz_input'):
        start_time = datetime.now()
        self.print_progress("Generating per-protein FASTA files", start_time)
        self.safe_create_output_folder(output_folder)

        records = []
        current_id = None
        current_seq = []
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

        self.print_progress(f"Parsed {len(records)} sequences")

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
        ligands = ligands[:4]

        for i, (prot_id, seq) in enumerate(records, 1):
            path = os.path.join(output_folder, f"{prot_id}.fasta")
            with open(path, 'w') as f:
                f.write(">A|protein\n" + seq + "\n")
                if ligands:
                    for j, smi in enumerate(ligands):
                        if j >= 4: break
                        f.write(f">{chr(66 + j)}|smiles\n{smi}\n")
            if i % 20 == 0 or i == len(records):
                self.print_progress(f"Generated {i}/{len(records)} FASTA files")

        self.print_success(f"Generated {len(records)} FASTA files (ligands: {bool(ligands)})", start_time)

    def perform_boltz_modeling(self, input_folder="02_boltz_input", boltz_results_folder='03_boltz_results',
                               protein_only_folder='04_structure_models'):
        start_time = datetime.now()
        input_folder = os.path.abspath(input_folder)
        self.print_progress(f"Starting Boltz modeling: {input_folder}", start_time)

        if not os.path.isdir(input_folder):
            self.print_error(f"Input folder missing: {input_folder}")
            return

        fasta_files = [f for f in os.listdir(input_folder) if f.endswith('.fasta')]
        if not fasta_files:
            self.print_progress("No FASTA files → skipping")
            return

        self.print_progress(f"Found {len(fasta_files)} FASTA files")
        self.safe_create_output_folder(boltz_results_folder)
        self.safe_create_output_folder(protein_only_folder)

        original_dir = os.getcwd()
        successful = 0

        for i, fn in enumerate(fasta_files, 1):
            if i % 10 == 0 or i == len(fasta_files):
                self.print_progress(f"Processing {i}/{len(fasta_files)}")
            src = os.path.join(input_folder, fn)
            dst = os.path.join(boltz_results_folder, fn)
            shutil.copy2(src, dst)
            os.chdir(boltz_results_folder)
            try:
                subprocess.run(f"boltz predict {fn} --use_msa_server", shell=True, check=True, capture_output=True, text=True)
                successful += 1
            except subprocess.CalledProcessError as e:
                self.print_error(f"Boltz failed {fn} (code {e.returncode})")
            finally:
                try:
                    os.remove(fn)
                except OSError:
                    pass
                os.chdir(original_dir)

        self.print_progress(f"Boltz: {successful}/{len(fasta_files)} successful")

        cif_files = glob.glob(os.path.join(boltz_results_folder, "**/*_model_0.cif"), recursive=True) or \
                    glob.glob(os.path.join(boltz_results_folder, "**/*.cif"), recursive=True)
        cif_files = list(dict.fromkeys(cif_files))

        if not cif_files:
            self.print_error("No CIF files found after modeling")
            return

        self.print_progress(f"Found {len(cif_files)} CIF files")

        processed = 0
        for i, cif in enumerate(cif_files, 1):
            if i % 10 == 0 or i == len(cif_files):
                self.print_progress(f"Converting {i}/{len(cif_files)}")
            pid = os.path.basename(cif).split("_model_0")[0]
            out_pdb = os.path.join(protein_only_folder, f"{pid}.pdb")
            traj = md.load(cif)
            traj.save_pdb(out_pdb)
            traj = md.load(out_pdb)
            sel = traj.topology.select("protein")
            protein_traj = traj.atom_slice(sel)
            protein_traj.save_pdb(out_pdb)
            processed += 1

        self.print_success(f"Generated {processed}/{len(cif_files)} protein-only PDBs", start_time)

    def calculate_structure_homology(self, input_csv_file='01_sequence_homology_results/sequence_homology.csv',
                                     input_folder='04_structure_models', output_folder='05_structure_homology_results'):
        start_time = datetime.now()
        self.print_progress("Starting structural homology (TM-align)", start_time)
        self.safe_create_output_folder(output_folder)

        if not os.path.exists(self.reference_pdb):
            self.print_error(f"Reference PDB missing: {self.reference_pdb}")
            return pd.DataFrame()

        df = pd.read_csv(input_csv_file)
        df['rmsd'] = np.nan
        df['tm_score'] = np.nan

        pdb_files = [f for f in os.listdir(input_folder) if f.endswith('.pdb')]
        pdb_dict = {Path(f).stem: os.path.join(input_folder, f) for f in pdb_files}

        for idx, row in df.iterrows():
            tid = row.get('homolog_id', row.get('reference_id', 'n.d.'))
            if tid == 'n.d.':
                continue
            acc = Path(str(tid)).stem
            pdb_path = pdb_dict.get(acc)
            if not pdb_path or not os.path.exists(pdb_path):
                continue
            try:
                res = subprocess.run(['TMalign', self.reference_pdb, pdb_path],
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
                rmsd = tm = None
                for line in res.stdout.splitlines():
                    if 'RMSD=' in line:
                        m = re.search(r'RMSD[=:]\s*([\d.]+)', line)
                        if m: rmsd = float(m.group(1))
                    if 'TM-score=' in line:
                        m = re.search(r'TM-score[=:]\s*([\d.]+)', line) or re.search(r'TM-score=\s*([\d.]+)', line)
                        if m: tm = float(m.group(1))
                if rmsd is not None and tm is not None:
                    df.at[idx, 'rmsd'] = rmsd
                    df.at[idx, 'tm_score'] = tm
            except subprocess.CalledProcessError:
                continue

        out_csv = os.path.join(output_folder, 'structure_homology.csv')
        df.to_csv(out_csv, index=False)
        self.print_success(f"Results saved → {out_csv}", start_time)
        return df

    def superpose_structures(self, input_folder='04_structure_models', output_folder='06_superposed_structures'):
        start_time = datetime.now()
        self.print_progress("Starting superposition to reference", start_time)
        self.safe_create_output_folder(output_folder)

        pdb_files = [f for f in os.listdir(input_folder) if f.endswith('.pdb')]
        self.print_progress(f"Found {len(pdb_files)} PDB files")

        from pymol import cmd
        successful = 0

        for i, fn in enumerate(pdb_files, 1):
            in_path = os.path.join(input_folder, fn)
            out_path = os.path.join(output_folder, f"{os.path.splitext(fn)[0]}_superposed.pdb")
            try:
                cmd.reinitialize()
                cmd.load(self.reference_pdb, "reference")
                cmd.load(in_path, "target")
                cmd.align("target and name CA", "reference and name CA")
                cmd.save(out_path, "target")
                cmd.delete("all")
                successful += 1
                if i % 10 == 0 or i == len(pdb_files):
                    self.print_progress(f"Superposed {i}/{len(pdb_files)}")
            except Exception as e:
                self.print_warning(f"Superposition failed {fn}: {type(e).__name__}")

        self.print_success(f"Superposed {successful}/{len(pdb_files)} structures", start_time)

    class SelectBindingSiteResidues:
        def __init__(self, binding_residues):
            self.binding_residues = binding_residues
        def accept_residue(self, residue):
            return residue in self.binding_residues
        def accept_chain(self, chain):
            return True
        def accept_model(self, model):
            return True
        def accept_atom(self, atom):
            return True

    def extract_reference_binding_site(self, distance_cutoff=6.0):
        start_time = datetime.now()
        self.print_progress(f"Extracting reference binding site ({distance_cutoff}Å)", start_time)
        parser = PDBParser(QUIET=True)
        ligand_struct = parser.get_structure("ligand", self.reference_ligands_pdb)
        ligand_atoms = list(ligand_struct.get_atoms())
        out_path = f"{os.path.splitext(self.reference_pdb)[0]}_bindingsite.pdb"
        prot_struct = parser.get_structure("protein", self.reference_pdb)
        prot_atoms = [a for a in prot_struct.get_atoms() if is_aa(a.parent)]
        ns = NeighborSearch(prot_atoms)
        binding_res = set()
        for la in ligand_atoms:
            for a in ns.search(la.coord, distance_cutoff):
                binding_res.add(a.get_parent())
        io = PDBIO()
        io.set_structure(prot_struct)
        io.save(out_path, self.SelectBindingSiteResidues(binding_res))
        self.print_success(f"Reference binding site saved → {out_path}", start_time)

    def extract_homolog_binding_sites(self, distance_cutoff=6.0, input_folder='06_superposed_structures',
                                      output_folder='07_homolog_bindingsites'):
        start_time = datetime.now()
        self.print_progress(f"Extracting homolog binding sites ({distance_cutoff}Å)", start_time)
        self.safe_create_output_folder(output_folder)
        parser = PDBParser(QUIET=True)
        ligand_struct = parser.get_structure("ligand", self.reference_ligands_pdb)
        ligand_atoms = list(ligand_struct.get_atoms())
        pdb_files = [f for f in os.listdir(input_folder) if f.endswith('.pdb')]
        self.print_progress(f"Found {len(pdb_files)} superposed structures")
        successful = 0
        for i, fn in enumerate(pdb_files, 1):
            in_path = os.path.join(input_folder, fn)
            out_path = os.path.join(output_folder, f"{os.path.splitext(fn)[0]}_bindingsite.pdb")
            try:
                prot_struct = parser.get_structure("protein", in_path)
                prot_atoms = [a for a in prot_struct.get_atoms() if is_aa(a.parent)]
                ns = NeighborSearch(prot_atoms)
                binding_res = set()
                for la in ligand_atoms:
                    for a in ns.search(la.coord, distance_cutoff):
                        binding_res.add(a.get_parent())
                if not binding_res:
                    continue
                io = PDBIO()
                io.set_structure(prot_struct)
                io.save(out_path, self.SelectBindingSiteResidues(binding_res))
                successful += 1
                if i % 10 == 0 or i == len(pdb_files):
                    self.print_progress(f"Processed {i}/{len(pdb_files)}")
            except Exception:
                continue
        self.print_success(f"Extracted {successful}/{len(pdb_files)} binding sites", start_time)

    def superpose_binding_sites(self, csv_file='05_structure_homology_results/structure_homology.csv',
                                input_folder='07_homolog_bindingsites', output_folder='08_homolog_bindingsites_superposed',
                                min_residues_total=5, min_aligned_atoms=5):
        start_time = datetime.now()
        self.print_progress("Starting refined binding site superposition", start_time)
        self.safe_create_output_folder(output_folder)

        df = pd.read_csv(csv_file)
        for col in ['rmsd_bindingsite', 'n_bindingsite_residues', 'n_aligned_residues', 'aligned_fraction']:
            if col not in df.columns:
                df[col] = np.nan

        ref_bs_pdb = f"{os.path.splitext(self.reference_pdb)[0]}_bindingsite.pdb"
        ref_stem = Path(self.reference_pdb).stem
        if not os.path.exists(ref_bs_pdb):
            self.print_error(f"Reference binding site missing: {ref_bs_pdb}")
            return

        from pymol import cmd
        processed = valid_rmsd = 0

        for idx, row in df.iterrows():
            tid = row.get('homolog_id', row.get('reference_id', 'n.d.'))
            if tid == 'n.d.':
                continue
            acc = Path(str(tid)).stem
            in_pdb = os.path.join(input_folder, f"{acc}_superposed_bindingsite.pdb")
            if not os.path.exists(in_pdb):
                continue
            processed += 1
            if acc == ref_stem:
                df.at[idx, 'rmsd_bindingsite'] = 0.0
                continue
            try:
                cmd.reinitialize()
                cmd.load(ref_bs_pdb, "ref_bs")
                cmd.load(in_pdb, "tgt_bs")
                cmd.select("ref_ca", "name CA and ref_bs")
                cmd.select("tgt_ca", "name CA and tgt_bs")
                n_tgt = cmd.count_atoms("tgt_ca")
                df.at[idx, 'n_bindingsite_residues'] = n_tgt
                if n_tgt < min_residues_total:
                    cmd.delete("all")
                    continue
                rmsd_res = cmd.super("tgt_ca", "ref_ca")
                rmsd = rmsd_res[0] if isinstance(rmsd_res, tuple) and len(rmsd_res) >= 2 else np.nan
                n_aln = rmsd_res[1] if isinstance(rmsd_res, tuple) and len(rmsd_res) >= 2 else 0
                df.at[idx, 'n_aligned_residues'] = n_aln
                df.at[idx, 'aligned_fraction'] = n_aln / n_tgt if n_tgt > 0 else np.nan
                if n_aln < min_aligned_atoms:
                    rmsd = np.nan
                out_fn = f"{acc}_superposed_bindingsite_refined.pdb"
                cmd.save(os.path.join(output_folder, out_fn), "tgt_bs")
                df.at[idx, 'rmsd_bindingsite'] = rmsd
                if not np.isnan(rmsd):
                    valid_rmsd += 1
                cmd.delete("all")
            except Exception:
                df.at[idx, 'rmsd_bindingsite'] = np.nan
                df.at[idx, 'n_aligned_residues'] = np.nan
                df.at[idx, 'aligned_fraction'] = np.nan

        df = df.fillna('n.d.')
        out_csv = os.path.join(output_folder, 'bindingsite_geometry_homology.csv')
        df.to_csv(out_csv, index=False)
        self.print_success(f"Refinement done → {processed} processed, {valid_rmsd} valid RMSDs", start_time)

    def analyze_binding_sites(self,
                              csv_file='08_homolog_bindingsites_superposed/bindingsite_geometry_homology.csv',
                              input_folder='08_homolog_bindingsites_superposed',
                              bindingsite_metasequence_folder='09_bindingsite_metasequences',
                              bindingsite_similarity_results='10_bindingsite_similarity_results',
                              boltz_results_base='03_boltz_results',
                              tolerated_misalignment=1,
                              substitution_matrix_name='BLOSUM62'):
        start_time = datetime.now()
        self.print_progress("Starting binding site metasequence & similarity analysis", start_time)

        ref_bs_pdb = f"{os.path.splitext(self.reference_pdb)[0]}_bindingsite.pdb"
        if not os.path.exists(ref_bs_pdb):
            self.print_error(f"Reference binding site missing: {ref_bs_pdb}")
            return

        aa_mapping = {
            'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L',
            'MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y'
        }

        def one_letter(res_dict):
            return ''.join(aa_mapping.get(str(r)[:3], '') for r in res_dict.values())

        def extract_metasequence(ref_pdb, tgt_pdb, max_dist):
            ref = md.load(ref_pdb)
            tgt = md.load(tgt_pdb)
            tgt_res = {i:res for i,res in enumerate(tgt.topology.residues)}
            orig_len = len(tgt_res)
            ref_ca_idx = [a.index for a in ref.topology.atoms if a.name == 'CA']
            tgt_ca_idx = [a.index for a in tgt.topology.atoms if a.name == 'CA']
            ref_xyz = ref.xyz[0][ref_ca_idx]
            tgt_xyz = tgt.xyz[0][tgt_ca_idx]
            dists = np.linalg.norm(ref_xyz[:,None] - tgt_xyz[None,:], axis=-1)
            closest = np.argmin(dists, axis=1)
            min_d = dists[np.arange(len(ref_ca_idx)), closest]
            matched_tgt = {tgt.topology.atom(tgt_ca_idx[i]).residue for i,d in enumerate(min_d) if d < max_dist}
            keep = [a.index for a in tgt.topology.atoms if a.residue in matched_tgt]
            tgt = tgt.atom_slice(keep)
            ref_res = {i:res for i,res in enumerate(ref.topology.residues)}
            tgt_res = {i:res for i,res in enumerate(tgt.topology.residues)}
            return one_letter(ref_res), one_letter(tgt_res), orig_len - len(tgt_res)

        def write_fasta(seq, path, sid):
            with open(path, 'w') as f:
                f.write(f">{sid}\n{seq}\n")

        def plddt_stats(values):
            if not values.size:
                return {k:np.nan for k in ['n_res_plddt','median_plddt','avg_plddt','std_plddt','iqr_plddt','p10_plddt']}
            return {
                'n_res_plddt': len(values),
                'median_plddt': np.median(values),
                'avg_plddt': np.mean(values),
                'std_plddt': np.std(values),
                'iqr_plddt': np.percentile(values,75) - np.percentile(values,25),
                'p10_plddt': np.percentile(values,10),
            }

        def get_plddt(target_id, mode='full'):
            pid = os.path.splitext(str(target_id))[0]
            plddt_path = os.path.join(boltz_results_base, f'boltz_results_{pid}', 'predictions', pid, f'plddt_{pid}_model_0.npz')
            if not os.path.exists(plddt_path):
                return None
            plddt = np.load(plddt_path)['plddt']
            if mode == 'full':
                return plddt_stats(plddt)
            bs_pdb = os.path.join(input_folder, f"{pid}_superposed_bindingsite_refined.pdb")
            if not os.path.exists(bs_pdb):
                return None
            struct = PDBParser(QUIET=True).get_structure("bs", bs_pdb)
            bs_plddt = []
            for model in struct:
                for chain in model:
                    for res in chain:
                        if not is_aa(res): continue
                        idx = res.id[1] - 1
                        if 0 <= idx < len(plddt):
                            bs_plddt.append(plddt[idx])
            return plddt_stats(np.array(bs_plddt)) if bs_plddt else None

        self.safe_create_output_folder(bindingsite_metasequence_folder)
        self.safe_create_output_folder(bindingsite_similarity_results)

        df = pd.read_csv(csv_file)
        plddt_cols = ['n_res_plddt','median_plddt','avg_plddt','std_plddt','iqr_plddt','p10_plddt']
        for c in plddt_cols:
            df[c] = np.nan
            df[f'{c}_bs'] = np.nan
        df['bindingsite_similarity_score'] = np.nan
        df['non_equivalent_residues'] = np.nan

        processed = 0
        for idx, row in df.iterrows():
            tid = row.get('homolog_id', row.get('reference_id', 'n.d.'))
            if tid == 'n.d.':
                continue
            tgt_pdb = os.path.join(input_folder, f"{tid}_superposed_bindingsite_refined.pdb")
            if not os.path.exists(tgt_pdb):
                continue
            try:
                ref_seq, tgt_seq, non_eq = extract_metasequence(ref_bs_pdb, tgt_pdb, tolerated_misalignment)
                write_fasta(tgt_seq, os.path.join(bindingsite_metasequence_folder, f"{tid}.fasta"), tid)
                stats = self.align_and_compute_identity_stats(Protein(ref_seq), Protein(tgt_seq), substitution_matrix_name, gap_open=10, gap_extend=1)
                if stats:
                    max_stats = self.align_and_compute_identity_stats(Protein(ref_seq), Protein(ref_seq), substitution_matrix_name, gap_open=10, gap_extend=1)
                    max_sc = max_stats['identity'] if max_stats else 1.0
                    df.at[idx, 'bindingsite_similarity_score'] = stats['identity'] / max_sc if max_sc else np.nan
                df.at[idx, 'non_equivalent_residues'] = non_eq
                full_s = get_plddt(tid, 'full')
                if full_s:
                    for k,v in full_s.items():
                        df.at[idx, k] = round(v,3) if pd.notna(v) else np.nan
                bs_s = get_plddt(tid, 'bindingsite')
                if bs_s:
                    for k,v in bs_s.items():
                        df.at[idx, f'{k}_bs'] = round(v,3) if pd.notna(v) else np.nan
                processed += 1
                if processed % 10 == 0:
                    self.print_progress(f"Analyzed {processed} binding sites")
            except Exception:
                continue

        out_csv = os.path.join(bindingsite_similarity_results, 'bindingsite_homology.csv')
        df.to_csv(out_csv, index=False)
        self.print_success(f"Processed {processed} binding sites → {out_csv}", start_time)
        return df

    def analyze_reference_cavity_properties(self, input='reference_bindingsite.pdb', output_results_folder='reference_cavity_analysis',
                                            probe_out=4.0, volume_cutoff=100.0):
        start_time = datetime.now()
        self.print_progress("Starting reference cavity analysis", start_time)
        self.safe_create_output_folder(output_results_folder)
        input_path = os.path.join(input)
        results = pyKVFinder.run_workflow(
            input_path, probe_out=probe_out, volume_cutoff=volume_cutoff,
            include_depth=True, include_hydropathy=True, ignore_backbone=True
        )
        results.export_all(
            fn=os.path.join(output_results_folder, "reference_results.toml"),
            output=os.path.join(output_results_folder, "reference_cavities.pdb"),
            include_frequencies_pdf=True,
            pdf=os.path.join(output_results_folder, "reference_barplots.pdf")
        )
        self.print_success("Reference cavity analysis completed", start_time)

    def analyze_cavity_properties(self,
                                  input_folder='08_homolog_bindingsites_superposed',
                                  output_cavities_folder='11_detected_cavities',
                                  output_results_folder='12_cavity_analysis_results',
                                  probe_out=4.0, volume_cutoff=100.0,
                                  csv_input='10_bindingsite_similarity_results/bindingsite_homology.csv',
                                  reference_pdb=None, reference_cavity=None):
        start_time = datetime.now()
        self.print_progress("Starting cavity detection + ICP", start_time)
        self.safe_create_output_folder(output_cavities_folder)
        self.safe_create_output_folder(output_results_folder)

        ref_pdb = reference_pdb or f"{os.path.splitext(self.reference_pdb)[0]}_bindingsite.pdb"
        ref_cav = reference_cavity or os.path.join("reference_cavity_analysis", "reference_cavities.pdb")

        def read_points(pdb_file):
            coords = []
            with open(pdb_file) as f:
                for line in f:
                    if line.startswith(("ATOM","HETATM")):
                        x,y,z = map(float, [line[30:38],line[38:46],line[46:54]])
                        coords.append([x,y,z])
            return np.array(coords)

        def run_icp(tgt_pdb, ref_pdb):
            ref_pts = read_points(ref_pdb)
            tgt_pts = read_points(tgt_pdb)
            if len(ref_pts) < 3 or len(tgt_pts) < 3:
                return np.nan, 0
            pcd_ref = o3d.geometry.PointCloud()
            pcd_ref.points = o3d.utility.Vector3dVector(ref_pts)
            pcd_tgt = o3d.geometry.PointCloud()
            pcd_tgt.points = o3d.utility.Vector3dVector(tgt_pts)
            res = o3d.pipelines.registration.registration_icp(
                pcd_tgt, pcd_ref, max_correspondence_distance=2.0,
                init=np.eye(4),
                estimation_method=o3d.pipelines.registration.TransformationEstimationPointToPoint()
            )
            return round(float(res.inlier_rmse),3), len(res.correspondence_set)

        df = pd.read_csv(csv_input) if os.path.exists(csv_input) else pd.DataFrame()

        new_cols = [
            'cavity_volume','cavity_area','cavity_avg_depth','cavity_avg_hydropathy',
            'cavity_frequency_Alipathic_apolar','cavity_frequency_Aromatic',
            'cavity_frequency_Polar_uncharged','cavity_frequency_Negatively_charged',
            'cavity_frequency_Positively_charged','cavity_frequency_Non-standard',
            'cavity_rmsd','cavity_n_points'
        ]
        for c in new_cols:
            if c not in df.columns:
                df[c] = np.nan

        pdb_files = [f for f in os.listdir(input_folder) if f.endswith('.pdb')]
        if not pdb_files:
            self.print_error(f"No PDB files in {input_folder}")
            df.fillna('n.d.').to_csv(os.path.join(output_results_folder, 'bindingsite_cavity_homology.csv'), index=False)
            return df

        self.print_progress(f"Found {len(pdb_files)} PDB files")
        processed = detected = 0

        for i, fn in enumerate(pdb_files, 1):
            tid = Path(fn).stem.split('_')[0]
            in_p = os.path.join(input_folder, fn)
            cav_pdb = os.path.join(output_cavities_folder, f"{tid}_cavities.pdb")
            if i % 5 == 0 or i == len(pdb_files):
                self.print_progress(f"Analyzing {i}/{len(pdb_files)}")
            try:
                res = pyKVFinder.run_workflow(
                    in_p, probe_out=probe_out, volume_cutoff=volume_cutoff,
                    include_depth=True, include_hydropathy=True, ignore_backbone=True
                )
                res.export_all(
                    fn=os.path.join(output_cavities_folder, f"{tid}_results.toml"),
                    output=cav_pdb,
                    include_frequencies_pdf=True,
                    pdf=os.path.join(output_cavities_folder, f"{tid}_barplots.pdf")
                )
                vol_dict = res.volume
                if not vol_dict:
                    continue
                detected += 1
                lbl = max(vol_dict, key=vol_dict.get)
                rmsd, ncorr = run_icp(cav_pdb, ref_cav)
                mask = df['homolog_id'].isin([tid, f"{tid}_superposed_bindingsite_refined"])
                if not mask.any():
                    continue
                ridx = df.index[mask].tolist()[0]
                df.loc[ridx, 'cavity_volume'] = vol_dict[lbl]
                df.loc[ridx, 'cavity_area'] = res.area.get(lbl, np.nan)
                df.loc[ridx, 'cavity_avg_depth'] = res.avg_depth.get(lbl, np.nan)
                df.loc[ridx, 'cavity_avg_hydropathy'] = res.avg_hydropathy.get(lbl, np.nan)
                freq = res.frequencies.get(lbl, {}).get("CLASS", {})
                df.loc[ridx, 'cavity_frequency_Alipathic_apolar'] = freq.get("R1", np.nan)
                df.loc[ridx, 'cavity_frequency_Aromatic'] = freq.get("R2", np.nan)
                df.loc[ridx, 'cavity_frequency_Polar_uncharged'] = freq.get("R3", np.nan)
                df.loc[ridx, 'cavity_frequency_Negatively_charged'] = freq.get("R4", np.nan)
                df.loc[ridx, 'cavity_frequency_Positively_charged'] = freq.get("R5", np.nan)
                df.loc[ridx, 'cavity_frequency_Non-standard'] = freq.get("RX", np.nan)
                df.loc[ridx, 'cavity_rmsd'] = rmsd
                df.loc[ridx, 'cavity_n_points'] = ncorr
                processed += 1
            except Exception:
                continue

        out_csv = os.path.join(output_results_folder, 'bindingsite_cavity_homology.csv')
        df.fillna('n.d.').to_csv(out_csv, index=False)
        self.print_success(f"Cavity+ICP done → {processed}/{len(pdb_files)} processed | {detected} cavities", start_time)
        return df