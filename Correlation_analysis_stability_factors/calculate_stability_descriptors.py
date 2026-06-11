import os
import sys
import pandas as pd
import numpy as np
import mdtraj as md
from pathlib import Path

def compute_features(pdb_path):
    try:
        traj = md.load(pdb_path)
        top = traj.topology
        xyz = traj.xyz[0]
        nres = top.n_residues
        natoms = top.n_atoms

        features = {}
        features["homolog_id"] = Path(pdb_path).stem

        features["num_residues"] = nres
        features["num_atoms"] = natoms

        rg = float(md.compute_rg(traj)[0])
        features["radius_of_gyration_nm"] = rg
        features["rg_per_residue"] = rg / nres if nres > 0 else 0.0

        sasa = md.shrake_rupley(traj, mode="residue")[0]
        total_sasa = float(np.sum(sasa))
        features["total_sasa_nm2"] = total_sasa
        features["avg_sasa_per_residue_nm2"] = float(np.mean(sasa)) if nres > 0 else 0.0

        hydrophobic_res = {"ALA", "VAL", "LEU", "ILE", "PHE", "TRP", "MET", "PRO", "TYR"}
        charged_res = {"ARG", "LYS", "HIS", "ASP", "GLU"}
        polar_res = {"SER", "THR", "CYS", "ASN", "GLN"}
        positive_res = {"ARG", "LYS", "HIS"}
        negative_res = {"ASP", "GLU"}

        hydrophobic_idx = [r.index for r in top.residues if r.name in hydrophobic_res]
        charged_idx = [r.index for r in top.residues if r.name in charged_res]
        polar_idx = [r.index for r in top.residues if r.name in polar_res]
        positive_idx = [r.index for r in top.residues if r.name in positive_res]
        negative_idx = [r.index for r in top.residues if r.name in negative_res]

        hydrophobic_sasa = float(np.sum(sasa[hydrophobic_idx])) if hydrophobic_idx else 0.0
        charged_sasa = float(np.sum(sasa[charged_idx])) if charged_idx else 0.0
        polar_sasa = float(np.sum(sasa[polar_idx])) if polar_idx else 0.0
        positive_sasa = float(np.sum(sasa[positive_idx])) if positive_idx else 0.0
        negative_sasa = float(np.sum(sasa[negative_idx])) if negative_idx else 0.0

        features["hydrophobic_sasa_nm2"] = hydrophobic_sasa
        features["charged_sasa_nm2"] = charged_sasa
        features["polar_sasa_nm2"] = polar_sasa
        features["positive_surface_area_nm2"] = positive_sasa
        features["negative_surface_area_nm2"] = negative_sasa

        features["hydrophobic_surface_fraction"] = hydrophobic_sasa / total_sasa if total_sasa > 0 else 0.0
        features["charged_surface_fraction"] = charged_sasa / total_sasa if total_sasa > 0 else 0.0
        features["polar_surface_fraction"] = polar_sasa / total_sasa if total_sasa > 0 else 0.0

        kd_values = {
            "ILE": 4.5, "VAL": 4.2, "LEU": 3.8, "PHE": 2.8, "CYS": 2.5, "MET": 1.9, "ALA": 1.8,
            "GLY": -0.4, "THR": -0.7, "SER": -0.8, "TRP": -0.9, "TYR": -1.3, "PRO": -1.6,
            "HIS": -3.2, "GLU": -3.5, "GLN": -3.5, "ASP": -3.5, "ASN": -3.5, "LYS": -3.9, "ARG": -4.5
        }
        
        surface_residues_kd = [kd_values.get(r.name, 0.0) for r in top.residues if sasa[r.index] > 0.05]
        features["mean_kd_hydrophobicity_surface"] = float(np.mean(surface_residues_kd)) if surface_residues_kd else 0.0

        hbonds = md.baker_hubbard(traj, periodic=False)
        num_hbonds = len(hbonds)
        features["num_hbonds_baker"] = num_hbonds
        features["hbonds_per_residue"] = num_hbonds / nres if nres > 0 else 0

        ks_hb = md.kabsch_sander(traj)
        try:
            total_ks_energy = float(ks_hb[0].sum())
        except Exception:
            total_ks_energy = float(np.sum(ks_hb[0]))
        features["total_ks_hbond_energy"] = total_ks_energy
        features["ks_hbond_energy_per_residue"] = total_ks_energy / nres if nres > 0 else 0.0

        dssp = md.compute_dssp(traj, simplified=False)[0]
        
        alpha_helix_mask = dssp == "H"
        beta_bridge_mask = dssp == "B"
        extended_strand_mask = dssp == "E"
        helix_3_10_mask = dssp == "G"
        pi_helix_mask = dssp == "I"
        turn_mask = dssp == "T"
        bend_mask = dssp == "S"
        loop_irregular_mask = (dssp == " ") | (dssp == "")

        features["alpha_helix_frac"] = float(np.mean(alpha_helix_mask))
        features["beta_bridge_frac"] = float(np.mean(beta_bridge_mask))
        features["extended_strand_frac"] = float(np.mean(extended_strand_mask))
        features["helix_3_10_frac"] = float(np.mean(helix_3_10_mask))
        features["pi_helix_frac"] = float(np.mean(pi_helix_mask))
        features["turn_frac"] = float(np.mean(turn_mask))
        features["bend_frac"] = float(np.mean(bend_mask))
        features["loop_irregular_frac"] = float(np.mean(loop_irregular_mask))
        
        features["ss_entropy"] = float(len(set(dssp)) / len(dssp)) if len(dssp) > 0 else 0.0

        features["alpha_helix_sasa_nm2"] = float(np.sum(sasa[alpha_helix_mask]))
        features["beta_bridge_sasa_nm2"] = float(np.sum(sasa[beta_bridge_mask]))
        features["extended_strand_sasa_nm2"] = float(np.sum(sasa[extended_strand_mask]))
        features["helix_3_10_sasa_nm2"] = float(np.sum(sasa[helix_3_10_mask]))
        features["pi_helix_sasa_nm2"] = float(np.sum(sasa[pi_helix_mask]))
        features["turn_sasa_nm2"] = float(np.sum(sasa[turn_mask]))
        features["bend_sasa_nm2"] = float(np.sum(sasa[bend_mask]))
        features["loop_irregular_sasa_nm2"] = float(np.sum(sasa[loop_irregular_mask]))

        features["fold_quality"] = (
            features["alpha_helix_frac"] 
            + features["extended_strand_frac"] 
            - features["loop_irregular_frac"]
        )

        disulfides = count_disulfides(top, xyz)
        features["num_disulfide_bonds"] = disulfides
        features["disulfides_per_residue"] = disulfides / nres if nres > 0 else 0
        features["disulfides_per_100_residues"] = (disulfides / nres) * 100 if nres > 0 else 0

        salt_bridges, max_cluster_size = count_salt_bridges_and_networks(top, xyz)
        features["num_salt_bridges"] = salt_bridges
        features["salt_bridges_per_residue"] = salt_bridges / nres if nres > 0 else 0
        features["largest_ionic_cluster_size"] = max_cluster_size

        hydrophobic_contacts, core_hydrophobic_contacts = count_hydrophobic_contacts_and_core(top, xyz, sasa)
        features["num_hydrophobic_contacts"] = hydrophobic_contacts
        features["hydrophobic_contacts_per_residue"] = hydrophobic_contacts / nres if nres > 0 else 0
        features["core_hydrophobic_contacts"] = core_hydrophobic_contacts
        features["core_hydrophobic_contact_fraction"] = core_hydrophobic_contacts / hydrophobic_contacts if hydrophobic_contacts > 0 else 0.0

        pi_contacts = count_pi_stacking(top, xyz)
        features["num_pi_pi_stacking"] = pi_contacts
        features["pi_pi_stacking_per_residue"] = pi_contacts / nres if nres > 0 else 0

        features["ca_clashes"] = count_ca_clashes(top, xyz)
        features["clash_energy"] = clash_energy_proxy(top, xyz)

        features["packing_density"] = float(natoms / (rg ** 3 + 1e-8))

        all_ca = [atom.index for atom in top.atoms if atom.name == "CA"]
        if len(all_ca) > 1:
            pairs_ca = [[all_ca[i], all_ca[j]] for i in range(len(all_ca)) for j in range(i + 1, len(all_ca))]
            dists_ca = md.compute_distances(traj, pairs_ca)[0]
            total_contacts = int(np.sum(dists_ca < 0.8))
            features["contacts_per_residue"] = total_contacts / nres if nres > 0 else 0.0
        else:
            features["contacts_per_residue"] = 0.0

        features["asphericity"] = float(md.asphericity(traj)[0])
        features["acylindricity"] = float(md.acylindricity(traj)[0])

        return features

    except Exception as e:
        print(f"Error processing {pdb_path}: {e}")
        return None


def count_disulfides(topology, xyz, cutoff_nm=0.23):
    sg = [a.index for a in topology.atoms if a.name == "SG" and a.residue.name == "CYS"]
    if len(sg) < 2:
        return 0
    pairs = [[i, j] for i in sg for j in sg if i < j]
    traj = md.Trajectory(xyz[np.newaxis], topology)
    d = md.compute_distances(traj, pairs)[0]
    return int(np.sum(d < cutoff_nm))


def count_salt_bridges_and_networks(topology, xyz, cutoff_nm=0.4):
    res_atoms = {}
    
    for res in topology.residues:
        if res.name in {"ARG", "LYS", "HIS"}:
            indices = [atom.index for atom in res.atoms if atom.name in {"NZ", "NH1", "NH2", "ND1", "NE2"}]
            if indices:
                res_atoms[res.index] = ("pos", indices)
        elif res.name in {"ASP", "GLU"}:
            indices = [atom.index for atom in res.atoms if atom.name in {"OD1", "OD2", "OE1", "OE2"}]
            if indices:
                res_atoms[res.index] = ("neg", indices)

    res_indices = list(res_atoms.keys())
    if len(res_indices) < 2:
        return 0, 0

    pairs = []
    pair_to_res = []
    for i in range(len(res_indices)):
        for j in range(i + 1, len(res_indices)):
            r1, r2 = res_indices[i], res_indices[j]
            t1, idxs1 = res_atoms[r1]
            t2, idxs2 = res_atoms[r2]
            if t1 != t2:
                for p in idxs1:
                    for n in idxs2:
                        pairs.append([p, n])
                        pair_to_res.append((r1, r2))

    if not pairs:
        return 0, 0

    traj = md.Trajectory(xyz[np.newaxis], topology)
    dists = md.compute_distances(traj, pairs)[0]
    
    adj = {r: set() for r in res_indices}
    salt_bridge_count = 0
    seen_pairs = set()

    for d, (r1, r2) in zip(dists, pair_to_res):
        if d < cutoff_nm:
            adj[r1].add(r2)
            adj[r2].add(r1)
            pair_key = tuple(sorted((r1, r2)))
            if pair_key not in seen_pairs:
                salt_bridge_count += 1
                seen_pairs.add(pair_key)

    visited = set()
    max_cluster_size = 0
    for res in res_indices:
        if res not in visited and (adj[res] or res_atoms[res][0] in res_atoms):
            cluster = []
            queue = [res]
            visited.add(res)
            while queue:
                curr = queue.pop(0)
                cluster.append(curr)
                for neighbor in adj[curr]:
                    if neighbor not in visited:
                        visited.add(neighbor)
                        queue.append(neighbor)
            if len(cluster) > 1:
                max_cluster_size = max(max_cluster_size, len(cluster))

    return salt_bridge_count, max_cluster_size


def count_hydrophobic_contacts_and_core(topology, xyz, sasa, cutoff_nm=0.5):
    hydrophobic = {"ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TRP", "TYR", "PRO"}

    ca_indices = []
    res_indices = []

    for res in topology.residues:
        if res.name in hydrophobic:
            for atom in res.atoms:
                if atom.name == "CA":
                    ca_indices.append(atom.index)
                    res_indices.append(res.index)

    if len(ca_indices) < 2:
        return 0, 0

    pairs = [
        [ca_indices[i], ca_indices[j]]
        for i in range(len(ca_indices))
        for j in range(i + 1, len(ca_indices))
    ]

    pair_to_res = [
        (res_indices[i], res_indices[j])
        for i in range(len(res_indices))
        for j in range(i + 1, len(res_indices))
    ]

    traj = md.Trajectory(xyz[np.newaxis], topology)
    dists = md.compute_distances(traj, pairs)[0]

    hydrophobic_contacts = 0
    core_hydrophobic_contacts = 0

    for d, (r1, r2) in zip(dists, pair_to_res):
        if d < cutoff_nm:
            hydrophobic_contacts += 1
            if sasa[r1] < 0.05 and sasa[r2] < 0.05:
                core_hydrophobic_contacts += 1

    return hydrophobic_contacts, core_hydrophobic_contacts


def count_pi_stacking(topology, xyz, cutoff_nm=0.6):
    aromatic = {"PHE", "TYR", "TRP", "HIS"}
    ca = [a.index for a in topology.atoms if a.name == "CA" and a.residue.name in aromatic]

    if len(ca) < 2:
        return 0

    pairs = [[i, j] for i in ca for j in ca if i < j]
    traj = md.Trajectory(xyz[np.newaxis], topology)
    d = md.compute_distances(traj, pairs)[0]
    return int(np.sum(d < cutoff_nm))


def count_ca_clashes(topology, xyz, cutoff_nm=0.45):
    ca = [a.index for a in topology.atoms if a.name == "CA"]
    if len(ca) < 2:
        return 0

    pairs = [[i, j] for i in ca for j in ca if i < j]
    traj = md.Trajectory(xyz[np.newaxis], topology)
    d = md.compute_distances(traj, pairs)[0]
    return int(np.sum(d < cutoff_nm))


def clash_energy_proxy(topology, xyz, cutoff_nm=0.25):
    traj = md.Trajectory(xyz[np.newaxis], topology)
    pairs = [[i, j] for i in range(topology.n_atoms) for j in range(i + 1, topology.n_atoms)]
    d = md.compute_distances(traj, pairs)[0]
    viol = cutoff_nm - d
    viol = viol[viol > 0]
    return float(np.sum(viol ** 2))


def main():
    if len(sys.argv) < 3:
        print("Usage: python calculate_stability_descriptors.py <input_folder> <output_csv>")
        sys.exit(1)

    input_folder = sys.argv[1]
    output_csv = sys.argv[2]

    if not os.path.isdir(input_folder):
        print(f"Error: {input_folder} is not a valid directory.")
        sys.exit(1)

    pdb_files = [
        f for f in os.listdir(input_folder)
        if f.lower().endswith(".pdb")
    ]

    results = []

    for pdb_file in sorted(pdb_files):
        pdb_path = os.path.join(input_folder, pdb_file)
        print(f"Processing {pdb_file}...")

        features = compute_features(pdb_path)

        if features is not None:
            results.append(features)

    if results:
        df = pd.DataFrame(results)

        cols = ["homolog_id"] + [c for c in df.columns if c != "homolog_id"]
        df = df[cols]

        df.to_csv(output_csv, index=False)

        print(f"\nAnalysis complete. {len(results)} structures processed.")
        print(f"Results saved to: {output_csv}")
    else:
        print("No structures were successfully processed.")


if __name__ == "__main__":
    main()