# 4_align_ligand.py
# Author: Iori Mochizuki
# Created: 2025-05-13
# Description: Align ligand output from AutoDock Vina to RDKit-style coordinates using Kabsch algorithm

import numpy as np
import subprocess

# === Configuration ===
docked_pdb = "output.pdb"                # From AutoDock Vina
original_pdb = "ligand.pdb"              # RDKit-generated full ligand
stripped_pdbqt = "ligand.pdbqt"          # Atom-style matching reference
aligned_pdb = "aligned_ligand_fixed.pdb" # Output (pre-clean)
final_pdb = "fixed_ligand.pdb"           # Output (final)
final_sdf = "ligand.sdf"                 # For OpenFF

# === Extract atomic coordinates ===
def extract_coordinates(pdb_file):
    atoms = []
    coords = []
    with open(pdb_file, "r") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    atom_name = line[12:16].strip()
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    atoms.append((atom_name, x, y, z, line))
                    coords.append([x, y, z])
                except ValueError:
                    print(f"⚠️ Skipping malformed line: {line.strip()}")
    return np.array(coords), atoms

# === Load ligand representations ===
docked_coords, docked_atoms = extract_coordinates(docked_pdb)
full_coords, full_atoms = extract_coordinates(original_pdb)
stripped_coords, _ = extract_coordinates(stripped_pdbqt)

print(f"✅ Loaded {len(docked_coords)} docked atoms, {len(full_coords)} full atoms, {len(stripped_coords)} stripped atoms")

# === Kabsch alignment ===
def kabsch(P, Q):
    P_centered = P - np.mean(P, axis=0)
    Q_centered = Q - np.mean(Q, axis=0)
    H = P_centered.T @ Q_centered
    U, _, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1] *= -1
        R = Vt.T @ U.T
    return R, np.mean(P, axis=0), np.mean(Q, axis=0)

# Sort atoms by distance to centroid (for alignment stability)
d1 = np.linalg.norm(docked_coords - np.mean(docked_coords, axis=0), axis=1)
d2 = np.linalg.norm(stripped_coords - np.mean(stripped_coords, axis=0), axis=1)
sorted_docked = docked_coords[np.argsort(d1)]
sorted_stripped = stripped_coords[np.argsort(d2)]

# Align coordinates
R, centroid_stripped, centroid_docked = kabsch(sorted_stripped, sorted_docked)
aligned_coords = (full_coords - centroid_stripped) @ R.T + centroid_docked

# === Write aligned structure ===
with open(aligned_pdb, "w") as f:
    for i, (_, _, _, _, line) in enumerate(full_atoms):
        x, y, z = aligned_coords[i]
        new_line = f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}"
        f.write(new_line)
print(f"✅ Saved aligned ligand to {aligned_pdb}")

# === Fix element column (OpenMM needs this) ===
def fix_element_column(input_pdb, output_pdb):
    element_map = {"C": " C", "N": " N", "O": " O", "H": " H", "S": " S", "P": " P"}
    corrected = []
    with open(input_pdb, "r") as infile:
        for line in infile:
            if line.startswith(("ATOM", "HETATM")):
                elem_guess = element_map.get(line[12:14].strip()[0], " X")
                line = line[:76] + f"{elem_guess:>2}" + line[78:]
            corrected.append(line)
    with open(output_pdb, "w") as out:
        out.writelines(corrected)
    print(f"✅ Saved element-corrected PDB to {output_pdb}")

fix_element_column(aligned_pdb, final_pdb)

# === Convert to SDF for OpenFF ===
subprocess.run(["obabel", final_pdb, "-O", final_sdf])
print(f"✅ Ligand exported to {final_sdf} for forcefield assignment")
