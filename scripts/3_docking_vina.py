import os
import subprocess
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

# Define paths
ligand_path = "ligand.pdbqt"
receptor_path = "receptor.pdbqt"
out_path = "vina_output.pdbqt"
log_path = "vina_log.txt"

# Calculate centroid from ligand SDF (used in Step 1)
mol = Chem.SDMolSupplier("ligand.sdf", removeHs=False)[0]
conf = mol.GetConformer()
coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
centroid = coords.mean(axis=0)
print("Centroid:", centroid)

# Run AutoDock Vina with CLI (binary assumed in repo root as vina_1.2.5_linux_x86_64)
vina_exe = "./vina_1.2.5_linux_x86_64"

vina_command = [
    vina_exe,
    "--receptor", receptor_path,
    "--ligand", ligand_path,
    "--center_x", str(centroid[0]),
    "--center_y", str(centroid[1]),
    "--center_z", str(centroid[2]),
    "--size_x", "20", "--size_y", "20", "--size_z", "20",
    "--out", out_path,
    "--log", log_path
]

print("Running Vina docking...")
subprocess.run(vina_command)

# Convert top-ranked pose from PDBQT to PDB
!obabel vina_output.pdbqt -O output.pdb
print("Docking completed. Output saved to output.pdb")
