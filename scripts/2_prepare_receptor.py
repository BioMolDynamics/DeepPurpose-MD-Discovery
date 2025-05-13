### 2_prepare_receptor.py
# Author: Iori Mochizuki
# Created: 2025-05-13
# Description: Fetch and prepare receptor PDB file for docking and MD

"""
Step 2: Receptor Preparation
- Download PDB structure from RCSB
- Remove water molecules
- Apply PDBFixer to patch missing atoms and residues
- Add hydrogens
- Output cleaned receptor PDB
"""

import os
import sys
from pdbfixer import PDBFixer
from openmm.app import PDBFile

# Example: download structure by PDB ID (use subprocess to call wget)
import subprocess

PDB_ID = "7KDT"  # Change to your target
receptor_raw = "receptor.pdb"
receptor_clean = "receptor_clean.pdb"
receptor_final = "receptor_fixed.pdb"

# Download PDB file
subprocess.run(["wget", f"https://files.rcsb.org/download/{PDB_ID}.pdb", "-O", receptor_raw])

# Remove water using Open Babel
subprocess.run(["obabel", receptor_raw, "-O", receptor_clean, "--delete", "HOH"])

# Apply PDBFixer
fixer = PDBFixer(filename=receptor_clean)
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(pH=7.0)

# Save fixed receptor
with open(receptor_final, "w") as f:
    PDBFile.writeFile(fixer.topology, fixer.positions, f)

print("✅ Receptor preparation complete → receptor_fixed.pdb")
