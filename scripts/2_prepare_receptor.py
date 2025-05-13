### 2_prepare_receptor.py
# Author: Iori Mochizuki
# Created: 2025-05-13
# Description: Fetch and prepare receptor PDB file for docking and MD

import os
from pdbfixer import PDBFixer
from openmm.app import PDBFile

# Define target protein PDB ID or local file
PDB_ID = "7JLT"  # Replace with your target ID or use a filepath
OUTPUT_FILE = "receptor_fixed.pdb"

# Step 1: Load receptor
if os.path.exists(PDB_ID):
    fixer = PDBFixer(filename=PDB_ID)
else:
    fixer = PDBFixer(pdbid=PDB_ID)

# Step 2: Apply PDBFixer operations
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(pH=7.4)

# Step 3: Save fixed receptor
with open(OUTPUT_FILE, "w") as f:
    PDBFile.writeFile(fixer.topology, fixer.positions, f)

print(f"✅ Receptor prepared and saved as {OUTPUT_FILE}")
