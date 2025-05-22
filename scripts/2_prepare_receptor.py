# 2_prepare_receptor.py
# Author: Iori Mochizuki
# Created: 2025-05-13
# Description: Fetch and prepare receptor PDB file for docking and MD

"""
Step 2: Receptor Preparation
- Download PDB structure from RCSB
- Remove water molecules
- Optionally skip PDBFixer for Meeko-sensitive structures
- Output ready-to-use receptor files for docking
"""

import os
import sys
import subprocess
import argparse

# === CLI Argument: --skip-fix ===
parser = argparse.ArgumentParser()
parser.add_argument("pdb_id", type=str, help="PDB ID to download and prepare")
parser.add_argument("--skip-fix", action="store_true", help="Skip PDBFixer (for structures like 7KDT that fail Meeko)")
args = parser.parse_args()

# === Configuration ===
PDB_ID = args.pdb_id  # Take from command-line
receptor_raw = "receptor.pdb"
receptor_clean = "receptor_clean.pdb"
receptor_fixed = "receptor_fixed.pdb"  # Will be overwritten if --skip-fix

# === Step 1: Download PDB ===
print(f"📦 Downloading {PDB_ID}...")
subprocess.run(["wget", f"https://files.rcsb.org/download/{PDB_ID}.pdb", "-O", receptor_raw])

# === Step 2: Remove water molecules with Open Babel ===
print("🧹 Removing water molecules...")
subprocess.run(["obabel", receptor_raw, "-O", receptor_clean, "--delete", "HOH"])

# === Step 3: Optional PDBFixer cleanup ===
if args.skip_fix:
    print("⚠️ Skipping PDBFixer cleanup (as requested via --skip-fix).")
    receptor_for_meeko = receptor_clean
else:
    print("🧬 Running PDBFixer to complete missing atoms/residues...")
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile

    fixer = PDBFixer(filename=receptor_clean)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.0)

    with open(receptor_fixed, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)

    print(f"✅ PDBFixer complete → {receptor_fixed}")
    receptor_for_meeko = receptor_fixed

# === Step 4: Run Meeko ===
print("🛠️ Preparing receptor with Meeko...")
meeko_command = (
    f"python3 /usr/local/envs/deeppurpose-md-env/lib/python3.11/site-packages/meeko/cli/mk_prepare_receptor.py "
    f"-i {receptor_for_meeko} -o receptor -p -j -v "
    "--box_size 20 20 20 --box_center 0 0 0 --allow_bad_res"
)
os.system(meeko_command)

print("✅ Receptor PDBQT preparation complete.")
