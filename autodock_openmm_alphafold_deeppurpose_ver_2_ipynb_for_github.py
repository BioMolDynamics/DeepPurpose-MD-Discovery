

## Install dependencies
"""

!pip install -q condacolab
import condacolab
condacolab.install()

"""condacolab https://pypi.org/project/condacolab/

## ***Wait for crash...***
"""

import condacolab
condacolab.check()

#!mamba create -n openff-toolkit -c conda-forge openff-toolkit
!mamba install -c conda-forge openff-toolkit -y
#!mamba install -c conda-forge openff-toolkit openff-units "numpy<1.24" -y

!mamba install -q openmm

from openff.toolkit import GLOBAL_TOOLKIT_REGISTRY
print(GLOBAL_TOOLKIT_REGISTRY.registered_toolkit_versions)

!mamba install -c conda-forge pdbfixer -y

from pdbfixer import PDBFixer

# Load a test PDB file
fixer = PDBFixer(pdbid="7JLT")  # Example: DNA structure
print("PDBFixer loaded and test PDB fetched successfully!")

"""https://autodock-vina.readthedocs.io/en/latest/docking_basic.html#preparing-the-receptor"""

!pip install meeko

!pip install -U numpy scipy rdkit vina meeko gemmi prody

import openmm
import openmm.app as app
import openmm.unit as unit
from openmm.app import PDBFile, Modeller, ForceField
from openmm.unit import nanometers

print("OpenMM version:", openmm.version.version)

!sudo apt-get install openbabel

!conda install --yes -c conda-forge openmmforcefields

!pip install prody

"""## *Reboot the runtime here.*"""

from google.colab import drive
drive.mount('/content/drive')

"""## Prepare ligands"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolfiles

# ✅ Define the SMILES string for the ligand
smiles = "CC1=C(C(=O)O[C@H](C1)[C@@H](C)[C@H]2CC[C@@H]3[C@@]2(CC[C@H]4[C@H]3C[C@@H]5[C@]6([C@@]4(C(=O)C=C[C@@H]6O)C)O5)C)CO"  #Withaferin A

# ✅ Convert SMILES to an RDKit molecule
mol = Chem.MolFromSmiles(smiles)

if mol:
    # ✅ Add explicit hydrogens (important for 3D structure & docking)
    mol = Chem.AddHs(mol)

    # ✅ Create ETKDG Embedding Parameters with a fixed random seed
    params = AllChem.ETKDG()
    params.randomSeed = 42  # Ensures reproducibility

    # ✅ Generate 3D coordinates
    status = AllChem.EmbedMolecule(mol, params)

    if status == -1:
        print("❌ ERROR: Embedding failed. Try another method.")
    else:
        print("✅ Successfully generated initial 3D coordinates!")

        # ✅ Minimize the structure to remove strain using **UFF**
        optimization_status = AllChem.UFFOptimizeMolecule(mol, maxIters=500)

        if optimization_status == 0:
            print("✅ Energy minimization successful!")
        else:
            print("⚠️ Warning: Energy minimization did not fully converge.")

        # ✅ Save the final, optimized ligand to a PDB file
        output_pdb = "ligand.pdb"
        rdmolfiles.MolToPDBFile(mol, output_pdb)
        print(f"📂 PDB file saved: {output_pdb}")

else:
    print("❌ ERROR: Invalid SMILES string! Check your input.")

!obabel /content/ligand.pdb -O /content/ligand.mol2

!obabel /content/ligand.mol2 -O /content/ligand.pdbqt

!cat /content/ligand.pdbqt

"""## Prepare the receptor"""

!wget https://files.rcsb.org/download/7KDT.pdb -O receptor.pdb #pdb ID needed

# A1 6D9H
# A2A 3EML
# A3 8X17
# POLG 3IKM
# sigma-1 receptor 5HK1
# ERK 4QTB, 5K4I, 5UMO
# NSP7-8 complex 7JLT
# NSP7-8-12 complex 7BV1
# NSP9 6W4B
# NSP9 complex 8GWB
# Cryo-EM Structure of Full-Length Drp1 Dimer (PDB ID: 8T1H)
# GMP-PNP Bound Dynamin-1-like Protein GTPase-GED Fusion (PDB ID: 4H1V)
# KEAP1 2FLU
# glucagon receptor bound to beta-arrestin 1 8JRU
# N-Protein 8R6E
# E-Protein 7K3G
# ORF3a ; 6XDC
# ORF9b ; 7KDT (complex with TOM70)
# S1 NTD ; 7B62

!obabel receptor.pdb -O receptor_clean.pdb --delete HOH

"""### Extra cleaning"""

import mdtraj as md

# ✅ Load receptor PDB
traj = md.load_pdb("receptor_clean.pdb")

# ✅ Remove non-protein residues (keep only the receptor)
filtered_traj = traj.atom_slice(traj.top.select("protein"))
filtered_traj.save_pdb("cleaned_receptor.pdb")

print("✅ Cleaned receptor saved as 'cleaned_receptor.pdb' (Non-protein residues removed).")

"""### Extra fixing"""

from pdbfixer import PDBFixer
from openmm.app import PDBFile

# ✅ Load cleaned receptor
fixer = PDBFixer("cleaned_receptor.pdb")

# ✅ Add missing atoms & hydrogens
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(pH=7.4)  # Adjust pH if needed

# ✅ Save the fixed structure
with open("fixed_receptor.pdb", "w") as f:
    PDBFile.writeFile(fixer.topology, fixer.positions, f)

print("✅ Fixed receptor saved as 'fixed_receptor.pdb'.")

!python3 /usr/local/lib/python3.11/site-packages/meeko/cli/mk_prepare_receptor.py \
-i receptor_clean.pdb -o receptor -p -j -v --box_size 20 20 20 --box_center 0 0 0 --allow_bad_res

"""## Calculate centroid"""

# Path for the receptor file
receptor_file_path = 'receptor_clean.pdb'  # Replace with your actual file path

# Define a list of target residue names or molecule identifiers to focus on
# Example: ["ZMA", "LIG"] (Adjust this list based on your receptor file)
target_residues = []

# Automatically detect residue names in the file if target_residues is empty
if not target_residues:
    with open(receptor_file_path, 'r') as file:
        for line in file:
            if line.startswith("HETATM"):
                residue_name = line[17:20].strip()  # Extract residue name
                if residue_name not in target_residues:
                    target_residues.append(residue_name)

print(f"Detected target residues: {target_residues}")

# Initialize variables for centroid calculation
x_sum = y_sum = z_sum = 0.0
atom_count = 0

# Process the receptor file to calculate centroid
with open(receptor_file_path, 'r') as file:
    for line in file:
        if line.startswith("HETATM") and any(res in line for res in target_residues):
            # Extract the X, Y, Z coordinates
            x_sum += float(line[30:38].strip())
            y_sum += float(line[38:46].strip())
            z_sum += float(line[46:54].strip())
            atom_count += 1

# Calculate the centroid
if atom_count > 0:
    x_center = x_sum / atom_count
    y_center = y_sum / atom_count
    z_center = z_sum / atom_count
    centroid_coordinates = {
        "x": round(x_center, 3),
        "y": round(y_center, 3),
        "z": round(z_center, 3)
    }
    print(f"Centroid of target residues: X={centroid_coordinates['x']}, Y={centroid_coordinates['y']}, Z={centroid_coordinates['z']}")
else:
    print("No target atoms found in the file.")
    centroid_coordinates = {"x": 0.0, "y": 0.0, "z": 0.0}

# Store the centroid coordinates for further use
print("Centroid coordinates for docking:", centroid_coordinates)

"""### When no target residues are detected and need to calculate the centroid of entire receptor"""

# Initialize variables for centroid calculation
x_sum = y_sum = z_sum = 0.0
atom_count = 0

# Calculate centroid from all ATOM entries
with open('receptor_clean.pdb', 'r') as file:
    for line in file:
        if line.startswith("ATOM"):  # Focus on ATOM entries
            x_sum += float(line[30:38].strip())
            y_sum += float(line[38:46].strip())
            z_sum += float(line[46:54].strip())
            atom_count += 1

# Calculate and print centroid
if atom_count > 0:
    x_center = x_sum / atom_count
    y_center = y_sum / atom_count
    z_center = z_sum / atom_count
    print(f"Centroid: X={x_center:.3f}, Y={y_center:.3f}, Z={z_center:.3f}")
else:
    print("No atoms found in the receptor file.")

"""### when to calculate the ATOM centroid instead"""

# Initialize variables for centroid calculation
x_sum = y_sum = z_sum = 0.0
atom_count = 0

# Path to the PDB file
receptor_file_path = "receptor_clean.pdb"  # Replace with your actual file name

# Open and process the PDB file
with open(receptor_file_path, 'r') as file:
    for line in file:
        if line.startswith("ATOM"):  # Focus on ATOM entries (standard residues)
            # Extract X, Y, Z coordinates
            x_sum += float(line[30:38].strip())
            y_sum += float(line[38:46].strip())
            z_sum += float(line[46:54].strip())
            atom_count += 1

# Calculate and display the centroid
if atom_count > 0:
    x_center = x_sum / atom_count
    y_center = y_sum / atom_count
    z_center = z_sum / atom_count
    print(f"Centroid of the protein: X={x_center:.3f}, Y={y_center:.3f}, Z={z_center:.3f}")
else:
    print("No atoms found in the file.")

"""## Install VINA"""

!wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64

!chmod +x vina_1.2.5_linux_x86_64

!./vina_1.2.5_linux_x86_64 --version

"""## Dock

### Manually enter x y z of centroid in the code below.
"""

!./vina_1.2.5_linux_x86_64 --receptor receptor.pdbqt --ligand ligand.pdbqt --out output.pdbqt --center_x 107.192 --center_y 107.312 --center_z 107.672 --size_x 30 --size_y 30 --size_z 30 --seed 12345 --exhaustiveness 20

"""## Extract only the best model from output.pdbqt"""

# ✅ Extract the first MODEL (Best Docking Pose) from PDBQT
input_pdbqt = "output.pdbqt"   # AutoDock output file
output_pdb = "best_docked_ligand.pdb"  # Output PDB file

with open(input_pdbqt, "r") as infile, open(output_pdb, "w") as outfile:
    inside_model = False
    for line in infile:
        if line.startswith("MODEL 2"):  # Stop reading after the first model
            break
        if line.startswith("MODEL 1"):  # Start recording only from MODEL 1
            inside_model = True
        if inside_model:
            if not line.startswith("MODEL") and not line.startswith("MODEL 2"):
                outfile.write(line)  # Save only relevant lines

print(f"✅ Best docking pose saved as {output_pdb}!")

"""## Clean output file and to in PDB."""

# Clean PDB file by removing docking-specific columns
with open("best_docked_ligand.pdb", "r") as infile, open("output.pdb", "w") as outfile:
    for line in infile:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            outfile.write(line)

print("✅ Cleaned output PDB saved as output.pdb")

"""# openMM"""

from pdbfixer import PDBFixer
from openmm.app import PDBFile

# Load the PDB file
fixer = PDBFixer(filename="receptor_clean.pdb")

# Find missing residues
fixer.findMissingResidues()

# If you do not want to add missing residues, clear the field
fixer.missingResidues = {}

# Continue with fixing nonstandard residues and adding atoms
fixer.findNonstandardResidues()
fixer.replaceNonstandardResidues()  # Replace nonstandard residues with standard ones
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(pH=7.0)  # Assume a neutral pH for hydrogens

# Save the fixed PDB file
with open("receptor_fixed.pdb", "w") as f:
    PDBFile.writeFile(fixer.topology, fixer.positions, f)

print("Protein preparation complete! Saved as receptor_fixed.pdb.")

from pdbfixer import PDBFixer
from openmm.app import PDBFile

# Load the PDB file
fixer = PDBFixer(filename="receptor_fixed.pdb")

# Remove all heterogens (non-standard residues) except water
fixer.removeHeterogens(keepWater=True)

# Save the cleaned PDB file
with open("receptor_cleaned.pdb", "w") as f:
    PDBFile.writeFile(fixer.topology, fixer.positions, f)

print("Non-standard residues removed. Saved as receptor_cleaned.pdb.")

"""## Minimize receptor"""

from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

# ✅ Load the receptor PDB
receptor_pdb = PDBFile("receptor_cleaned.pdb")

# ✅ Use the modern Amber14 force field
forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")

# ✅ Create the system (Amber topology)
system = forcefield.createSystem(
    receptor_pdb.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=1.0 * nanometer,
    constraints="HBonds"
)

# ✅ Set up the integrator (Langevin for stability)
integrator = LangevinMiddleIntegrator(
    300 * kelvin,   # Temperature
    1.0 / picosecond,  # Friction coefficient
    0.004 * picoseconds  # Step size
)

# ✅ Create a simulation context
simulation = Simulation(receptor_pdb.topology, system, integrator)
simulation.context.setPositions(receptor_pdb.positions)

# ✅ Perform energy minimization
print("Minimizing energy...")
simulation.minimizeEnergy(maxIterations=1000)  # You can increase iterations if needed

# ✅ Save minimized receptor structure
with open("minimized_receptor.pdb", "w") as f:
    PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), f)

print("Minimization complete! Saved as 'minimized_receptor.pdb'")

"""## Align the ligand position with docking result

### Strip down the original ligand file to pdbqt to extract the centroid
"""

!obabel ligand.pdb -O ligand.pdbqt

"""### Align"""

import numpy as np
from scipy.spatial.distance import cdist

# Input files
docked_pdb = "output.pdb"         # Docked ligand (AutoDock result)
original_pdb = "ligand.pdb"       # Full ligand with hydrogens
stripped_pdbqt = "ligand.pdbqt"   # PDBQT-like version for proper centroid calculation
aligned_pdb = "aligned_ligand_fixed.pdb"  # Output file

# Function to extract atomic coordinates & element types from PDB/PDBQT
def extract_coordinates(pdb_file):
    atoms = []
    coords = []
    with open(pdb_file, "r") as infile:
        for line in infile:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    atom_name = line[12:16].strip()
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    atoms.append((atom_name, x, y, z, line))
                    coords.append([x, y, z])
                except ValueError:
                    print(f"⚠️ Skipping malformed line: {line.strip()}")
                    continue  # Skip malformed lines
    return np.array(coords, dtype=float), atoms

# Extract ligand coordinates from 3 input files
docked_coords, docked_atoms = extract_coordinates(docked_pdb)
full_original_coords, full_original_atoms = extract_coordinates(original_pdb)
stripped_coords, stripped_atoms = extract_coordinates(stripped_pdbqt)

print(f"✅ Docked ligand has {len(docked_coords)} atoms")
print(f"✅ Full original ligand has {len(full_original_coords)} atoms")
print(f"✅ Stripped original ligand (PDBQT) has {len(stripped_coords)} atoms")

# **Step 1: Compute Centroids from Stripped PDBQT**
docked_centroid = np.mean(docked_coords, axis=0)
stripped_centroid = np.mean(stripped_coords, axis=0)  # Use PDBQT version for consistency

# **Step 2: Match Atoms by Distance from the Centroid**
docked_distances = np.linalg.norm(docked_coords - docked_centroid, axis=1)
stripped_distances = np.linalg.norm(stripped_coords - stripped_centroid, axis=1)

docked_sorted_indices = np.argsort(docked_distances)
stripped_sorted_indices = np.argsort(stripped_distances)

docked_coords_sorted = docked_coords[docked_sorted_indices]
stripped_coords_sorted = stripped_coords[stripped_sorted_indices]

# **Step 3: Kabsch Algorithm for Rotation & Translation**
def kabsch(P, Q):
    """Kabsch algorithm: finds optimal rotation to align P onto Q"""
    centroid_P = np.mean(P, axis=0)
    centroid_Q = np.mean(Q, axis=0)
    P_centered = P - centroid_P
    Q_centered = Q - centroid_Q

    H = P_centered.T @ Q_centered
    U, _, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T

    # Prevent reflections (flipping the molecule)
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T

    return R, centroid_P, centroid_Q

# Compute transformation using all heavy atoms from stripped PDBQT
R, stripped_centroid, docked_centroid = kabsch(stripped_coords_sorted, docked_coords_sorted)

# Apply transformation **without scaling**
transformed_coords = (full_original_coords - stripped_centroid) @ R.T + docked_centroid

# **Step 4: Map Back to All Atoms in the Full Structure**
aligned_atom_lines = []
for i, (name, _, _, _, original_line) in enumerate(full_original_atoms):
    x, y, z = transformed_coords[i]
    new_line = f"{original_line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{original_line[54:]}"
    aligned_atom_lines.append(new_line)

# **Step 5: Save the Corrected Ligand PDB**
with open(aligned_pdb, "w") as outfile:
    outfile.writelines(aligned_atom_lines)

print(f"✅ Aligned ligand saved as {aligned_pdb} with {len(aligned_atom_lines)} atoms!")

"""https://github.com/openmm/openmmforcefields"""

def fix_pdb_element_column(input_pdb, output_pdb):
    corrected_lines = []

    # Define element dictionary for common atoms
    element_map = {
        "C": " C", "N": " N", "O": " O", "H": " H", "S": " S", "P": " P"
    }

    with open(input_pdb, "r") as infile:
        for line in infile:
            if line.startswith(("ATOM", "HETATM")):
                atom_name = line[76:78].strip()  # Extract element (if present)
                if not atom_name:  # If empty, infer from atom type
                    atom_type = line[12:14].strip()[0]  # Extract first character
                    atom_name = element_map.get(atom_type, " X")  # Default to 'X' if unknown
                corrected_line = line[:76] + f"{atom_name:>2}" + line[78:]  # Fix element column
                corrected_lines.append(corrected_line)
            else:
                corrected_lines.append(line)

    with open(output_pdb, "w") as outfile:
        outfile.writelines(corrected_lines)

    print(f"✅ Fixed PDB saved as {output_pdb}")

# Run the fix
fix_pdb_element_column("aligned_ligand_fixed.pdb", "fixed_ligand.pdb")

!obabel fixed_ligand.pdb -O ligand.sdf

"""## MD 2ns (65mins)"""

from openmm.app import *
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from openff.toolkit import Molecule, Topology as offTopology
from openff.units.openmm import to_openmm as offquantity_to_openmm
import openmm.unit as unit
import openmm.app as app
import openmm as mm
import numpy as np
from sys import stdout
from openmm import MonteCarloBarostat
import mdtraj as md
import matplotlib.pyplot as plt

# ✅ Load minimized Sigma-1 receptor
# receptor_pdb = PDBFile("receptor_cleaned.pdb")
receptor_pdb = PDBFile("receptor_cleaned.pdb") # for all proteins other than 7BV1 or 8T1H or 8JRU
# <-- Tried without pre-minimization

# ✅ Load ligand from SDF
ligand_molecule = Molecule.from_file("ligand.sdf")

# ✅ Use OpenFF SMIRNOFF for ligand force field
smirnoff = SMIRNOFFTemplateGenerator(molecules=[ligand_molecule])
print(f"Using OpenFF Force Field: {smirnoff.smirnoff_filename}")

# ✅ Load Amber14 force field for receptor
amber_forcefield = app.ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
amber_forcefield.registerTemplateGenerator(smirnoff.generator)

# ✅ Convert ligand to OpenMM Topology
ligand_off_topology = offTopology.from_molecules(molecules=[ligand_molecule])
ligand_omm_topology = ligand_off_topology.to_openmm()
ligand_positions = offquantity_to_openmm(ligand_molecule.conformers[0])

# ✅ Merge receptor and ligand
modeller = app.Modeller(receptor_pdb.topology, receptor_pdb.positions)
modeller.add(ligand_omm_topology, ligand_positions)
modeller.addHydrogens(amber_forcefield)

# ✅ Save combined system with hydrogens
with open("combined_receptor_ligand.pdb", "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)

print("✅ Successfully merged receptor and ligand with OpenFF force field + Hydrogens!")

# ✅ Solvate with a smaller box (default 7.0 nm)
modeller.addSolvent(
    amber_forcefield, model="tip3p",
    boxSize=(7.0, 7.0, 7.0) * unit.nanometer,
    ionicStrength=0.15 * unit.molar, neutralize=True
)

# ✅ Save solvated system
with open("solvated_receptor_ligand.pdb", "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)

print("🔥 Running Quick MD Simulation (2 ns) 🔥")

# ✅ Create System with PME and constraints
system = amber_forcefield.createSystem(
    modeller.topology, nonbondedMethod=app.PME,
    nonbondedCutoff=1 * unit.nanometer, constraints=HBonds
)

# ✅ Add Barostat for Pressure Control
system.addForce(MonteCarloBarostat(1 * unit.bar, 300 * unit.kelvin, 25))

# ✅ Add Langevin Integrator (NVT)
integrator = mm.LangevinMiddleIntegrator(
    300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds
)

# ✅ Set up Simulation
platform = mm.Platform.getPlatformByName("CUDA")
simulation = app.Simulation(modeller.topology, system, integrator, platform)
simulation.context.setPositions(modeller.positions)

# ------------------------------
# 🔹 **Step 1: Energy Minimization**
# ------------------------------
print("🔹 Running Energy Minimization...")
simulation.minimizeEnergy()
simulation.reporters.append(PDBReporter("minimized.pdb", 100))

# ------------------------------
# 🔹 **Step 2: NVT Equilibration**
# ------------------------------
print("🔹 Running NVT Equilibration (1 ps)...")
simulation.context.setVelocitiesToTemperature(300 * unit.kelvin)
simulation.reporters.append(PDBReporter("nvt_equilibrated.pdb", 100))
simulation.reporters.append(StateDataReporter(stdout, 100, step=True, potentialEnergy=True, temperature=True))
simulation.step(500)  # 1 ps

# ------------------------------
# 🔹 **Step 3: NPT Equilibration**
# ------------------------------
print("🔹 Running NPT Equilibration (5 ps)...")
simulation.reporters.append(PDBReporter("npt_equilibrated.pdb", 500))
simulation.step(2500)  # 5 ps

# ------------------------------
# 🔹 **Step 4: Production MD (2 ns)**
# ------------------------------
print("🔥 Running 2 ns Production MD...")
simulation.reporters.append(DCDReporter("production_md.dcd", 1000))
simulation.reporters.append(StateDataReporter("production_md.log", 1000, step=True, potentialEnergy=True, temperature=True))
simulation.step(1000000)  # 2 ns originally 1000000. Don't forget to put it back to the original.

# ✅ Save final snapshot
with open("final_structure.pdb", "w") as f:
    PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), f)

print("🎉 Quick MD Simulation Complete!")

# -----------------------------------------------
# 📊 **RMSD & RMSF Calculation**
# -----------------------------------------------

print("📊 Calculating RMSD...")

traj = md.load("production_md.dcd", top="solvated_receptor_ligand.pdb")
ref = md.load("npt_equilibrated.pdb")

rmsd = md.rmsd(traj, ref, frame=0)

plt.figure()
plt.plot(rmsd, label="RMSD (nm)")
plt.xlabel("Frame")
plt.ylabel("RMSD (nm)")
plt.legend()
plt.savefig("quick_rmsd_plot.png")
print("✅ Quick RMSD Plot Saved!")

print("📊 Calculating RMSF...")

topology = traj.topology
ca_indices = topology.select("name CA")

rmsf = md.rmsf(traj, ref, atom_indices=ca_indices)

plt.figure()
plt.plot(rmsf, label="RMSF (nm)")
plt.xlabel("Residue Index")
plt.ylabel("RMSF (nm)")
plt.legend()
plt.savefig("quick_rmsf_plot.png")
print("✅ Quick RMSF Plot Saved!")

"""### MD just to protein"""

from openmm.app import *
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from openff.toolkit import Molecule, Topology as offTopology
from openff.units.openmm import to_openmm as offquantity_to_openmm
import openmm.unit as unit
import openmm.app as app
import openmm as mm
import numpy as np
from sys import stdout
from openmm import MonteCarloBarostat
import mdtraj as md
import matplotlib.pyplot as plt

# ✅ Load minimized Sigma-1 receptor
# receptor_pdb = PDBFile("receptor_cleaned.pdb")
receptor_pdb = PDBFile("receptor_cleaned.pdb") # for all proteins other than 7BV1 or 8T1H or 8JRU

# ✅ Load ligand from SDF
#ligand_molecule = Molecule.from_file("ligand.sdf")

# ✅ Use OpenFF SMIRNOFF for ligand force field
#smirnoff = SMIRNOFFTemplateGenerator(molecules=[ligand_molecule])
#print(f"Using OpenFF Force Field: {smirnoff.smirnoff_filename}")

# ✅ Load Amber14 force field for receptor
amber_forcefield = app.ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
#amber_forcefield.registerTemplateGenerator(smirnoff.generator)

# ✅ Convert ligand to OpenMM Topology
#ligand_off_topology = offTopology.from_molecules(molecules=[ligand_molecule])
#ligand_omm_topology = ligand_off_topology.to_openmm()
#ligand_positions = offquantity_to_openmm(ligand_molecule.conformers[0])

# ✅ Merge receptor and ligand
modeller = app.Modeller(receptor_pdb.topology, receptor_pdb.positions)
#modeller.add(ligand_omm_topology, ligand_positions)
modeller.addHydrogens(amber_forcefield)

# ✅ Save combined system with hydrogens
with open("combined_receptor_ligand.pdb", "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)

print("✅ Successfully merged receptor and ligand with OpenFF force field + Hydrogens!")

# ✅ Solvate with a smaller box (7.0 nm)
modeller.addSolvent(
    amber_forcefield, model="tip3p",
    boxSize=(7.0, 7.0, 7.0) * unit.nanometer, # 8T1H takes 9.0.
    ionicStrength=0.15 * unit.molar, neutralize=True
)

# ✅ Save solvated system
with open("solvated_receptor_ligand.pdb", "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)

print("🔥 Running Quick MD Simulation (2 ns) 🔥")

# ✅ Create System with PME and constraints
system = amber_forcefield.createSystem(
    modeller.topology, nonbondedMethod=app.PME,
    nonbondedCutoff=1 * unit.nanometer, constraints=HBonds
)

# ✅ Add Barostat for Pressure Control
system.addForce(MonteCarloBarostat(1 * unit.bar, 300 * unit.kelvin, 25))

# ✅ Add Langevin Integrator (NVT)
integrator = mm.LangevinMiddleIntegrator(
    300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds
)

# ✅ Set up Simulation
platform = mm.Platform.getPlatformByName("CUDA")
simulation = app.Simulation(modeller.topology, system, integrator, platform)
simulation.context.setPositions(modeller.positions)

# ------------------------------
# 🔹 **Step 1: Energy Minimization**
# ------------------------------
print("🔹 Running Energy Minimization...")
simulation.minimizeEnergy()
simulation.reporters.append(PDBReporter("minimized.pdb", 100))

# ------------------------------
# 🔹 **Step 2: NVT Equilibration**
# ------------------------------
print("🔹 Running NVT Equilibration (1 ps)...")
simulation.context.setVelocitiesToTemperature(300 * unit.kelvin)
simulation.reporters.append(PDBReporter("nvt_equilibrated.pdb", 100))
simulation.reporters.append(StateDataReporter(stdout, 100, step=True, potentialEnergy=True, temperature=True))
simulation.step(500)  # 1 ps

# ------------------------------
# 🔹 **Step 3: NPT Equilibration**
# ------------------------------
print("🔹 Running NPT Equilibration (5 ps)...")
simulation.reporters.append(PDBReporter("npt_equilibrated.pdb", 500))
simulation.step(2500)  # 5 ps

# ------------------------------
# 🔹 **Step 4: Production MD (? ns)**
# ------------------------------
print("🔥 Running 2 ns Production MD...")
simulation.reporters.append(DCDReporter("production_md.dcd", 1000))
simulation.reporters.append(StateDataReporter("production_md.log", 1000, step=True, potentialEnergy=True, temperature=True))
simulation.step(1000000)  # ? ns
#simulation.step(50000)  # 2 ns originally 1000000. Don't forget to put it back to the original.


# ✅ Save final snapshot
with open("final_structure.pdb", "w") as f:
    PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), f)

print("🎉 Quick MD Simulation Complete!")

"""## 🔬Check simulation by final_structure.pdb + dcd files on chimera x

## This will check H-bonds over the entire trajectory (DCD file).
"""

import mdtraj as md

# ✅ Load the trajectory and topology
traj = md.load_dcd("production_md.dcd", top="final_structure.pdb")

# ✅ Define Hydrogen Bond Criteria
h_bonds = md.baker_hubbard(traj, freq=0.1, distance_cutoff=0.35, angle_cutoff=120)

# ✅ Convert Atom Indices to Residue Names
topology = traj.topology
h_bond_list = []
for hbond in h_bonds:
    donor_idx, hydrogen_idx, acceptor_idx = hbond
    donor_residue = topology.atom(donor_idx).residue
    acceptor_residue = topology.atom(acceptor_idx).residue
    h_bond_list.append((donor_residue, acceptor_residue))

# ✅ Print H-Bond Results
print("🔹 Hydrogen Bonds Found (Donor → Acceptor):")
for donor, acceptor in h_bond_list:
    print(f"{donor} → {acceptor}")

import mdtraj as md
import numpy as np

# ✅ Load Trajectory and Topology
traj = md.load_dcd("production_md.dcd", top="final_structure.pdb")

# ✅ Identify DHEA Residue
ligand_residue_name = "UNK"  # <- Check in your PDB file!
ligand_residues = {a.residue.index for a in traj.topology.atoms if a.residue.name == ligand_residue_name}

if not ligand_residues:
    raise ValueError("❌ Ligand residue not found! Check residue name in your PDB.")

# ✅ Identify Protein Residues
protein_residues = {a.residue.index for a in traj.topology.atoms if a.residue.is_protein}

# ✅ Define Contact Criteria
cutoff_distance = 0.4  # 4.0 Å (standard for hydrophobic interactions)

# ✅ Generate Residue Pairs (DHEA vs Protein)
residue_pairs = np.array([(lig, prot) for lig in ligand_residues for prot in protein_residues])

# ✅ Compute Pairwise Contacts
contacts, _ = md.compute_contacts(traj, contacts=residue_pairs)

# ✅ Identify Most Frequent Residue Contacts
residue_contact_counts = np.sum(contacts < cutoff_distance, axis=0)

# ✅ Extract Residues in Contact
contact_residues = {}
for i, count in enumerate(residue_contact_counts):
    if count > 0:  # Contact occurred at least once
        residue = traj.topology.residue(residue_pairs[i, 1])  # Get the protein residue
        if residue not in contact_residues:
            contact_residues[residue] = 0
        contact_residues[residue] += count

# ✅ Sort by Frequency
contact_residues = sorted(contact_residues.items(), key=lambda x: x[1], reverse=True)

# ✅ Print Top 10 Interacting Residues
print("🔹 **Most Frequent Residue Contacts with Withaferin A:**")
for res, freq in contact_residues[:10]:  # Show top 10 interacting residues
    print(f"{res}: {freq} frames")

"""## 1️⃣ Residue-Ligand Interaction Frequency (%)
How often DHEA interacts with specific receptor residues over time.

👉 What it tells us:

Identifies the key binding residues stabilizing DHEA in the pocket.
Helps explain whether a specific residue is essential for interaction.
"""

import mdtraj as md
import numpy as np
from collections import defaultdict

# ✅ Load trajectory
traj = md.load_dcd("production_md.dcd", top="final_structure.pdb")

# ✅ Select ligand and protein residues
ligand_resname = "UNK"
protein_atoms = traj.topology.select("protein")
ligand_atoms = traj.topology.select(f"resname {ligand_resname}")

# ✅ Compute pairwise distances between ligand & protein
atom_pairs = np.array([(lig, prot) for lig in ligand_atoms for prot in protein_atoms])
distances = md.compute_distances(traj, atom_pairs)


# ✅ Identify residue contacts (Threshold: 0.4 nm)
contact_threshold = 0.4  # 4 Å cutoff for binding interaction
contact_frames = distances < contact_threshold  # Shape: (n_frames, n_atom_pairs)

# ✅ Count **total contacts per residue** over all frames
residue_contact_counts = defaultdict(int)

for frame_idx in range(traj.n_frames):
    for pair_idx in np.where(contact_frames[frame_idx])[0]:  # Atom pairs in contact
        protein_residue = traj.topology.atom(atom_pairs[pair_idx][1]).residue
        residue_contact_counts[protein_residue.name] += 1  # Count each residue contact

# ✅ Normalize by total contacts across all residues
total_contacts = sum(residue_contact_counts.values())
interaction_frequencies = {res: (count / total_contacts) * 100 for res, count in residue_contact_counts.items()}

# ✅ Print the correctly normalized frequencies
print("🔹 **Residue-Ligand Interaction Frequencies (%)**")
for res, freq in sorted(interaction_frequencies.items(), key=lambda x: x[1], reverse=True):
    print(f"{res}: {freq:.2f}%")

"""## 2️⃣ Ligand RMSD Over Time (Structural Stability)
How stable DHEA stays inside the binding pocket.

👉 What it tells us:

If DHEA stays in the binding site or diffuses away over time.
If the binding mode fluctuates a lot, indicating weak binding.
"""

import matplotlib.pyplot as plt
import mdtraj as md

# ✅ Load trajectory
traj = md.load_dcd("production_md.dcd", top="final_structure.pdb")

# ✅ Compute RMSD of ligand relative to frame 0
ligand_rmsd = md.rmsd(traj, traj, frame=0, atom_indices=ligand_atoms)

# ✅ Plot RMSD over time
plt.figure(figsize=(8, 5))
plt.plot(traj.time / 1000, ligand_rmsd, label="Ligand RMSD (nm)")
plt.xlabel("Time (ns)")
plt.ylabel("RMSD (nm)")
plt.title("Ligand RMSD Over Time")
plt.legend()
plt.grid()
plt.show()

"""## 3️⃣ Ligand Residence Time (How Long It Stays Bound)
How long DHEA remains within 4 Å of key residues to get binding stability insights.

👉 What it tells us:

If the ligand stays bound throughout the simulation (strong binding).
If it frequently dissociates (weak or transient binding).
"""

import mdtraj as md
import numpy as np

# ✅ Load trajectory
traj = md.load_dcd("production_md.dcd", top="final_structure.pdb")

# ✅ Select ligand and protein atoms
ligand_resname = "UNK"
ligand_atoms = traj.topology.select(f"resname {ligand_resname}")

# ✅ Correct LYS125 Selection (MDTraj indexes from 0)
correct_trp_resid = 124  # LYS125 in PDB is LYS124 in MDTraj. Check final_structure.pdb
trp_atoms = traj.topology.select(f"resname LYS and resid {correct_trp_resid}")

if len(trp_atoms) == 0:
    raise ValueError(f"❌ No atoms found for LYS125 (resid {correct_trp_resid}). Check topology!")

# ✅ Compute pairwise distances between ligand and LYS125
atom_pairs = np.array([(lig, trp) for lig in ligand_atoms for trp in trp_atoms])

if len(atom_pairs) == 0:
    raise ValueError("❌ No valid atom pairs found between ligand and LYS125.")

ligand_distances = md.compute_distances(traj, atom_pairs)

# ✅ Identify time points where ligand is within 0.4 nm (4 Å)
binding_threshold = 0.4  # 4 Å cutoff for binding interaction
bound_frames = np.any(ligand_distances < binding_threshold, axis=1)

# ✅ Compute % of time ligand is bound
binding_percentage = (np.sum(bound_frames) / traj.n_frames) * 100
print(f"🔹 Ligand Residence Time near LYS125: {binding_percentage:.2f}% of the simulation")

"""## 4️⃣ Water Network Analysis (Hydration Stability)
If your ligand binds through water bridges, we can check how often water molecules stabilize it.

👉 What it tells us:

If DHEA requires water molecules for binding.
If it gets displaced by solvent over time.
"""

# ✅ Compute water bridges (Hydration shell around ligand)
water_atoms = traj.topology.select("water")
ligand_water_distances = md.compute_distances(traj, np.array([(lig, wat) for lig in ligand_atoms for wat in water_atoms]))

# ✅ Count number of water molecules within 0.35 nm (3.5 Å)
hydration_counts = np.sum(ligand_water_distances < 0.35, axis=1)

# ✅ Plot hydration over time
plt.figure(figsize=(8, 5))
plt.plot(traj.time / 1000, hydration_counts, label="Hydration Count (Water Bridges)")
plt.xlabel("Time (ns)")
plt.ylabel("Number of Water Molecules")
plt.title("Hydration Shell Around Ligand Over Time")
plt.legend()
plt.grid()
plt.show()

"""# PCA

## PCA with ligand.
"""

import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

# === Load the MD trajectory ===
traj = md.load_dcd('production_md.dcd', top='final_structure.pdb')

# === Align trajectory to the first frame to remove overall rotation/translation ===
traj.superpose(traj, 0)

# === Extract C-alpha atom indices (focus on backbone motion) ===
ca_atoms = traj.topology.select("name CA")

# === Flatten coordinates: shape (n_frames, n_atoms * 3) ===
coords = traj.xyz[:, ca_atoms, :].reshape(traj.n_frames, -1)

# === Apply PCA ===
pca = PCA(n_components=2)
pca_result = pca.fit_transform(coords)

# === Plot first 2 principal components ===
plt.figure(figsize=(8, 5))
plt.plot(pca_result[:, 0], pca_result[:, 1], alpha=0.7)
plt.title('Principal Component Analysis of Protein Motion')
plt.xlabel('PC1 (%.2f%% variance)' % (pca.explained_variance_ratio_[0]*100))
plt.ylabel('PC2 (%.2f%% variance)' % (pca.explained_variance_ratio_[1]*100))
plt.grid(True)
plt.tight_layout()
plt.savefig("Principal Component Analysis of Protein Motion.png")
plt.show()

"""## PCA with ligand vs control"""

import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

# Load both trajectories
traj_control = md.load("production_md_no_ligand.dcd", top="final_structure_no_ligand.pdb")
traj_dhea = md.load("production_md.dcd", top="final_structure.pdb")

# Align atom selection: Select only protein backbone (avoids water/ions issue)
atom_selection = traj_control.top.select("protein and backbone")  # Cα, N, C, O atoms

traj_control = traj_control.atom_slice(atom_selection)
traj_dhea = traj_dhea.atom_slice(atom_selection)

# Ensure the same shape
xyz_control = traj_control.xyz.reshape(traj_control.n_frames, -1)
xyz_dhea = traj_dhea.xyz.reshape(traj_dhea.n_frames, -1)

# Perform PCA using control as the reference
pca = PCA(n_components=2)
pca.fit(xyz_control)

# Project both datasets onto the same PCA space
pc_control = pca.transform(xyz_control)
pc_dhea = pca.transform(xyz_dhea)

# Plot PCA projection
plt.figure(figsize=(8, 5))
plt.scatter(pc_control[:, 0], pc_control[:, 1], alpha=0.5, label="Control", color="blue")
plt.scatter(pc_dhea[:, 0], pc_dhea[:, 1], alpha=0.5, label="Withaferin A", color="red")
plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.2f}% variance)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.2f}% variance)")
plt.legend()
plt.title("PCA Projection: Withaferin A vs Control")




plt.savefig("Principal Component Analysis with ligand vs control.png")
plt.show()

"""## FEL Calculation & Visualization"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import Boltzmann

# Define constants
T = 310.15  # Temperature in Kelvin (human body ~37°C)
kB = Boltzmann  # Boltzmann constant in J/K

def compute_fel(pc1, pc2, bins=100):
    """
    Compute the Free Energy Landscape (FEL) from PCA components.
    """
    H, xedges, yedges = np.histogram2d(pc1, pc2, bins=bins, density=True)  # 2D histogram
    P = H.T  # Transpose for correct visualization
    P[P == 0] = 1e-12  # Avoid log(0) errors
    G = -kB * T * np.log(P)  # Compute free energy

    return G, xedges, yedges

def plot_fel(G, xedges, yedges, filename, title="Free Energy Landscape"):
    """
    Plot FEL as a contour map and save as PNG.
    """
    X, Y = np.meshgrid(xedges[:-1], yedges[:-1])  # Create grid from bin edges

    plt.figure(figsize=(8, 6))
    plt.contourf(X, Y, G, levels=30, cmap="viridis")
    plt.colorbar(label="Free Energy (J)")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title(title)

    # Save as PNG
    plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.show()

# Compute Free Energy Landscapes
G_control, xedges_c, yedges_c = compute_fel(pc_control[:, 0], pc_control[:, 1])
G_dhea, xedges_d, yedges_d = compute_fel(pc_dhea[:, 0], pc_dhea[:, 1])

# Save and Download Free Energy Landscapes as PNG
plot_fel(G_control, xedges_c, yedges_c, "FEL_Control.png", title="Free Energy Landscape - Control")
plot_fel(G_dhea, xedges_d, yedges_d, "FEL_Withaferin A.png", title="Free Energy Landscape - Withaferin A")

"""# DeepPurpose

### Install dependencies
"""

!pip install git+https://github.com/bp-kelley/descriptastorus
!pip install pandas-flavor
!pip install subword-nmt

#This has to be installed separately for some reason.
!pip install rdkit-pypi

!git clone https://github.com/mosmos6/Deeppurpose
!pip install ./Deeppurpose

#must run this again for some reason
!pip install rdkit-pypi

# This errors at first, so reinstall rdkit again and come back here.
import DeepPurpose.oneliner as oneliner
from DeepPurpose import dataset
import time

"""### Data processing from the original BindingDB tsv

#### Cleaning
"""

import pandas as pd

tsv_path = '/content/BindingDB_Covid-19.tsv'

# ✅ Load raw lines
with open(tsv_path, 'r', encoding='utf-8', errors='ignore') as f:
    raw_lines = f.readlines()

# ✅ Keep the first line as header
header_line = raw_lines[0].strip()

# ✅ Filter: keep the first line + all other lines that are not repeated headers
cleaned_lines = [raw_lines[0]] + [line for line in raw_lines[1:] if line.strip() != header_line]

# ✅ Save to a new file
cleaned_path = "cleaned_bindingdb.tsv"
with open(cleaned_path, 'w', encoding='utf-8') as f:
    f.writelines(cleaned_lines)

print(f"✅ Cleaned TSV saved as {cleaned_path}.")


df_preview = pd.read_csv(cleaned_path, sep='\t', engine='python', nrows=1)
print(df_preview.columns.tolist())

"""#### Rearrange columns"""

import pandas as pd

# Only select key columns (some rows may still contain junk fields)
usecols = ['Ligand SMILES', 'BindingDB Target Chain Sequence', 'IC50 (nM)']

df = pd.read_csv(cleaned_path, sep='\t', usecols=usecols, engine='python')
df = df.dropna()

# Optional: Clean and convert IC50 to molar
df = df[df['Ligand SMILES'].apply(lambda x: isinstance(x, str) and len(x) > 0)]
df = df[df['BindingDB Target Chain Sequence'].apply(lambda x: isinstance(x, str) and len(x) > 0)]
df['IC50 (nM)'] = df['IC50 (nM)'].astype(str).str.extract(r'([\d.]+)').astype(float)
df = df.dropna()
df['Affinity'] = df['IC50 (nM)'] * 1e-9  # convert to molar
df = df[['Ligand SMILES', 'BindingDB Target Chain Sequence', 'Affinity']]
df.columns = ['SMILES', 'Target Sequence', 'Affinity']

df.to_csv("processed_bindingdb.tsv", sep='\t', index=False)


print(f"✅ Cleaned and loaded {len(df)} compound-target pairs.")

"""#### Inspect dataset"""

from DeepPurpose.dataset import load_bindingdb_covid_tsv

X = load_bindingdb_covid_tsv("/content/processed_bindingdb.tsv")

"""## ProtTrans

### Install ProtTrans
"""

from transformers import T5Tokenizer, T5EncoderModel

tokenizer = T5Tokenizer.from_pretrained("Rostlab/prot_t5_xl_uniref50", do_lower_case=False)
model = T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_uniref50")
model = model.eval()  # turn off dropout

"""### Inspect dataset"""

import pandas as pd

df = pd.read_csv('/content/processed_bindingdb.tsv', sep='\t')
unique_targets = df['Target Sequence'].unique()
print(f"🔎 Unique protein sequences: {len(unique_targets)}")

"""### Embed proteins"""

# This needs A100

from tqdm import tqdm
import numpy as np
import torch

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model = model.to(device)  # move model to GPU

protein_embeddings = {}

with torch.no_grad():
    for seq in tqdm(unique_targets):
        seq_clean = seq.replace(" ", "")
        seq_tok = " ".join(list(seq_clean))

        ids = tokenizer(seq_tok, return_tensors="pt", padding=True).to(device)  # move inputs to GPU
        embedding = model(**ids).last_hidden_state  # shape: [1, L, 1024]
        pooled = torch.mean(embedding, dim=1)       # shape: [1, 1024]

        protein_embeddings[seq_clean] = pooled.squeeze().cpu().numpy()  # move back to CPU before saving

"""### Save embedded proteins"""

import pickle

# Save to file
with open("protein_embeddings.pkl", "wb") as f:
    pickle.dump(protein_embeddings, f)

print("✅ Embeddings saved to protein_embeddings.pkl")

"""### Add embedded proteins to the dataset"""

import pandas as pd
import pickle

# Load processed BindingDB data
df = pd.read_csv("processed_bindingdb.tsv", sep='\t')

# Load protein embeddings
with open("protein_embeddings.pkl", "rb") as f:
    protein_embeddings = pickle.load(f)

# Map embedding vector to each row
df["ProtTrans"] = df["Target Sequence"].map(protein_embeddings)

# Save extended TSV for backup
df.to_pickle("embedded_bindingdb.pkl")
print("✅ Added ProtTrans embeddings to dataframe.")

"""## Extract top results from the dataset"""

import pandas as pd
import matplotlib.pyplot as plt

# Load your filtered strong binders
df = pd.read_csv("strong_binders_cleaned.csv")

# Count protein frequency
top_targets = df["Target Sequence"].value_counts()

# Print quick stats
print(f"🧬 Total unique proteins: {len(top_targets)}")
print("\n🔝 Top proteins by occurrence:")
print(top_targets.head(10))

# Save to CSV
top_targets.to_csv("top_protein_targets.csv", header=["Count"])
print("✅ Saved protein frequency table to top_protein_targets.csv")

# Optional: Visualize top 10
plt.figure(figsize=(10, 4))
top_targets.head(10).plot(kind="barh", title="Top 10 Protein Targets")
plt.xlabel("Count")
plt.gca().invert_yaxis()
plt.grid(True)
plt.tight_layout()
plt.show()

import pandas as pd

# Load the cleaned strong binders (this should have many entries per protein)
df = pd.read_csv("strong_binders_cleaned.csv")

# Group by Target Sequence and sample up to 200 per target
df_sampled = (
    df.groupby("Target Sequence", group_keys=False)
    .apply(lambda g: g.sample(n=min(len(g), 150), random_state=42))
)

# Save the result
df_sampled.to_csv("strong_binders_top150_per_protein.csv", index=False)
print(f"✅ Final dataset shape: {df_sampled.shape}")

# Verify top counts
top_counts = df_sampled["Target Sequence"].value_counts()
print("✅ Top 10 target protein sample counts:\n")
print(top_counts.head(10))

"""## Model Training"""

from DeepPurpose import utils, CompoundPred
import pandas as pd
from rdkit import Chem
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import numpy as np
import re
import pickle

# Load dataset
df = pd.read_csv("/content/strong_binders_top200_per_protein.csv")

# 🧪 Step 1: Clean SMILES
def is_valid_smiles(smi):
    try:
        mol = Chem.MolFromSmiles(smi)
        return mol is not None
    except:
        return False

print(f"🧪 Original size: {len(df)}")
df = df[df["SMILES"].apply(is_valid_smiles)].reset_index(drop=True)
print(f"✅ After SMILES cleaning: {len(df)}")

# 🧪 Step 2: Clean missing values
df = df.dropna(subset=["SMILES", "Target Sequence", "Affinity", "ProtTrans"])
df['Affinity'] = df['Affinity'] * 1e+9  # convert to nM
# After loading and scaling to nM
df['Affinity'] = 9 - np.log10(df['Affinity'])  # convert to pIC50
scaler = StandardScaler()
df['Affinity'] = scaler.fit_transform(df['Affinity'].values.reshape(-1, 1))



# 🧪 Step 3: Parse embedding column (from string to np.array)
# Load original embedding dictionary
with open("protein_embeddings.pkl", "rb") as f:
    prot_embed_dict = pickle.load(f)

# Add back correct 1024-dim vectors (not corrupted by CSV)
df["ProtTrans"] = df["Target Sequence"].map(prot_embed_dict)

# Drop invalid
df = df[df["ProtTrans"].apply(lambda x: isinstance(x, np.ndarray) and len(x) == 1024)].reset_index(drop=True)
print(f"✅ Successfully re-injected {len(df)} true 1024-dim embeddings.")

# 🧪 Step 4: Reduce ProtTrans dimension to 512 with PCA
all_embeddings = np.stack(df["ProtTrans"].values)
pca = PCA(n_components=512, random_state=42)
reduced = pca.fit_transform(all_embeddings)

# Replace original embeddings with reduced ones
df["ProtTrans"] = list(reduced)
print("✅ PCA reduction complete. New shape:", reduced.shape)

# 🔁 Step 4-2: Augment AFTER PCA so ProtTrans column exists
top_15_thresh = df["Affinity"].quantile(0.85)
high_df = df[df["Affinity"] > top_15_thresh]

df["source"] = "original"
high_df["source"] = "boosted"

df = pd.concat([df, high_df, high_df]).reset_index(drop=True)
df["ProtTrans"] = list(pca.transform(np.stack(df["Target Sequence"].map(prot_embed_dict))))
df = df.sample(frac=1.0, random_state=42).reset_index(drop=True)
df.loc[913:, "source"] = "boosted"



# 🧪 Step 5: Proceed to DeepPurpose training
train, val, test = utils.data_process(
    X_drug = df["SMILES"].values,
    X_target = df["Target Sequence"].values,
    y = df["Affinity"].values,
    drug_encoding = "Transformer",
    target_encoding = "ProtTrans",  # <- using reduced vectors
    split_method = "random",
    frac = [0.7, 0.1, 0.2]
)

config = utils.generate_config(
    drug_encoding="Transformer",
    target_encoding="ProtTrans",
    cls_hidden_dims=[512, 256, 128, 128, 128, 128, 256],  # 512, 256, 128, 128, 128, 128, 256
    transformer_emb_size_drug = 512,
    transformer_intermediate_size_drug = 2048,
    transformer_num_attention_heads_drug = 8,
    transformer_n_layer_drug = 4,
    train_epoch=50,
    LR=1e-4,
    batch_size=128
)

model = CompoundPred.model_initialize(**config)
model.train(train, val, test)

test_df = df.iloc[test.index]  # use test.index from data_process()
print("🔍 Test set breakdown:")
print(test_df["source"].value_counts())

"""## Save model"""

# 🧠 Save trained DeepPurpose model to a folder
import os

save_path = "/content/deeppurpose_model_saved"
os.makedirs(save_path, exist_ok=True)

model.save_model(save_path)
print(f"✅ Model saved to: {save_path}")

"""## Save model folder"""

import os
import shutil

# ✅ Define the source and destination
src_folder = "/content/deeppurpose_model_saved"
dst_folder = "/content/drive/My Drive/DeepPurpose/deeppurpose_model_saved"

# ✅ If the destination folder exists, remove it to avoid overwrite issues
if os.path.exists(dst_folder):
    shutil.rmtree(dst_folder)

# ✅ Copy the entire folder to Google Drive
shutil.copytree(src_folder, dst_folder)
print(f"✅ Folder uploaded successfully to: {dst_folder}")

"""## Load model"""

from DeepPurpose import CompoundPred

# ✅ Path to your saved model
model_path = "/content/drive/My Drive/DeepPurpose_model/deeppurpose_model_saved"

# ✅ Load the model
model = CompoundPred.model_pretrained(path_dir=model_path)
print("✅ Model loaded successfully!")

"""# Deeppurpose on SARS2 proteins

## STEP 1: Embed all FASTA sequences using ProtTrans

### FASTA extraction
"""

import pandas as pd

fasta_path = "/content/protein.faa"

# Manual parsing of .faa FASTA format
sequences = []
with open(fasta_path, "r") as f:
    current_seq = ""
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if current_seq:
                sequences.append(current_seq)
                current_seq = ""
        else:
            current_seq += line
    if current_seq:
        sequences.append(current_seq)

# Convert to DataFrame and save
df = pd.DataFrame(sequences)
csv_path = "/content/fasta_sequences.csv"
df.to_csv(csv_path, index=False)

csv_path

import pandas as pd

# 1. Load your existing CSV file (created from the .faa file)
existing_csv_path = "/content/fasta_sequences.csv"
df_existing = pd.read_csv(existing_csv_path, header=None)


# 2. Your manually collected FASTA list
manual_fasta = ["MESLVPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHLKDGTCGLVEVEKGVLPQLEQPYVFIKRSDARTAPHGHVMVELVAELEGIQYGRSGETLGVLVPHVGEIPVAYRKVLLRKNGNKGAGGHSYGADLKSFDLGDELGTDPYEDFQENWNTKHSSGVTRELMRELNGG",
"AYTRYVDNNFCGPDGYPLECIKDLLARAGKASCTLSEQLDFIDTKRGVYCCREHEHEIAWYTERSEKSYELQTPFEIKLAKKFDTFNGECPNFVFPLNSIIKTIQPRVEKKKLDGFMGRIRSVYPVASPNECNQMCLSTLMKCDHCGETSWQTGDFVKATCEFCGTENLTKEGATTCGYLPQNAVVKIYCPACHNSEVGPEHSLAEYHNESGLKTILRKGGRTIAFGGCVFSYVGCHNKCAYWVPRASANIGCNHTGVVGEGSEGLNDNLLEILQKEKVNINIVGDFKLNEEIAIILASFSASTSAFVETVKGLDYKAFKQIVESCGNFKVTKGKAKKGAWNIGEQKSILSPLYAFASEAARVVRSIFSRTLETAQNSVRVLQKAAITILDGISQYSLRLIDAMMFTSDLATNNLVVMAYITGGVVQLTSQWLTNIFGTVYEKLKPVLDWLEEKFKEGVEFLRDGWEIVKFISTCACEIVGGQIVTCAKEIKESVQTFFKLVNKFLALCADSIIIGGAKLKALNLGETFVTHSKGLYRKCVKSREETGLLMPLKAPKEIIFLEGETLPTEVLTEEVVLKTGDLQPLEQPTSEAVEAPLVGTPVCINGLMLLEIKDTEKYCALAPNMMVTNNTFTLKGG",
"APTKVTFGDDTVIEVQGYKSVNITFELDERIDKVLNEKCSAYTVELGTEVNEFACVVADAVIKTLQPVSELLTPLGIDLDEWSMATYYLFDESGEFKLASHMYCSFYPPDEDEEEGDCEEEEFEPSTQYEYGTEDDYQGKPLEFGATSAALQPEEEQEEDWLDDDSQQTVGQQDGSEDNQTTTIQTIVEVQPQLEMELTPVVQTIEVNSFSGYLKLTDNVYIKNADIVEEAKKVKPTVVVNAANVYLKHGGGVAGALNKATNNAMQVESDDYIATNGPLKVGGSCVLSGHNLAKHCLHVVGPNVNKGEDIQLLKSAYENFNQHEVLLAPLLSAGIFGADPIHSLRVCVDTVRTNVYLAVFDKNLYDKLVSSFLEMKSEKQVEQKIAEIPKEEVKPFITESKPSVEQRKQDDKKIKACVEEVTTTLEETKFLTENLLLYIDINGNLHPDSATLVSDIDITFLKKDAPYIVGDVVQEGVLTAVVIPTKKAGGTTEMLAKALRKVPTDNYITTYPGQGLNGYTVEEAKTVLKKCKSAFYILPSIISNEKQEILGTVSWNLREMLAHAEETRKLMPVCVETKAIVSTIQRKYKGIKIQEGVVDYGARFYFYTSKTTVASLINTLNDLNETLVTMPLGYVTHGLNLEEAARYMRSLKVPATVSVSSPDAVTAYNGYLTSSSKTPEEHFIETISLAGSYKDWSYSGQSTQLGIEFLKRGDKSVYYTSNPTTFHLDGEVITFDNLKTLLSLREVRTIKVFTTVDNINLHTQVVDMSMTYGQQFGPTYLDGADVTKIKPHNSHEGKTFYVLPNDDTLRVEAFEYYHTTDPSFLGRYMSALNHTKKWKYPQVNGLTSIKWADNNCYLATALLTLQQIELKFNPPALQDAYYRARAGEAANFCALILAYCNKTVGELGDVRETMSYLFQHANLDSCKRVLNVVCKTCGQQQTTLKGVEAVMYMGTLSYEQFKKGVQIPCTCGKQATKYLVQQESPFVMMSAPPAQYELKHGTFTCASEYTGNYQCGHYKHITSKETLYCIDGALLTKSSEYKGPITDVFYKENSYTTTIKPVTYKLDGVVCTEIDPKLDNYYKKDNSYFTEQPIDLVPNQPYPNASFDNFKFVCDNIKFADDLNQLTGYKKPASRELKVTFFPDLNGDVVAIDYKHYTPSFKKGAKLLHKPIVWHVNNATNKATYKPNTWCIRCLWSTKPVETSNSFDVLKSEDAQGMDNLACEDLKPVSEEVVENPTIQKDVLECNVKTTEVVGDIILKPANNSLKITEEVGHTDLMAAYVDNSSLTIKKPNELSRVLGLKTLATHGLAAVNSVPWDTIANYAKPFLNKVVSTTTNIVTRCLNRVCTNYMPYFFTLLLQLCTFTRSTNSRIKASMPTTIAKNTVKSVGKFCLEASFNYLKSPNFSKLINIIIWFLLLSVCLGSLIYSTAALGVLMSNLGMPSYCTGYREGYLNSTNVTIATYCTGSIPCSVCLSGLDSLDTYPSLETIQITISSFKWDLTAFGLVAEWFLAYILFTRFFYVLGLAAIMQLFFSYFAVHFISNSWLMWLIINLVQMAPISAMVRMYIFFASFYYVWKSYVHVVDGCNSSTCMMCYKRNRATRVECTTIVNGVRRSFYVYANGGKGFCKLHNWNCVNCDTFCAGSTFISDEVARDLSLQFKRPINPTDQSSYIVDSVTVKNGSIHLYFDKAGQKTYERHSLSHFVNLDNLRANNTKGSLPINVIVFDGKSKCEESSAKSASVYYSQLMCQPILLLDQALVSDVGDSAEVAVKMFDAYVNTFSSTFNVPMEKLKTLVATAEAELAKNVSLDNVLSTFISAARQGFVDSDVETKDVVECLKLSHQSDIEVTGDSCNNYMLTYNKVENMTPRDLGACIDCSARHINAQVAKSHNIALIWNVKDFMSLSEQLRKQIRSAAKKNNLPFKLTCATTRQVVNVVTTKIALKGG",
"KIVNNWLKQLIKVTLVFLFVAAIFYLITPVHVMSKHTDFSSEIIGYKAIDGGVTRDIASTDTCFANKHADFDTWFSQRGGSYTNDKACPLIAAVITREVGFVVPGLPGTILRTTNGDFLHFLPRVFSAVGNICYTPSKLIEYTDFATSACVLAAECTIFKDASGKPVPYCYDTNVLEGSVAYESLRPDTRYVLMDGSIIQFPNTYLEGSVRVVTTFDSEYCRHGTCERSEAGVCVSTSGRWVLNNDYYRSLPGVFCGVDAVNLLTNMFTPLIQPIGALDISASIVAGGIVAIVVTCLAYYFMRFRRAFGEYSHVVAFNTLLFLMSFTVLCLTPVYSFLPGVYSVIYLYLTFYLTNDVSFLAHIQWMVMFTPLVPFWITIAYIICISTKHFYWFFSNYLKRRVVFNGVSFSTFEEAALCTFLLNKEMYLKLRSDVLLPLTQYNRYLALYNKYKYFSGAMDTTSYREAACCHLAKALNDFSNSGSDVLYQPPQTSITSAVLQ",
"SGFRKMAFPSGKVEGCMVQVTCGTTTLNGLWLDDVVYCPRHVICTSEDMLNPNYEDLLIRKSNHNFLVQAGNVQLRVIGHSMQNCVLKLKVDTANPKTPKYKFVRIQPGQTFSVLACYNGSPSGVYQCAMRPNFTIKGSFLNGSCGSVGFNIDYDCVSFCYMHHMELPTGVHAGTDLEGNFYGPFVDRQTAQAAGTDTTITVNVLAWLYAAVINGDRWFLNRFTTTLNDFNLVAMKYNYEPLTQDHVDILGPLSAQTGIAVLDMCASLKELLQNGMNGRTILGSALLEDEFTPFDVVRQCSGVTFQ",
"SAVKRTIKGTHHWLLLTILTSLLVLVQSTQWSLFFFLYENAFLPFAMGIIAMSAFAMMFVKHKHAFLCLFLLPSLATVAYFNMVYMPASWVMRIMTWLDMVDTSLSGFKLKDCVMYASAVVLLILMTARTVYDDGARRVWTLMNVLTLVYKVYYGNALDQAISMWALIISVTSNYSGVVTTVMFLARGIVFMCVEYCPIFFITGNTLQCIMLVYCFLGYFCTCYFGLFCLLNRYFRLTLGVYDYLVSTQEFRYMNSQGLLPPKNSIDAFKLNIKLLGVGGKPCIKVATVQ",
"SKMSDVKCTSVVLLSVLQQLRVESSSKLWAQCVQLHNDILLAKDTTEAFEKMVSLLSVLLSMQGAVDINKLCEEMLDNRATLQ",
"AIASEFSSLPSYAAFATAQEAYEQAVANGDSEVVLKKLKKSLNVAKSEFDRDAAMQRKLEKMADQAMTQMYKQARSEDKRAKVTSAMQTMLFTMLRKLDNDALNNIINNARDGCVPLNIIPLTTAAKLMVVIPDYNTYKNTCDGTTFTYASALWEIQQVVDADSKIVQLSEISMDNSPNLAWPLIVTALRANSAVKLQ",
"NNELSPVALRQMSCAAGTTQTACTDDNALAYYNTTKGGRFVLALLSDLQDLKWARFPKSDGTGTIYTELEPPCRFVTDTPKGPKVKYLYFIKGLNNLNRGMVLGSLAATVRLQ",
"AGNATEVPANSTVLSFCAFAVDAAKAYKDYLASGGQPITNCVKMLCTHTGTGQAITVTPEANMDQESFGGASCCLYCRCHIDHPNPKGFCDLKGKYVQIPTTCANDPVGFTLKNTVCTVCGMWKGYGCSCDQLREPMLQ",
"SADAQSFLNRVCGVSAARLTPCGTGTSTDVVYRAFDIYNDKVAGFAKFLKTNCCRFQEKDEDDNLIDSYFVVKRHTFSNYQHEETIYNLLKDCPAVAKHDFFKFRIDGDMVPHISRQRLTKYTMADLVYALRHFDEGNCDTLKEILVTYNCCDDDYFNKKDWYDFVENPDILRVYANLGERVRQALLKTVQFCDAMRNAGIVGVLTLDNQDLNGNWYDFGDFIQTTPGSGVPVVDSYYSLLMPILTLTRALTAESHVDTDLTKPYIKWDLLKYDFTEERLKLFDRYFKYWDQTYHPNCVNCLDDRCILHCANFNVLFSTVFPPTSFGPLVRKIFVDGVPFVVSTGYHFRELGVVHNQDVNLHSSRLSFKELLVYAADPAMHAASGNLLLDKRTTCFSVAALTNNVAFQTVKPGNFNKDFYDFAVSKGFFKEGSSVELKHFFFAQDGNAAISDYDYYRYNLPTMCDIRQLLFVVEVVDKYFDCYDGGCINANQVIVNNLDKSAGFPFNKWGKARLYYDSMSYEDQDALFAYTKRNVIPTITQMNLKYAISAKNRARTVAGVSICSTMTNRQFHQKLLKSIAATRGATVVIGTSKFYGGWHNMLKTVYSDVENPHLMGWDYPKCDRAMPNMLRIMASLVLARKHTTCCSLSHRFYRLANECAQVLSEMVMCGGSLYVKPGGTSSGDATTAYANSVFNICQAVTANVNALLSTDGNKIADKYVRNLQHRLYECLYRNRDVDTDFVNEFYAYLRKHFSMMILSDDAVVCFNSTYASQGLVASIKNFKSVLYYQNNVFMSEAKCWTETDLTKGPHEFCSQHTMLVKQGDDYVYLPYPDPSRILGAGCFVDDIVKTDGTLMIERFVSLAIDAYPLTKHPNQEYADVFHLYLQYIRKLHDELTGHMLDMYSVMLTNDNTSRYWEPEFYEAMYTPHTVLQ",
"AVGACVLCNSQTSLRCGACIRRPFLCCKCCYDHVISTSHKLVLSVNPYVCNAPGCDVTDVTQLYLGGMSYYCKSHKPPISFPLCANGQVFGLYKNTCVGSDNVTDFNAIATCDWTNAGDYILANTCTERLKLFAAETLKATEETFKLSYGIATVREVLSDRELHLSWEVGKPRPPLNRNYVFTGYRVTKNSKVQIGEYTFEKGDYGDAVVYRGTTTYKLNVGDYFVLTSHTVMPLSAPTLVPQEHYVRITGLYPTLNISDEFSSNVANYQKVGMQKYSTLQGPPGTGKSHFAIGLALYYPSARIVYTACSHAAVDALCEKALKYLPIDKCSRIIPARARVECFDKFKVNSTLEQYVFCTVNALPETTADIVVFDEISMATNYDLSVVNARLRAKHYVYIGDPAQLPAPRTLLTKGTLEPEYFNSVCRLMKTIGPDMFLGTCRRCPAEIVDTVSALVYDNKLKAHKDKSAQCFKMFYKGVITHDVSSAINRPQIGVVREFLTRNPAWRKAVFISPYNSQNAVASKILGLPTQTVDSSQGSEYDYVIFTQTTETAHSCNVNRFNVAITRAKVGILCIMSDRDLYDKLQFTSLEIPRRNVATLQ",
"AENVTGLFKDCSKVITGLHPTQAPTHLSVDTKFKTEGLCVDIPGIPKDMTYRRLISMMGFKMNYQVNGYPNMFITREEAIRHVRAWIGFDVEGCHATREAVGTNLPLQLGFSTGVNLVAVPTGYVDTPNNTDFSRVSAKPPPGDQFKHLIPLMYKGLPWNVVRIKIVQMLSDTLKNLSDRVVFVLWAHGFELTSMKYFVKIGPERTCCLCDRRATCFSTASDTYACWHHSIGFDYVYNPFMIDVQQWGFTGNLQSNHDLYCQVHGNAHVASCDAIMTRCLAVHECFVKRVDWTIEYPIIGDELKINAACRKVQHMVVKAALLADKFPVLHDIGNPKAIKCVPQADVEWKFYDAQPCSDKAYKIEELFYSYATHSDKFTDGVCLFWNCNVDRYPANSIVCRFDTRVLSNLNLPGCDGGSLYVNKHAFHTPAFDKSAFVNLKQLPFFYYSDSPCESHGKQVVSDIDYVPLKSATCITRCNLGGAVCRHHANEYRLYLDAYNMMISAGFSLWVYKQFDTYNLWNTFTRLQ",
"SLENVAFNVVNKGHFDGQQGEVPVSIINNTVYTKVDGVDVELFENKTTLPVNVAFELWAKRNIKPVPEVKILNNLGVDIAANTVIWDYKRDAPAHISTIGVCSMTDIAKKPTETICAPLTVFFDGRVDGQVDLFRNARNGVLITEGSVKGLQPSVGPKQASLNGVTLIGEAVKTQFNYYKKVDGVVQQLPETYFTQSRNLQEFKPRSQMEIDFLELAMDEFIERYKLEGYAFEHIVYGDFSHSQLGGLHLLIGLAKRFKESPFELEDFIPMDSTVKNYFITDAQTGSSKCVCSVIDLLLDDFVEIIKSQDLSVVSKVVKVTIDYTEISFMLWCKDGHVETFYPKLQ",
"SSQAWQPGVAMPNLYKMQRMLLEKCDLQNYGDSATLPKGIMMNVAKYTQLCQYLNTLTLAVPYNMRVIHFGAGSDKGVAPGTAVLRQWLPTGTLLVDSDLNDFVSDADSTLIGDCATVHTANKWDLIISDMYDPKTKNVTKENDSKEGFFTYICGFIQQKLALGGSVAIKITEHSWNADLYKLMGHFAWWTAFVTNVNASSSEAFLIGCNYLGKPREQIDGYVMHANYIFWRNTNPIQLSSYSLFDMSKFPLKLRGTAVMSLKEGQINDMILSLLSKGRLIIRENNRVVISSDVLVNN",
"MDPKISEMHPALRLVDPQIQLAVTRMENAVGRDQNNVGPKVYPIILRLGSPLSLNMARKTLNSLEDKAFQLTPIAVQMTKLATTEELPDEFVVVTVK",
"MLQSCYNFLKEQHCQKASTQKGAEAAVKPLLVPHHVVATVQEIQLQAAVGELLLLEWLAMAVMLLLLCCCLTD"]

# 3. Convert manual list to DataFrame
df_manual = pd.DataFrame(manual_fasta)

# 4. Combine both
df_combined = pd.concat([df_existing, df_manual], ignore_index=True)
df_combined.columns = ["Target Sequence"]
df_combined = df_combined.drop(index=0)

# 5. Save to new CSV
output_path = "/content/combined_fasta_targets.csv"
df_combined.to_csv(output_path, index=False)

output_path

df=pd.read_csv("/content/combined_fasta_targets.csv")
print(df)

"""### faa file inspection"""

!pip install biopython

from Bio import SeqIO

fasta_path = "/content/protein.faa"

print("🧬 FASTA entries found in protein.faa:")
for record in SeqIO.parse(fasta_path, "fasta"):
    print(f"- {record.id} | {record.description}")

"""### Applying ProtTrans"""

from transformers import T5Tokenizer, T5EncoderModel
import torch
import numpy as np
import pandas as pd
from tqdm import tqdm
import pickle

# Load ProtTrans model and tokenizer
tokenizer = T5Tokenizer.from_pretrained("Rostlab/prot_t5_xl_uniref50", do_lower_case=False)
model_pt = T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_uniref50")
model_pt = model_pt.eval().to("cuda" if torch.cuda.is_available() else "cpu")

# Load sequences from CSV (skip header)
df = pd.read_csv("/content/combined_fasta_targets.csv")
sequences = df.iloc[:, 0].tolist()  # assumes first column contains sequences

# Embed each sequence
prot_embeddings = {}
for i, seq in tqdm(enumerate(sequences)):
    clean_seq = seq.replace(" ", "")
    spaced_seq = " ".join(list(clean_seq))
    ids = tokenizer(spaced_seq, return_tensors="pt", padding=True)
    ids = {k: v.to(model_pt.device) for k, v in ids.items()}
    with torch.no_grad():
        embedding = model_pt(**ids).last_hidden_state.mean(1).squeeze().cpu().numpy()
    prot_embeddings[f"prot_{i:02d}"] = embedding  # e.g., prot_00, prot_01, ...

# Save the embeddings
with open("/content/sars2_prottrans_embeddings.pkl", "wb") as f:
    pickle.dump(prot_embeddings, f)

print("✅ Embeddings saved as sars2_prottrans_embeddings.pkl")

"""## STEP 2: Applying PCA and make df_infer"""

import pickle
import numpy as np
import pandas as pd

# Step 1: Load embedded ProtTrans vectors (1024 → 512 after PCA)
with open("/content/sars2_prottrans_embeddings.pkl", "rb") as f:
    prot_embeddings = pickle.load(f)

with open("/content/pca_model.pkl", "rb") as f:
    pca_model = pickle.load(f)

# Step 2: Load the raw sequences (assumed in same order)
raw_sequences = pd.read_csv("/content/combined_fasta_targets.csv").iloc[:, 0].str.replace(" ", "").tolist()

# Step 3: Match the embeddings and apply PCA
sorted_keys = sorted(prot_embeddings.keys())  # should match prot_00 to prot_28
raw_vectors = np.stack([prot_embeddings[k] for k in sorted_keys])
reduced_vectors = pca_model.transform(raw_vectors)

# Step 4: Build aligned df
df_infer = pd.DataFrame({
    "Target Sequence": raw_sequences,
    "ProtTrans": list(reduced_vectors)
})

# Step 5: Save
df_infer.to_csv("/content/df_infer.csv", index=False)
print("✅ df_infer.csv saved with real Target Sequences and reduced embeddings.")

df_infer.head(1)

"""## Predict (MAKE SURE TO LOAD MY MODEL)

### Reformat
"""

import re
import json
import pandas as pd

# Step 1: Load as raw text
df_infer = pd.read_csv("/content/df_infer.csv")

# Step 2: Fix each malformed array string (regex inserts commas)
def fix_array_string(s):
    s = re.sub(r'(?<=\d)\s+(?=[\-]?\d)', ', ', s.strip())  # insert comma between numbers
    return json.dumps(eval(s))  # now it's safe to eval and dump to clean JSON

# Apply only to strings
df_infer["ProtTrans"] = df_infer["ProtTrans"].apply(lambda x: fix_array_string(x) if isinstance(x, str) else json.dumps(x))

# Save to a new file
df_infer.to_csv("/content/df_infer_fixed.csv", index=False)
print("✅ Saved: df_infer_fixed.csv (JSON-safe format)")

"""### Prediction"""

from DeepPurpose import utils
import pandas as pd
import numpy as np
import json
from DeepPurpose import CompoundPred
#model = CompoundPred.model_initialize(**config) # <- Untrained model
import pickle

with open("scaler.pkl", "rb") as f:
    scaler = pickle.load(f)


# 🔥 Load your FASTA-to-ProteinName mapping
df_mapping = pd.read_csv("/content/metrics - SARS2 FASTA.csv")

# 🔥 Load df_infer and merge with names
df_infer = pd.read_csv("/content/df_infer_fixed.csv", converters={"ProtTrans": json.loads})
df_infer = df_infer.merge(df_mapping, on="Target Sequence", how="left")


# Ligands to test
ligand_smiles_list = {
    "13-cis-Retinoic Acid": r"CC1=C(C(CCC1)(C)C)/C=C/C(=C/C=C/C(=C\C(=O)O)/C)/C",
    "4-Hydroxybenzaldehyde": r"C1=CC(=CC=C1C=O)O",
    "8-Prenylnaringenin": r"CC(=CCC1=C2C(=C(C=C1O)O)C(=O)C[C@H](O2)C3=CC=C(C=C3)O)C",
    "All-trans Retinoic Acid": r"CC1=C(C(CCC1)(C)C)/C=C/C(=C/C=C/C(=C/C(=O)O)/C)/C",
    "Apigenin": r"C1=CC(=CC=C1C2=CC(=O)C3=C(C=C(C=C3O2)O)O)O",
    "Artemisinin": r"C[C@@H]1CC[C@H]2[C@H](C(=O)O[C@H]3[C@@]24[C@H]1CC[C@](O3)(OO4)C)C",
    "Auraptene": r"CC(=CCC/C(=C/COC1=CC2=C(C=C1)C=CC(=O)O2)/C)C",
    "Avermectin B1a": r"CC[C@H](C)[C@@H]1[C@H](C=C[C@@]2(O1)C[C@@H]3C[C@H](O2)C/C=C(/[C@H]([C@H](/C=C/C=C/4\CO[C@H]5[C@@]4([C@@H](C=C([C@H]5O)C)C(=O)O3)O)C)O[C@H]6C[C@@H]([C@H]([C@@H](O6)C)O[C@H]7C[C@@H]([C@H]([C@@H](O7)C)O)OC)OC)\C)C",
    "Beta-glucan": r"C([C@@H]1[C@H]([C@@H]([C@H]([C@@H](O1)OC2[C@H](O[C@H]([C@@H]([C@H]2O)O)OC3[C@H](O[C@@H]([C@@H]([C@H]3O)O)O)CO)CO)O)O)O)O)O",
    "Bryostatin-1": r"CCC/C=C/C=C/C(=O)O[C@H]1/C(=C/C(=O)OC)/C[C@H]2C[C@@H](OC(=O)C[C@@H](C[C@@H]3C[C@@H](C([C@@](O3)(C[C@@H]4C/C(=C/C(=O)OC)/C[C@@H](O4)/C=C/C([C@@]1(O2)O)(C)C)O)(C)C)OC(=O)C)O)[C@@H](C)O",
    "Butyrate": r"CCCC(=O)[O-]",
    "Caffeic acid": r"C1=CC(=C(C=C1/C=C/C(=O)O)O)O",
    "Caffeic acid phenethyl ester": r"C1=CC=C(C=C1)CCOC(=O)/C=C/C2=CC(=C(C=C2)O)O",
    "Catechin": r"C1[C@@H]([C@H](OC2=CC(=CC(=C21)O)O)C3=CC(=C(C=C3)O)O)O",
    "Celastrol": r"CC1=C(C(=O)C=C2C1=CC=C3[C@]2(CC[C@@]4([C@@]3(CC[C@@]5([C@H]4C[C@](CC5)(C)C(=O)O)C)C)C)C)O",
    "Coumarin": r"C1=CC=C2C(=C1)C=CC(=O)O2",
    "Eckol": r"C1=C(C=C(C=C1O)OC2=C(C=C(C3=C2OC4=C(C=C(C=C4O3)O)O)O)O)O",
    "Ergothioneine": r"C[N+](C)(C)[C@@H](CC1=CNC(=S)N1)C(=O)[O-]",
    "Erythronolide B": r"CC[C@@H]1[C@@H]([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O)C)O)(C)O)C)C)O)C",
    "Eugenol": r"COC1=C(C=CC(=C1)CC=C)O",
    "Ferulic acid": r"COC1=C(C=CC(=C1)/C=C/C(=O)O)O",
    "Fisetin": r"C1=CC(=C(C=C1C2=C(C(=O)C3=C(O2)C=C(C=C3)O)O)O)O",
    "Fucoxanthinol": r"C/C(=C\\C=C\\C=C(/C)\\C=C\\C=C(/C)\\C(=O)C[C@]12[C@](O1)(C[C@H](CC2(C)C)O)C)/C=C/C=C(\\C)/C=C=C3[C@](C[C@H](CC3(C)C)O)(C)O",
    "Gallic acid": r"C1=C(C=C(C(=C1O)O)O)C(=O)O",
    "Guanosine": r"C1=NC2=C(N1[C@H]3[C@@H]([C@@H]([C@H](O3)CO)O)O)N=C(NC2=O)N",
    "Guaiacol": r"COC1=CC=CC=C1O",
    "Halocynthiaxanthin": r"CC1=C(C(C[C@@H](C1)O)(C)C)C#C/C(=C/C=C/C(=C/C=C/C=C(\C)/C=C/C=C(\C)/C(=O)C[C@]23[C@](O2)(C[C@H](CC3(C)C)O)C)/C)/C",
    "Hesperidin": r"C[C@H]1[C@@H]([C@H]([C@H]([C@@H](O1)OC[C@@H]2[C@H]([C@@H]([C@H]([C@@H](O2)OC3=CC(=C4C(=O)C[C@H](OC4=C3)C5=CC(=C(C=C5)OC)O)O)O)O)O)O)O)O",
    "Hydroxyproline": r"C1[C@H](CN[C@@H]1C(=O)O)O",
    "Licoflavone C": r"CC(=CCC1=C2C(=C(C=C1O)O)C(=O)C=C(O2)C3=CC=C(C=C3)O)C",
    "Luteolin": r"C1=CC(=C(C=C1C2=CC(=O)C3=C(C=C(C=C3O2)O)O)O)O",
    "Lycopene": r"CC(=CCC/C(=C/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C=C(/CCC=C(C)C)\C)\C)\C)/C)/C)C)C",
    "Melatonin": r"CC(=O)NCCC1=CNC2=C1C=C(C=C2)OC",
    "Monensin": r"CC[C@]1(CC[C@@H](O1)[C@@]2(CC[C@@]3(O2)C[C@@H]([C@H]([C@H](O3)[C@@H](C)[C@H]([C@H](C)C(=O)O)OC)C)O)C)[C@H]4[C@H](C[C@@H](O4)[C@@H]5[C@H](C[C@H]([C@@](O5)(CO)O)C)C)C",
    "Myricetin": r"C1=C(C=C(C(=C1O)O)O)C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O",
    "Myriocin": r"CCCCCCC(=O)CCCCCC/C=C/C[C@H]([C@@H]([C@@](CO)(C(=O)O)N)O)O",
    "Nigericin": r"C[C@H]1CC[C@@H](O[C@H]1[C@@H](C)C(=O)O)C[C@@H]2C[C@H]([C@H]([C@@]3(O2)[C@@H](C[C@@](O3)(C)[C@H]4CC[C@@](O4)(C)[C@H]5[C@H](C[C@@H](O5)[C@@H]6[C@H](C[C@H]([C@@](O6)(CO)O)C)C)C)C)C)OC",
    "Niphimycin": r"CC1CCC(C(C(CC(C(/C=C/C(C(C(=O)O[C@@H](C(/C=C/C=C/C(C(C(CC(CC(CC2CC(C(C(O2)(CC1O)O)O)O)OC(=O)CC(=O)O)O)O)C)O)C)[C@@H](C)C[C@@H](C)CCC/C=C/CCCNC(=NC)N)C)O)C)O)O)C)O",
    "Phloroglucinol": r"C1=C(C=C(C=C1O)O)O",
    "Plakortide F Acid": r"CC/C=C/C(CC)CCC[C@]1(C[C@@H]([C@@H](OO1)CC(=O)O)CC)CC",
    "Pterostilbene": r"COC1=CC(=CC(=C1)/C=C/C2=CC=C(C=C2)O)OC",
    "Resorcinol": r"C1=CC(=CC(=C1)O)O",
    "Resveratrol": r"C1=CC(=CC=C1/C=C/C2=CC(=CC(=C2)O)O)O",
    "Shikonin": r"CC(=CC[C@H](C1=CC(=O)C2=C(C=CC(=C2C1=O)O)O)O)C",
    "Sophoraflavanone G": r"CC(=CC[C@H](CC1=C2C(=C(C=C1O)O)C(=O)C[C@H](O2)C3=C(C=C(C=C3)O)O)C(=C)C)C",
    "Uridine": r"C1=CN(C(=O)NC1=O)[C@H]2[C@@H]([C@@H]([C@H](O2)CO)O)O",
    "Vanillin": r"COC1=C(C=CC(=C1)C=O)O",
    "Withaferin A": r"CC1=C(C(=O)O[C@H](C1)[C@@H](C)[C@H]2CC[C@@H]3[C@@]2(CC[C@H]4[C@H]3C[C@@H]5[C@]6([C@@]4(C(=O)C=C[C@@H]6O)C)O5)C)CO",
    "Withanolide A": r"CC1=C(C(=O)O[C@H](C1)[C@@](C)([C@H]2CC[C@@H]3[C@@]2(CC[C@H]4[C@H]3[C@H]5[C@H](O5)[C@@]6([C@@]4(C(=O)C=CC6)C)O)C)O)C",
    "Xanthohumol": r"CC(=CCC1=C(C(=C(C=C1O)OC)C(=O)/C=C/C2=CC=C(C=C2)O)O)C",
    "Zerumbone": r"C/C/1=C\\CC(/C=C/C(=O)/C(=C/CC1)/C)(C)C"
}


# Step 1: Repeat df_infer once per ligand
df_list = []
for ligand_name, smiles in ligand_smiles_list.items():
    df_ligand = df_infer.copy()
    df_ligand["SMILES"] = smiles
    df_ligand["Label"] = 0  # dummy value
    df_ligand["Ligand"] = ligand_name
    df_list.append(df_ligand)

df = pd.concat(df_list).reset_index(drop=True)


# Step 2: Run DeepPurpose-compatible preprocessing
processed = utils.data_process(
    drug_encoding = "Transformer",
    target_encoding = "ProtTrans",
    X_drug = df["SMILES"].values,
    X_target = df["Target Sequence"].values,  # dummy string field, required
    y = df["Label"].values,
    split_method = "no_split",
    X_target_ = df["ProtTrans"].tolist(),  # actual embedded vectors go here
    mode = "DTI"
)


# Step 3: Predict pIC50
pred = model.predict(processed)

# Step 4: Inverse transform (back to real pIC50/IC50 units)
df["pIC50_scaled"] = pred
df["pIC50"] = scaler.inverse_transform(np.array(pred).reshape(-1, 1))
df["IC50_nM"] = 10 ** (9 - df["pIC50"])

# Step 5: Tidy result
result = df[["Ligand", "Protein Names", "pIC50", "IC50_nM"]].copy()

# ✅ Show top results
print(result.head(50))

# --- collapse the 29-protein fan-out to a single row per ligand -------------
# (they all have identical pIC50 / IC50_nM, so .first() is fine)
ligand_summary = (
    result
    .groupby("Ligand", as_index=False)
    .agg({"pIC50": "first", "IC50_nM": "first"})
)

# --- pick the 20 most potent ligands (highest pIC50 == lowest IC50) --------
top20 = (
    ligand_summary
    .sort_values("pIC50", ascending=False)   # or "IC50_nM", ascending=True
    .head(20)
    .reset_index(drop=True)
)

print(top20)

result.to_csv("/content/deeppurpose_result.csv")

print("📉 Prediction StdDev:", np.std(pred))

"""## Re-clone new repo"""

# 🧹 Remove existing repo if it exists
!rm -rf Deeppurpose

# 📥 Clone fresh copy
!git clone https://github.com/mosmos6/Deeppurpose

# ⚙️ Reinstall the package
!pip install ./Deeppurpose

"""# SMILES Generator"""

!pip install pubchempy

from pubchempy import get_compounds

compound_names = [
    "13-cis-Retinoic Acid",
    "4-Hydroxybenzaldehyde",
    "8-Prenylnaringenin",
    "All-trans Retinoic Acid",
    "Amphidinolide",
    "Apigenin",
    "Artemisinin",
    "Astaxanthin",
    "Auraptene",
    "Avermectin B1a",
    "Beta-glucan",
    "Bisphenol A",
    "Bisresorcinol",
    "Bromophenols",
    "Bryostatin-1",
    "Butyrate",
    "Caffeic acid",
    "Caffeic acid phenethyl ester",
    "Catechin",
    "Celastrol",
    "Conotoxin peptide fragments",
    "Coumarin",
    "Curculigosides",
    "Desmethylol-karitegeisin",
    "Eckol",
    "Epothilone B",
    "Ergothioneine",
    "Erythronolide B",
    "Eugenol",
    "Ferulic acid",
    "Fisetin",
    "Fish mucus peptides",
    "Fucoidan",
    "Fucoxanthinol",
    "Gallic acid",
    "Gracilamide",
    "Guanosine",
    "Guaiacol",
    "Halocynthiaxanthin",
    "Hesperidin",
    "Hydroxyproline",
    "Isobavachalcone",
    "Laminarin",
    "Laeva D",
    "Licoflavone C",
    "Luteolin",
    "Lycopene",
    "Melatonin",
    "Methylhentriol",
    "Monensin",
    "Myricetin",
    "Myriocin",
    "Nigericin",
    "Niphimycin",
    "Pentacyclic Triterpenoids",
    "Phloroglucinol",
    "Phloroglucinol derivatives",
    "Plakortide F Acid",
    "Polyether ionophores",
    "Pterostilbene",
    "Reupefek F",
    "Resorcinol",
    "Resveratrol",
    "Rhodiola rosea",
    "Rosavin",
    "Shikonin",
    "Sophoraflavanone G",
    "Spinosyn A",
    "Spirulina platensis polyphenols",
    "Sulfated polysaccharides",
    "Syringaldehyde",
    "Tetrahydrocurcumin",
    "Thymol",
    "Tylactone derivatives",
    "Tylosin",
    "Ulva lactuca",
    "Uridine",
    "Vanillin",
    "Withaferin A",
    "Withanolide A",
    "Xanthohumol",
    "Zerumbone"
]


smiles_results = {}
failed = []

for name in compound_names:
    try:
        result = get_compounds(name, 'name')
        if result and result[0].isomeric_smiles:
            smiles_results[name] = result[0].isomeric_smiles
        else:
            failed.append(name)
    except Exception as e:
        failed.append(name)

# 結果表示
print("✅ 成功した化合物数:", len(smiles_results))
print("❌ 失敗した化合物数:", len(failed))
print("\n--- 成功したSMILES ---")
for name, smiles in smiles_results.items():
    print(f"{name}: {smiles}")

if failed:
    print("\n--- 見つからなかった化合物 ---")
    for name in failed:
        print(name)





"""# Store all the file outputs

## Save files
"""

import os
import shutil

# ✅ Define the destination folder in Google Drive
save_dir = "/content/drive/My Drive/MD_Simulation/Withaferin A_VS_ORF9b"

# ✅ Create the directory if it doesn't exist
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

# ✅ Get a list of all files in the current working directory
all_files = [f for f in os.listdir("/content") if os.path.isfile(f"/content/{f}")]

# ✅ Save all detected files to Google Drive
for file in all_files:
    src = f"/content/{file}"
    dst = f"{save_dir}/{file}"
    shutil.copy(src, dst)
    print(f"✅ Saved {file} to Google Drive!")

"""## Load files"""

import os
import shutil

save_dir = "/content/drive/My Drive/MD_Simulation/Withaferin A_VS_ORF9b"

# ✅ Get a list of all saved files in Google Drive
saved_files = os.listdir(save_dir)

# ✅ Copy them back to Colab
for file in saved_files:
    src = f"{save_dir}/{file}"
    dst = f"/content/{file}"
    shutil.copy(src, dst)
    print(f"✅ Restored {file} from Google Drive!")

