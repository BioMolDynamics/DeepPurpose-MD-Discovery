### 3_docking_vina.py
# Author: Iori Mochizuki
# Created: 2025-05-13
# Description: Run docking using AutoDock Vina and extract best pose

"""
Step 3: Docking with AutoDock Vina
- Calculate centroid (box center) from receptor structure
- Run AutoDock Vina with specified box center and size
- Extract the best scoring pose (MODEL 1)
- Save to clean PDB file for MD alignment
"""

import subprocess

# Step 1: Calculate centroid from receptor ATOM entries
x_sum = y_sum = z_sum = 0.0
atom_count = 0

with open("receptor_clean.pdb", "r") as f:
    for line in f:
        if line.startswith("ATOM"):
            x_sum += float(line[30:38].strip())
            y_sum += float(line[38:46].strip())
            z_sum += float(line[46:54].strip())
            atom_count += 1

if atom_count == 0:
    raise ValueError("No ATOM entries found in receptor file.")

x_center = round(x_sum / atom_count, 3)
y_center = round(y_sum / atom_count, 3)
z_center = round(z_sum / atom_count, 3)

print(f"📍 Docking box center: x={x_center}, y={y_center}, z={z_center}")

# Step 2: Run AutoDock Vina
vina_command = f"""
./vina_1.2.5_linux_x86_64 \
--receptor receptor.pdbqt \
--ligand ligand.pdbqt \
--out output.pdbqt \
--center_x {x_center} \
--center_y {y_center} \
--center_z {z_center} \
--size_x 30 --size_y 30 --size_z 30 \
--seed 12345 --exhaustiveness 20
"""
print("Running AutoDock Vina...")
subprocess.run(vina_command, shell=True, check=True)

# Step 3: Extract MODEL 1 as best pose
with open("output.pdbqt", "r") as infile, open("best_docked_ligand.pdb", "w") as outfile:
    inside_model = False
    for line in infile:
        if line.startswith("MODEL 2"):
            break
        if line.startswith("MODEL 1"):
            inside_model = True
        if inside_model and not line.startswith("MODEL"):
            outfile.write(line)

print("✅ Best docking pose extracted → best_docked_ligand.pdb")

# Step 4: Clean up docking-specific columns and write final PDB
with open("best_docked_ligand.pdb", "r") as infile, open("output.pdb", "w") as outfile:
    for line in infile:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            outfile.write(line)

print("✅ Cleaned ligand pose written to output.pdb")
