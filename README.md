# DeepPurpose-MD-Discovery Pipeline

This repository provides an end-to-end pipeline for identifying, docking, simulating, and analyzing potential small molecule inhibitors targeting SARS-CoV-2 viral proteins (e.g., ORF3a, ORF9b, Spike S1-NTD) using deep learning and molecular dynamics.

### Features
- 🧠 **DeepPurpose-based DTI prediction**, trained on COVID-specific BindingDB data (UC San Diego)
- 🔍 **AutoDock Vina** docking
- 🌊 **OpenMM GPU-accelerated MD simulations**
- 📈 PCA, FEL, RMSD/RMSF, hydrogen bonding, residence time, water shell
- 🔬 Deep ligand-protein interaction analysis (per-residue)

### Requirements

Install the dependencies using:

```bash
pip install -r requirements.txt
