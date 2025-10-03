### 🔄 Updates — 2025‑09 → 2025‑10 (v2.5)

**No breaking changes.** Existing notebooks and CLI flows continue to work.  
This release adds determinism, a quality‑of‑life docking control, and a safety fix in receptor prep.

1) Deterministic MD runs (`scripts/5_md_simulation.py`)
   - New flag: `--seed <int>`  
     Seeds both the Langevin integrator and the Monte Carlo barostat for reproducible trajectories.
   - Example:  
     ```bash
     conda run python scripts/5_md_simulation.py --seed 42
     ```
   - The chosen seed is echoed in the log.

2) Manual docking‑box center override (`scripts/3_docking_vina.py`)
   - New flag: `--box-center X Y Z` (in Å)  
     Allows you to override the automatically detected box center for Vina.

3) Duplicate‑atom cleanup in receptor prep (`scripts/2_prepare_receptor.py`)
   - Detects and removes **exact duplicate atoms** in input PDBs by `--resolve-altdups`
     (common artifact from some structure writers; duplicates share the same residue/name/element
     with indistinguishable coordinates).  
   - Prints a small summary and preserves indexing for unique atoms.  
   - Prevents downstream selection/pipeline mismatches.

**Also in this train**
- Minor doc touch‑ups and clearer runtime messages.  
- The core end‑to‑end usage (DeepPurpose → Vina → OpenMM) is unchanged.



# 🧬 DeepPurpose-MD-Discovery

A streamlined, Colab-optimized drug discovery pipeline integrating:

- ✅ Ligand-target prediction with a custom-trained [DeepPurpose](https://github.com/BioMolDynamics/Deeppurpose) fork
- ✅ Structural docking using [AutoDock Vina](http://vina.scripps.edu/)
- ✅ GPU-accelerated Molecular Dynamics with [OpenMM](https://openmm.org/) and [OpenFF](https://openforcefield.org/)
- ✅ **RNA-enabled**: analyze protein and RNA as targets, including viral subgenomic RNA (SARS-CoV-2 case studies)
- ✅ Mechanistic analyses: **PCA, FEL, RMSD, H-bonds, water networks, allostery (Δexposure), π–π stacking**, and more—now via an interactive Colab UI.

This repository demonstrates simulation and evaluation of ligand–protein/RNA interactions, with robust, **user-guided exploratory analysis in Google Colab.**

---

## 🔧 Setup Instructions (Google Colab)

This pipeline is designed for use in **Google Colab**, with full support for `condacolab`. A demo notebook is included in this repository to reproduce all steps.
All advanced trajectory analyses—including RNA-specific and allosteric functions—are performed via a **UI-driven Colab panel** (see below).

---

### ✅ Step 1: Enable Conda in Colab

Paste the following **at the very top of your Colab notebook**:

```python
!pip install -q condacolab
import condacolab
condacolab.install()
```



🔄 NOTE: This will crash your runtime once. That's expected.

After Colab restarts, rerun the following cell:

```python
import condacolab
condacolab.check()
```

### ✅ Step 2: Clone Required Repositories

```python
# Main pipeline repo (this one)
!git clone https://github.com/BioMolDynamics/DeepPurpose-MD-Discovery.git

# Custom fork of DeepPurpose (installed later)
!git clone https://github.com/BioMolDynamics/Deeppurpose
```

### ✅ Step 3: Download AutoDock Vina Binary

```python
!wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64
!chmod +x vina_1.2.5_linux_x86_64
!./vina_1.2.5_linux_x86_64 --version
```

### ✅ Step 4: Install Conda Environment

```python
!mamba env create -f environment.yml
```

### ✅ Step 5: Finalize Setup

```python
# Install custom fork of DeepPurpose without overwriting key dependencies
!conda run -n deeppurpose-md-env pip install --no-deps ./Deeppurpose

# Install optional dependencies like Open Babel
!conda run -n deeppurpose-md-env python scripts/install_optional.py
```

### ✅ Step 6: Run the Pipeline

```python
!conda run -n deeppurpose-md-env python scripts/1_prepare_ligand.py "$ligand_smiles"
!conda run -n deeppurpose-md-env python scripts/2_prepare_receptor.py "$pdb_id" --strict-protein  # or --rna for RNA
!conda run -n deeppurpose-md-env python scripts/3_docking_vina.py --use-residue-centroid
!conda run -n deeppurpose-md-env python scripts/3b_prepare_protein.py  # For proteins
# RNA users: follow README/Colab for RNA-specific prep before MD
!conda run -n deeppurpose-md-env python scripts/4_align_ligand.py
!conda run -n deeppurpose-md-env python scripts/5_md_simulation.py --protein/--rna [--no-ligand]  # See notebook for options
!conda run -n deeppurpose-md-env python scripts/7_deeppurpose_training.py
!conda run -n deeppurpose-md-env python scripts/8_deeppurpose_prediction.py
```
Each script corresponds to a specific stage in the full drug discovery pipeline — from ligand design to MD simulation to deep learning prediction.

**Script 5b (old analysis script) is deprecated.
Script 6 is now performed directly via a UI analysis panel in Colab, using interactive checkboxes to launch all trajectory analyses (RMSD, PCA, H-bonds, RNA allostery, π–π stacking, etc.) with no additional scripts required.**

### 🖥️ Interactive Analysis: UI-Driven in Colab
All trajectory analysis is now performed via a **UI panel in your Colab notebook**:
- Select and run only the analyses you need (checkbox UI).


- All backend code is visible (for transparency/extensibility), but users only interact with the UI.


- **RNA-specialized analyses are available** (Watson–Crick pairs, backbone Δexposure, ligand–RNA atom contacts, π–π stacking, and more).


- Outputs are auto-saved as images, CSVs, and PDBs for downstream reporting or 3D visualization.


(See the demo Colab notebook for details and examples.)





### 🧪 Dataset Information

This demo uses a COVID-19 specific subset of BindingDB, available from UC San Diego.  
To simplify setup, we provide pre-cleaned versions of this dataset:

- `BindingDB_Covid-19.tsv` (214MB, hosted via Zenodo)
- `strong_binders_cleaned.csv` (optional for filtering)
- `protein.faa` (optional for filtering)
- `metrics - SARS2 FASTA.csv` (matching data of SARS-CoV-2 proteins and FASTA)

📎 Dataset download link: [https://doi.org/10.5281/zenodo.15613825)

You are welcome to use your own SMILES/FASTA data by modifying `7_deeppurpose_training.py`.

### 🧬 Features
- **Fully RNA-capable**: Run all stages (including MD, contact analysis, allostery, π–π stacking) for RNA targets.


- Interactive, modular trajectory analysis: UI lets you select any combination of analyses, including custom, RNA-specific ones.


- Robust error handling: clear feedback if input files/definitions are missing.


- All outputs are automatically saved (figures, CSV, PDB).


- **Colab-native**: No installation required outside the notebook.

### 📜 License
MIT License. Please cite this repository if used in academic work.

### 📖 Citation
If you use **DeepPurpose-MD** in your work, please cite:

Mochizuki, I. (2025). DeepPurpose-MD: An End-to-End Colab-Based Drug Discovery Pipeline Integrating Docking, Molecular Dynamics, and Deep Learning. Zenodo. https://doi.org/10.5281/zenodo.15613825

### 📫 Contact
Maintained by BioMolDynamics
For academic inquiries, collaboration, or feedback, please open an issue.


