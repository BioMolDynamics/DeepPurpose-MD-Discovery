# 🧬 DeepPurpose-MD-Discovery

A streamlined, Colab-optimized drug discovery pipeline integrating:

- ✅ Ligand-target prediction with a custom-trained [DeepPurpose](https://github.com/mosmos6/Deeppurpose) fork
- ✅ Structural docking using [AutoDock Vina](http://vina.scripps.edu/)
- ✅ GPU-accelerated Molecular Dynamics with [OpenMM](https://openmm.org/) and [OpenFF](https://openforcefield.org/)
- ✅ Stability and mechanistic analyses via PCA, FEL, RMSD, H-bonding, water networks, and more.

This repo demonstrates how to simulate and evaluate ligand–protein interactions from end-to-end using SARS-CoV-2 viral proteins as case studies.

---

## 🔧 Setup Instructions (Google Colab)

This pipeline is designed for use in **Google Colab**, with full support for `condacolab`.

---

### ✅ Step 1: Enable Conda in Colab

Paste the following **at the very top of your Colab notebook**:

```python
!pip install -q condacolab
import condacolab
condacolab.install()```



🔄 NOTE: This will crash your runtime once. That's expected.

After Colab restarts, rerun the following cell:

```import condacolab
condacolab.check()
!mamba env update -n base -f environment.yml```


