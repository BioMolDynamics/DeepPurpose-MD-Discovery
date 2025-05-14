### 7_deeppurpose_training.py
# Author: Iori Mochizuki
# Step 7: Dataset Cleaning and Preparation for DeepPurpose Training

import os
import pandas as pd
from DeepPurpose.dataset import load_bindingdb_covid_tsv

# ✅ Paths
raw_tsv_path = "BindingDB_Covid-19.tsv"  # Make sure this file is uploaded in advance
cleaned_tsv_path = "cleaned_bindingdb.tsv"
processed_tsv_path = "processed_bindingdb.tsv"

# ✅ Step 1: Clean Header Duplications
with open(raw_tsv_path, 'r', encoding='utf-8', errors='ignore') as f:
    raw_lines = f.readlines()

header_line = raw_lines[0].strip()
cleaned_lines = [raw_lines[0]] + [line for line in raw_lines[1:] if line.strip() != header_line]

with open(cleaned_tsv_path, 'w', encoding='utf-8') as f:
    f.writelines(cleaned_lines)

print(f"✅ Cleaned TSV saved as {cleaned_tsv_path}")

# ✅ Step 2: Select Required Columns and Clean Values
usecols = ['Ligand SMILES', 'BindingDB Target Chain Sequence', 'IC50 (nM)']
df = pd.read_csv(cleaned_tsv_path, sep='\t', usecols=usecols, engine='python')
df = df.dropna()

# ✅ Clean SMILES and Sequence
df = df[df['Ligand SMILES'].apply(lambda x: isinstance(x, str) and len(x) > 0)]
df = df[df['BindingDB Target Chain Sequence'].apply(lambda x: isinstance(x, str) and len(x) > 0)]

# ✅ Clean numeric IC50 values and convert to molar
df['IC50 (nM)'] = df['IC50 (nM)'].astype(str).str.extract(r'([\d.]+)').astype(float)
df = df.dropna()
df['Affinity'] = df['IC50 (nM)'] * 1e-9  # convert to molar units

df = df[['Ligand SMILES', 'BindingDB Target Chain Sequence', 'Affinity']]
df.columns = ['SMILES', 'Target Sequence', 'Affinity']
df.to_csv(processed_tsv_path, sep='\t', index=False)

print(f"✅ Cleaned and saved {len(df)} compound-target pairs to {processed_tsv_path}")

# ✅ Step 3: Load into DeepPurpose Dataset Format
X = load_bindingdb_covid_tsv(processed_tsv_path)
