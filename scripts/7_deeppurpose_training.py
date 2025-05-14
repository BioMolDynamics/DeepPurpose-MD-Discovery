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


### 7_deeppurpose_training.py (Part 2)
# Author: Iori Mochizuki
# Step 7b: ProtTrans Embedding and Data Augmentation

import pandas as pd
import numpy as np
import torch
import pickle
from tqdm import tqdm
from transformers import T5Tokenizer, T5EncoderModel
import matplotlib.pyplot as plt

# === Load processed dataset ===
df = pd.read_csv("processed_bindingdb.tsv", sep='\t')
unique_targets = df['Target Sequence'].unique()
print(f"🔎 Unique protein sequences: {len(unique_targets)}")

# === Load ProtTrans model ===
tokenizer = T5Tokenizer.from_pretrained("Rostlab/prot_t5_xl_uniref50", do_lower_case=False)
model = T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_uniref50").eval()

# === Move model to GPU if available ===
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model = model.to(device)

# === Generate ProtTrans embeddings ===
protein_embeddings = {}

with torch.no_grad():
    for seq in tqdm(unique_targets):
        seq_clean = seq.replace(" ", "")
        tokenized = tokenizer(" ".join(list(seq_clean)), return_tensors="pt", padding=True).to(device)
        embedding = model(**tokenized).last_hidden_state
        pooled = torch.mean(embedding, dim=1)
        protein_embeddings[seq_clean] = pooled.squeeze().cpu().numpy()

# === Save embeddings to disk ===
with open("protein_embeddings.pkl", "wb") as f:
    pickle.dump(protein_embeddings, f)
print("✅ Embeddings saved to protein_embeddings.pkl")

# === Merge embeddings with dataframe ===
df["ProtTrans"] = df["Target Sequence"].map(lambda s: protein_embeddings.get(s.replace(" ", ""), None))
df.to_pickle("embedded_bindingdb.pkl")
print("✅ Merged ProtTrans embeddings with dataset → embedded_bindingdb.pkl")

# === Load filtered strong binders (manually prepared) ===
df_strong = pd.read_csv("strong_binders_cleaned.csv")
top_targets = df_strong["Target Sequence"].value_counts()

# === Save and plot top proteins ===
top_targets.to_csv("top_protein_targets.csv", header=["Count"])
print("✅ Saved protein frequency table to top_protein_targets.csv")

plt.figure(figsize=(10, 4))
top_targets.head(10).plot(kind="barh", title="Top 10 Protein Targets")
plt.xlabel("Count")
plt.gca().invert_yaxis()
plt.grid(True)
plt.tight_layout()
plt.show()

# === Sample top binders: max 150 ligands per protein ===
df_sampled = (
    df_strong.groupby("Target Sequence", group_keys=False)
    .apply(lambda g: g.sample(n=min(len(g), 150), random_state=42))
)

df_sampled.to_csv("strong_binders_top150_per_protein.csv", index=False)
print(f"✅ Final dataset shape: {df_sampled.shape}")
