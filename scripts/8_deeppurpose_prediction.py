### 8_deeppurpose_prediction.py (Part 1)
# Author: Iori Mochizuki
# Created: 2025-05-13
# Description: Embed SARS-CoV-2 protein FASTA sequences for DeepPurpose prediction

import os
import pandas as pd

# === Step 1: Load and clean FASTA data ===

# Manually parse .faa file
fasta_path = "protein.faa"
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

# Save parsed FASTA to CSV
fasta_df = pd.DataFrame(sequences, columns=["Target Sequence"])
fasta_df.to_csv("parsed_fasta.csv", index=False)
print(f"✅ Parsed {len(fasta_df)} sequences from protein.faa")

# === Step 2: Combine with manually curated SARS-CoV-2 FASTA ===
manual_fasta = [
    "MESLVPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHLKDGTCGLVEVEKGVLPQLEQPYVFIKRSDARTAPHGHVMVELVAELEGIQYGRSGETLGVLVPHVGEIPVAYRKVLLRKNGNKGAGGHSYGADLKSFDLGDELGTDPYEDFQENWNTKHSSGVTRELMRELNGG",
    "AYTRYVDNNFCGPDGYPLECIKDLLARAGKASCTLSEQLDFIDTKRGVYCCREHEHEIAWYTERSEKSYELQTPFEIKLAKKFDTFNGECPNFVFPLNSIIKTIQPRVEKKKLDGFMGRIRSVYPVASPNECNQMCLSTLMKCDHCGETSWQTGDFVKATCEFCGTENLTKEGATTCGYLPQNAVVKIYCPACHNSEVGPEHSLAEYHNESGLKTILRKGGRTIAFGGCVFSYVGCHNKCAYWVPRASANIGCNHTGVVGEGSEGLNDNLLEILQKEKVNINIVGDFKLNEEIAIILASFSASTSAFVETVKGLDYKAFKQIVESCGNFKVTKGKAKKGAWNIGEQKSILSPLYAFASEAARVVRSIFSRTLETAQNSVRVLQKAAITILDGISQYSLRLIDAMMFTSDLATNNLVVMAYITGGVVQLTSQWLTNIFGTVYEKLKPVLDWLEEKFKEGVEFLRDGWEIVKFISTCACEIVGGQIVTCAKEIKESVQTFFKLVNKFLALCADSIIIGGAKLKALNLGETFVTHSKGLYRKCVKSREETGLLMPLKAPKEIIFLEGETLPTEVLTEEVVLKTGDLQPLEQPTSEAVEAPLVGTPVCINGLMLLEIKDTEKYCALAPNMMVTNNTFTLKGG",
    # Add more if needed
]
manual_df = pd.DataFrame(manual_fasta, columns=["Target Sequence"])

# Combine auto-parsed + manual
combined_df = pd.concat([fasta_df, manual_df], ignore_index=True)
combined_df = combined_df.drop_duplicates().reset_index(drop=True)

# Save unified FASTA
combined_df.to_csv("combined_fasta_targets.csv", index=False)
print(f"✅ Combined total: {len(combined_df)} sequences → combined_fasta_targets.csv")

# === Step 3: Inspect protein.faa entries (optional) ===
from Bio import SeqIO

print("🧬 FASTA entries found in protein.faa:")
for record in SeqIO.parse(fasta_path, "fasta"):
    print(f"- {record.id} | {record.description}")


### 8_deeppurpose_prediction.py (Part 2)
# Author: Iori Mochizuki
# Created: 2025-05-13
# Description: Apply ProtTrans and PCA to embed SARS-CoV-2 FASTA sequences

import pandas as pd
import numpy as np
import torch
import pickle
from transformers import T5Tokenizer, T5EncoderModel
from tqdm import tqdm

# === Step 1: Load ProtTrans model ===
print("🔄 Loading ProtTrans model and tokenizer...")
tokenizer = T5Tokenizer.from_pretrained("Rostlab/prot_t5_xl_uniref50", do_lower_case=False)
model_pt = T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_uniref50")
model_pt = model_pt.eval().to("cuda" if torch.cuda.is_available() else "cpu")
print("✅ ProtTrans model loaded.")

# === Step 2: Load FASTA sequences ===
df = pd.read_csv("combined_fasta_targets.csv")
sequences = df["Target Sequence"].tolist()

# === Step 3: Generate ProtTrans embeddings ===
prot_embeddings = {}
for i, seq in tqdm(enumerate(sequences), total=len(sequences), desc="🔬 Embedding FASTA"):
    clean_seq = seq.replace(" ", "")
    spaced_seq = " ".join(list(clean_seq))  # Insert spaces between amino acids
    ids = tokenizer(spaced_seq, return_tensors="pt", padding=True).to(model_pt.device)

    with torch.no_grad():
        emb = model_pt(**ids).last_hidden_state.mean(dim=1).squeeze().cpu().numpy()
        prot_embeddings[f"prot_{i:02d}"] = emb

# Save full 1024-dim embeddings
with open("sars2_prottrans_embeddings.pkl", "wb") as f:
    pickle.dump(prot_embeddings, f)
print("✅ ProtTrans embeddings saved → sars2_prottrans_embeddings.pkl")

# === Step 4: Apply PCA dimensionality reduction ===
with open("pca_model.pkl", "rb") as f:
    pca_model = pickle.load(f)

keys_sorted = sorted(prot_embeddings.keys())
matrix = np.stack([prot_embeddings[k] for k in keys_sorted])
reduced = pca_model.transform(matrix)

# === Step 5: Create final DataFrame for prediction ===
df_infer = pd.DataFrame({
    "Target Sequence": sequences,
    "ProtTrans": list(reduced)
})

df_infer.to_csv("df_infer.csv", index=False)
print("✅ Reduced embeddings saved → df_infer.csv")


### 8_deeppurpose_prediction.py (Part 3)
# Author: Iori Mochizuki
# Created: 2025-05-13
# Description: Run ligand screening against SARS-CoV-2 FASTA-embedded proteins using trained DeepPurpose model

import pandas as pd
import numpy as np
import json
import pickle
import re
from DeepPurpose import CompoundPred, utils

# === Step 1: Load trained DeepPurpose model ===
model_path = "deeppurpose_model_saved"
model = CompoundPred.model_pretrained(path_dir=model_path)
print("✅ DeepPurpose model loaded.")

# === Step 2: Fix JSON string format for ProtTrans column ===
df_infer = pd.read_csv("df_infer.csv")

def fix_array_string(s):
    s = re.sub(r'(?<=\d)\s+(?=[\-]?\d)', ', ', s.strip())
    return json.dumps(eval(s))

df_infer["ProtTrans"] = df_infer["ProtTrans"].apply(lambda x: fix_array_string(x) if isinstance(x, str) else json.dumps(x))
df_infer.to_csv("df_infer_fixed.csv", index=False)
print("✅ Cleaned ProtTrans array saved as df_infer_fixed.csv")

# === Step 3: Load scaler and protein names ===
with open("scaler.pkl", "rb") as f:
    scaler = pickle.load(f)

df_mapping = pd.read_csv("metrics - SARS2 FASTA.csv")  # Optional name map
df_infer = pd.read_csv("df_infer_fixed.csv", converters={"ProtTrans": json.loads})
df_infer = df_infer.merge(df_mapping, on="Target Sequence", how="left")

# === Step 4: Ligand SMILES list ===
ligand_smiles_list = {
    "Ferulic Acid": "COC1=C(C=CC(=C1)/C=C/C(=O)O)O",
    "Xanthohumol": "CC(=CCC1=C(C(=C(C=C1O)OC)C(=O)/C=C/C2=CC=C(C=C2)O)O)C",
    "Withaferin A": "CC1=C(C(=O)O[C@H](C1)[C@@H](C)[C@H]2CC[C@@H]3[C@@]2(CC[C@H]4[C@H]3C[C@@H]5[C@]6([C@@]4(C(=O)C=C[C@@H]6O)C)O5)C)CO",
    # Add others if needed
}

# === Step 5: Expand df_infer × ligand fanout ===
df_list = []
for ligand, smiles in ligand_smiles_list.items():
    temp = df_infer.copy()
    temp["Ligand"] = ligand
    temp["SMILES"] = smiles
    temp["Label"] = 0
    df_list.append(temp)

df = pd.concat(df_list).reset_index(drop=True)

# === Step 6: Prepare input for DeepPurpose inference ===
processed = utils.data_process(
    X_drug = df["SMILES"].values,
    X_target = df["Target Sequence"].values,
    y = df["Label"].values,
    drug_encoding = "Transformer",
    target_encoding = "ProtTrans",
    X_target_ = df["ProtTrans"].tolist(),
    split_method = "no_split",
    mode = "DTI"
)

# === Step 7: Predict pIC50 scores ===
pred = model.predict(processed)
df["pIC50_scaled"] = pred
df["pIC50"] = scaler.inverse_transform(np.array(pred).reshape(-1, 1))
df["IC50_nM"] = 10 ** (9 - df["pIC50"])

# === Step 8: Save full results ===
result = df[["Ligand", "Protein Names", "pIC50", "IC50_nM"]]
result.to_csv("deeppurpose_result.csv", index=False)
print("✅ Full prediction results saved to deeppurpose_result.csv")

# === Step 9: Collapse to top ligand scores ===
ligand_summary = result.groupby("Ligand", as_index=False).agg({"pIC50": "first", "IC50_nM": "first"})
top20 = ligand_summary.sort_values("pIC50", ascending=False).head(20).reset_index(drop=True)
print("🔝 Top ligands by predicted binding affinity:")
print(top20)
