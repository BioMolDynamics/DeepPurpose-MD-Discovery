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
