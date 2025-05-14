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
