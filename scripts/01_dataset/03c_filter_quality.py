# Verified script: copied unchanged from scripts_logs/filter_acb.py (header comment added).
# Result on the study data: 897 pass (identical to highquality_accessions.txt), 359 fail (1,256 candidates).
# Accessions of the 897 genomes: column "Name" of acb_highquality.tsv.
import pandas as pd

INPUT       = "output_assacineto/merged_final.tsv"
OUTPUT_PASS = "output_assacineto/acb_highquality.tsv"
OUTPUT_FAIL = "output_assacineto/acb_lowquality.tsv"

# Load table
df = pd.read_csv(INPUT, sep="\t")
print(f"Total genomes loaded: {len(df)}")

# Calculate total N's from rate
df["Total_Ns"] = (df["# N's per 100 kbp"] / 100000) * df["Total length"]
df["Total_Ns"] = df["Total_Ns"].round(0).astype(int)

# Define filters
mask_completeness  = df["Completeness"] >= 95
mask_contamination = df["Contamination"] <= 5
mask_contigs       = df["Total_Contigs"] <= 300
mask_ns            = df["Total_Ns"] <= 500

# Combined filter
mask_pass = mask_completeness & mask_contamination & mask_contigs & mask_ns

# Split
df_pass = df[mask_pass].copy()
df_fail = df[~mask_pass].copy()

# Save
df_pass.to_csv(OUTPUT_PASS, sep="\t", index=False)
df_fail.to_csv(OUTPUT_FAIL, sep="\t", index=False)

print(f"High quality genomes: {len(df_pass)}")
print(f"Low quality genomes:  {len(df_fail)}")
