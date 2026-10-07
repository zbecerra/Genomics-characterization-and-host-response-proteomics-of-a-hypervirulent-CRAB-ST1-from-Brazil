# Verified script: copied unchanged from scripts_logs/merge_quast_checkm2.py (header comment added).
# Run from the project root. Merges QUAST, CheckM2 and NCBI metadata (metadata_acb.tsv) -> merged_final.tsv
import pandas as pd

# Load files
quast = pd.read_csv("output_assacineto/transposed_report.tsv", sep="\t")
checkm2 = pd.read_csv("output_assacineto/checkm2_output/quality_report.tsv", sep="\t")
metadata = pd.read_csv("metadata_acb.tsv", sep="\t")

# Rename key columns to match
quast = quast.rename(columns={"Assembly": "Name"})
metadata = metadata.rename(columns={"Assembly Accession": "Name"})

# Merge all three
merged = pd.merge(quast, checkm2, on="Name", how="inner")
merged = pd.merge(merged, metadata, on="Name", how="left")

# Save
merged.to_csv("output_assacineto/merged_final.tsv", sep="\t", index=False)

print(f"QUAST genomes: {len(quast)}")
print(f"CheckM2 genomes: {len(checkm2)}")
print(f"Metadata genomes: {len(metadata)}")
print(f"Merged genomes: {len(merged)}")
print(f"Total columns: {len(merged.columns)}")
print("Saved to output_assacineto/merged_final.tsv")
