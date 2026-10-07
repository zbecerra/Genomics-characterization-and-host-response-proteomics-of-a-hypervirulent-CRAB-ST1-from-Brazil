# Reproducibility notes

## Output files

| File | Produced by | Used for |
|---|---|---|
| `ACB_ST1_pipeline/kept_ST1/*.fna`, `logs/mlst_results_all.tsv` | Step 1 | candidate ST1 genomes |
| `logs/reverify/mlst_reverify_full.tsv` | Step 2 | audit trail: 1256/1256 ST1 |
| `highquality_accessions.txt` | Step 3 | the 897 genomes (Table S2) |
| `results/pangenome/panaroo_results/core_gene_alignment.aln` | Step 5 | recombination, phylogeny |
| `*.per_branch_statistics.csv`, `recombination_summary.csv` | Step 6a | Gubbins r/m, % recombinant sites |
| `*.em.txt`, `*.importation_status.txt` | Step 6b | R/θ, δ, ν; imported blocks |
| `acb_ST1_snp_distance_matrix.tsv`, `*.treefile` | Step 6c | Figures 1–3, SNP distances |
| `results/amr/amrfinder_results/*.tsv` | Step 7 | Table S4 (AMR matrix) |
| `oxa23_context/oxa23_context_classification.tsv` | Step 8d | Table S6 (848 copies) |
| `databases/abar/AbaR_combined_v4.fna` | Step 8c | curated reference database |
| `vfdb_rebuilt_results/*_acb_vfdb_rebuilt.tsv` | Step 10 | Table S7 (virulome) |

## Troubleshooting

**`mlst` fails with "needs BLAST+ 2.9.0 or later" although BLAST+ 2.16 is installed.** `mlst` 2.17.6 compares the BLAST+ version as a decimal, so `2.16 >= 2.9` is false in Perl. Patch the installed script (keeps a backup):

```bash
bash scripts/00_setup/patch_mlst_blast_check.sh        # idempotent: patches once, keeps a .bak, does nothing if already patched
```

(Check first whether a newer `mlst` release has fixed the comparison; if so, upgrade instead of patching.)

The patch lives in the environment's `bin/mlst`; **cloning or renaming a conda environment does not carry it** — re-apply it afterwards. Perl/bioperl conflicts prevented upgrading the original environment, so a fresh environment with `mlst=2.17.6` was used.

**Truncated virulence hits (all sequences exactly 60 bp).** The first custom Acinetobacter VFDB databases were corrupted; always rebuild from the raw `VFDB_setB_nt.fas` (Step 10) and check `abricate --list` reports 966 sequences.

**Two copies of the merged virulence table.** Keep only `results/acb_ST1_complete_table_HQ_allVF.tsv`; a root-level copy from before the VFDB fix produced mismatches and was renamed `*.STALE_*`.

**AMRFinder output has no coordinates.** It was run on proteins; use Bakta locus tags to recover coordinates (Step 8a).

**BLAST classification is biased towards large references.** Never use whole chromosomes/plasmids as classification references; extract the isolated region (Step 8c) and use normalised coverage (Step 8d).

**Terminal output cut off (tmux scroll-back).** Redirect long outputs to a file (`> out.txt`) and inspect the file instead of pasting from the terminal.

**ClonalFrameML values look implausible (R/θ ~ 10⁻³).** The SNP-only alignment was used; rerun on the complete core-gene alignment (Step 6b).

## Key tool references

- Croucher NJ *et al.* (2015) Gubbins. *Nucleic Acids Res* 43:e15.
- Didelot X, Wilson DJ (2015) ClonalFrameML. *PLoS Comput Biol* 11:e1004041.
- Tonkin-Hill G *et al.* (2020) Panaroo. *Genome Biol* 21:180.
- Hoang DT *et al.* (2018) UFBoot2. *Mol Biol Evol* 35:518–522.
- Schwengers O *et al.* (2021) Bakta. *Microb Genom* 7:000685.
- Feldgarden M *et al.* (2021) AMRFinderPlus. *Sci Rep* 11:12728.
- Robertson J, Nash JHE (2018) MOB-suite. *Microb Genom* 4:e000206.
- Fisher RA (1922) On the interpretation of χ² from contingency tables, and the calculation of P. *J R Stat Soc* 85:87–94.

The complete reference list is in the manuscript.
