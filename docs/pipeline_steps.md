# Step-by-step pipeline

## How to read this repository

Each step names the script(s) in `scripts/` that run it. Four kinds of commands appear below, and they are labelled explicitly:

- **Verified script** — a file in `scripts/`, run on the cluster during the study (or the corrected version of it), including all coordinates and thresholds.
- **Command from log** — the command line was recovered from a run log or the shell history of the study (the history only keeps the last 1,000 lines, so older commands are missing).
- **Reconstructed** — assembled from the tool's own log (version, threads, input names) and the Methods; flags not visible in a log are the tool defaults.
- **Command template** — the tool was run for the manuscript, but its exact command line was not recorded in a script that is part of this repository. The template reproduces the *parameters stated in the Methods*; paths and thread counts are placeholders. Check these against your own server logs before citing them as the exact command.


### Step 1 — Download the ST1 genome collection from NCBI (verified script)

`scripts/01_dataset/01_download_filter_ST1.sh` (run under `nohup`, resumable):

1. `datasets download genome taxon 470 --dehydrated` retrieves the *A. baumannii* manifest (Taxon ID 470, including genomes assigned to strain-level TaxIDs); all `GCA_/GCF_` accessions are extracted from `fetch.txt`.
2. **De-duplication:** where the same assembly is present as RefSeq (`GCF_`) and GenBank (`GCA_`), the **RefSeq accession is kept**; where several versions of an assembly exist (`.1`, `.2`, …), **only the latest version is kept**.
3. Each remaining accession is **downloaded one at a time** (`datasets download genome accession … --include genome`), typed with `mlst`, and the FASTA is **kept only if ST = 1** (Pasteur scheme, `abaumannii_2`). The temporary download is removed after each genome. Progress is written to `logs/processed.txt` (resume log) and `logs/mlst_results_all.tsv`.

```bash
conda activate mlst            # mlst 2.17.6 (patched, see [reproducibility notes](reproducibility_notes.md#troubleshooting)); `datasets` must also be on PATH (the `mlst` env does not contain it)
nohup bash scripts/01_dataset/01_download_filter_ST1.sh > logs_download.out 2>&1 &
# Result: ACB_ST1_pipeline/kept_ST1/<accession>.fna   (1,256 ST1 candidate genomes)
```

> The original run of this script used the scheme name `acinetobacter`; the corrected script uses `abaumannii_2`, the scheme used for the verification in Step 2 and for all reported STs. Every kept genome was re-typed with `abaumannii_2` in Step 2 (1,256/1,256 ST1), so the final collection does not depend on the scheme name used in the first pass. Oxford-scheme STs (`abaumannii`) are **not** reported (BR5 lacks *gyrB* in the assembly, so no Oxford ST can be called).

### Step 2 — Re-verify every kept genome (verified script)

The first MLST log was incomplete (it recorded far fewer lines than genomes), so every genome in `kept_ST1/` was re-typed with full logging:

```bash
conda activate mlst
bash scripts/01_dataset/02_reverify_mlst_pasteur.sh
# Result: ACB_ST1_pipeline/logs/reverify/mlst_reverify_full.tsv
#         1256/1256 genomes confirmed ST1 (Pasteur), 0 failures
```

### Step 3 — Quality filtering: 1,256 candidates → 897 genomes (verified scripts)

Run in this order (CheckM2 v1.1.0, QUAST v5.2.0):

```bash
bash   scripts/01_dataset/03a_run_checkm2_quast.sh      # QUAST (command from quast.log) + CheckM2 (call reconstructed)
python scripts/01_dataset/03b_merge_quast_checkm2.py    # verified original: QUAST + CheckM2 + metadata -> merged_final.tsv
python scripts/01_dataset/03c_filter_quality.py         # verified original: the filter below
```

Genomes are retained when **predicted completeness ≥ 95 %, contamination ≤ 5 %, ≤ 300 contigs (CheckM2 `Total_Contigs`) and ≤ 500 ambiguous bases (N)**, where the N count is derived from the QUAST column `# N's per 100 kbp` × total length (rounded). On the study data this passes **897 genomes (identical to `highquality_accessions.txt`) and fails 359**; no genome has exactly 300 contigs, so `≤ 300` and `< 300` give the same set. The accessions of the 897 genomes used in all analyses are in `highquality_accessions.txt` (Table S2).

`extras/qc_independent_recheck_N_from_fasta.sh` is an independent re-implementation (per-genome QUAST, N counted from the FASTA, `< 300` contigs) that returns the same 897 genomes; it is a cross-check, not the original run. The 897 genomes represent **892 distinct BioSamples**: five BioSamples are represented by two independent assemblies, and these pairs were not collapsed. The 1,256 candidates were de-duplicated beforehand with `scripts_logs/deduplicate_genomes.py` (GCF preferred over GCA, latest version kept).

BR5 itself was assembled from Illumina NextSeq reads (Trimmomatic → SPAdes v3.15; the NCBI assembly record lists the assembler as "SPAdes v. 0.4.8", a submitter typo — no SPAdes release has that number), and species identity was confirmed by ANI against *A. baumannii* ATCC 19606 (`CP045110.1`; ANI 97.84 %, EzBioCloud calculator).

### Step 4 — Annotation with Bakta (verified script)

```bash
bash scripts/02_annotation/01_run_bakta.sh
# Script copied unchanged from the cluster. Bakta v1.11.3, database v6.0 (full, 2025-02-24).
# 7 genomes in parallel x 8 threads, up to 3 attempts per genome, genomes with an existing .json are skipped.
# Used downstream: .gff3 (Panaroo), .faa (AMRFinderPlus, RGI), .tsv (locus tags + coordinates for the
# blaOXA-23 context), .fna (contig sequences).
```

The annotation folder held 898 genomes because one extra FASTA (`GCF_009035845.1`) was still in `highquality_fastas/`; it was excluded afterwards (it is absent from `highquality_accessions.txt` and from Tables S2–S7), and Panaroo and all downstream steps used the 897 genomes in `highquality_accessions.txt`.

### Step 5 — Pangenome and core-genome alignment (reconstructed from outputs; no original script or log survives)

Panaroo v1.7.0 was run on the Bakta GFF3 files of the **897 genomes** (`highquality_accessions.txt`). Verified from the output files:

| Item | Value | Source |
|---|---|---|
| Genes in total | 10,727 | `summary_statistics.txt` |
| Core genes (99–100 % of strains) | 2,246 | `summary_statistics.txt`, recount of `gene_presence_absence.csv` |
| Genes in all 897 genomes | 1,561 | recount of `gene_presence_absence.csv` |
| Soft-core (95–<99 %) / shell / cloud | 613 / 1,534 / 6,334 | `summary_statistics.txt` |
| **Genes in the core-gene alignment (≥ 98 % of strains)** | **2,519** | `aligned_gene_sequences/` (2,519 files) and `alignment_resume_state.json` (`core_threshold: 0.98`, `aligner: mafft`, `alignment: core`) |
| Alignment aligner | MAFFT v7.526 | conda env `panaroo` |
| Alignment length after Panaroo filtering | 2,136,645 positions × 897 sequences | `core_gene_alignment_filtered.aln` (60-column FASTA) |

Note that the **alignment threshold (98 %) is lower than the 99 % used for the "core genes" statistic**: the alignment contains 2,519 genes, whereas 2,246 genes are present in ≥ 99 % of the genomes.

```bash
# RECONSTRUCTED (the original command was not recorded; the alignment step was resumed on 2026-07-01
# with these settings according to alignment_resume_state.json). The clean-mode and thread settings of the
# first pangenome run are NOT recorded and must not be quoted as verified.
panaroo -i bakta_output/*/*.gff3 -o results/pangenome/panaroo_results \
        -a core --aligner mafft --core_threshold 0.98 -t 16
# outputs: gene_presence_absence.csv, summary_statistics.txt, core_gene_alignment.aln,
#          core_gene_alignment_filtered.aln (Panaroo's own filtered alignment, used downstream)

# core_gene_alignment_filtered_clean.aln = the filtered alignment with 5 lowercase IUPAC ambiguity
# characters (w, s, m, m, k) replaced by N; no other position differs (verified with cmp -l).
# All downstream tools (Gubbins, IQ-TREE, ClonalFrameML) used the _clean file.
```

### Step 6 — Recombination and recombination-corrected phylogeny (6a Gubbins reconstructed from logs; 6b IQ-TREE from log, ClonalFrameML settings from log; 6c IQ-TREE verbatim from log, SNP-sites flags reconstructed)

**6a. Gubbins** (v3.4.3) on the core-gene alignment gives recombinant blocks, the proportion of SNPs due to recombination, and **per-branch statistics** (the *mean per-branch r/m* is the average of the per-branch r/m column; the *pooled r/m* is calculated from summed recombinant vs. point-mutation substitutions).

```bash
# RECONSTRUCTED from results/gubbins/gubbins_clean_run.log and results/phylogeny/acb_ST1.log:
# input = core_gene_alignment_filtered_clean.aln, RAxML 8.2.12 GTRGAMMA (20 threads),
# sequences with excess missing data removed by Gubbins' default filter.
run_gubbins.py --prefix acb_ST1 --tree-builder raxml --threads 20 core_gene_alignment_filtered_clean.aln
# acb_ST1.per_branch_statistics.csv, acb_ST1.recombination_predictions.gff,
# acb_ST1.filtered_polymorphic_sites.fasta, acb_ST1.summary_of_snp_distribution.vcf, ...

# mean per-branch r/m (verified on the per-branch table of this study: 2.93304 over 1,793 branches)
awk -F',' 'NR>1{sum+=$2; n++} END{print "Mean r/m:", sum/n, "(n="n" branches)"}' recombination_per_branch.csv
```

**6b. ClonalFrameML** (v1.20) was run independently on the **complete core-gene alignment (2,136,645 positions, invariant sites retained)** with an ML tree inferred from that same alignment (IQ-TREE, GTR+G); uncertainty from **100 simulations** (`-emsim 100`). Because the whole alignment is used, ClonalFrameML block coordinates are **directly positions in the core-gene alignment** — no coordinate conversion is needed.

```bash
# IQ-TREE 3.1.2 -- command taken verbatim from results/clonalframeml_corrected/core_897.log
iqtree3 -s core_gene_alignment_filtered_clean.aln -m GTR+G -T 36 -pre core_897 -fast
# ClonalFrameML 1.20 -- settings read from cfml_897_newtree.log: emsim = 100, show_progress = true, 38 threads,
# 897 sequences x 2,136,645 sites, 58,420 variable sites analysed (thread-option spelling not recorded)
ClonalFrameML core_897.treefile core_gene_alignment_filtered_clean.aln ST1_cfml_897_newtree -emsim 100 -show_progress true
# *.em.txt: global R/theta, 1/delta, nu (posterior mean, variance, a_post, b_post) + per-branch R/theta
# *.emsim.txt: 100 simulated (R/theta, delta, nu) triplets -> 95 % intervals
# Study run (5 Oct 2026): R/theta 0.618892, 1/delta 0.00202071 (delta 494.9 bp), nu 0.0237018
# -> r/m = 0.618892 x 494.9 x 0.0237018 = 7.26
# ClonalFrameML r/m = (R/theta) x delta x nu
```

> **Do not run ClonalFrameML on the SNP-only alignment.** Without invariant sites the R/θ, δ and ν estimates are on the wrong scale (R/θ ≈ 1.3 × 10⁻³, r/m ≈ 4.4 × 10⁻³ were obtained that way and are not comparable with the reported values). The Gubbins r/m and the ClonalFrameML r/m are derived differently and are **not numerically equivalent**; they are reported side by side, each attributed to its own tool.

**6c. Recombination-filtered cgSNP phylogeny** (used for Figures 1–3 and the SNP distances) — `scripts/03_pangenome_phylogeny/04_cgsnp_tree_snpdists.sh`:

```bash
# IQ-TREE 3.1.2 -- command taken verbatim from iqtree_run.logn (run on the Gubbins recombination-filtered polymorphic sites)
iqtree -s acb_ST1.filtered_polymorphic_sites.fasta -m GTR+G -fconst 619984,389653,445493,622792 \
       -bb 1000 -nt 38 -pre iqtree_results/acb_ST1
# -fconst = constant-site counts (A,C,G,T); -bb 1000 = UFBoot2, 1000 replicates; nodes with >= 95 % support are
# considered well supported -> iTOL v6 with metadata.

# Verified from the output files: SNP-sites gave acb_ST1_snpsites.{fasta,phylip,vcf} (897 sequences x 15,282 SNP sites);
# the distance matrix acb_ST1_snp_distance_matrix.tsv (897 x 897) was written by SNP-dists 1.2.0 from acb_ST1_snpsites.fasta
# (three distances recounted from the FASTA equal the matrix).
# NOT RECORDED: the SNP-sites flags and input file; reconstructed in the script as
#   snp-sites -m -v -p -o acb_ST1_snpsites acb_ST1.filtered_polymorphic_sites.fasta
#   snp-dists acb_ST1_snpsites.fasta > acb_ST1_snp_distance_matrix.tsv
```

**6d. (Optional, exploratory) recombination hotspots** — `extras/recombination_hotspots_legacy/`. These scripts translated ClonalFrameML blocks from an **SNP-alignment** run into positions using the Gubbins VCF, binned them into 20-kb windows, mapped windows to Panaroo core genes (`core_alignment_header.embl`) and flagged AMR / virulence genes. Positions are in **Panaroo core-alignment concatenation order, not chromosomal coordinates** (the pseudo-contig is a 2.14 Mb concatenation). They were written for the earlier SNP-alignment run; with the corrected run in 6b the translation step is unnecessary (block coordinates are already alignment positions), and hotspot results are **not part of the reported manuscript results**.

### Step 7 — Resistome (AMRFinderPlus: command from run log; RGI: reconstructed from outputs)

AMRFinderPlus was run on the **Bakta protein FASTA** (BLASTP/EXACTP/ALLELEP/POINTP methods in the output; no contig/start/stop columns), organism `Acinetobacter_baumannii` for point mutations, `--plus` for stress-response / metal / biocide genes, **default curated gene-family thresholds** (fallback 90 % identity / 50 % coverage when no curated cutoff exists). RGI/CARD was used as a cross-check on the same Bakta protein sequences, retaining Perfect and Strict hits only.

```bash
# Command taken from results/amr/amrfinder_run.log (v4.2.7, DB 2026-05-15.1)
amrfinder -p bakta_output/${ACC}/${ACC}.faa --organism Acinetobacter_baumannii --plus \
          --output amrfinder_results/${ACC}.tsv --threads 2
# RECONSTRUCTED from the RGI output tables in results/amr/card_results/ (v6.0.5; local CARD database).
# The tables have blank Contig/Start/Stop and Predicted_DNA columns and ORF ids that are Bakta locus tags,
# i.e. RGI was run in PROTEIN mode on the Bakta .faa files. Only Perfect and Strict hits are present
# (no --include_loose). 886 of 897 genomes returned >= 1 hit; 11 returned an empty result ({}).
# Checked: results/amr/card_results holds exactly one .json + one .txt per genome (1,794 files; no temporary files, as with
# --clean). The alignment tool (-a DIAMOND|BLAST) and thread count were not recorded: the JSON _metadata is empty and
# card_run.log only lists skipped genomes. The options above are valid in `rgi main --help`.
rgi main --input_sequence bakta_output/${ACC}/${ACC}.faa --input_type protein \
         --output_file results/amr/card_results/${ACC} --local --clean
```

The *bla*<sub>OXA-23</sub> local context in BR5 (Figure S1) was annotated in Geneious Prime v2025.0 and plotted in RStudio.

### Step 8 — Genetic context of every *bla*<sub>OXA-23</sub> copy (verified scripts)

```bash
# 8a. Hits -> genomic coordinates -> ±5 kb flanks -> CD-HIT-EST (95 % identity, 90 % coverage)
bash scripts/05_oxa23_context/01_extract_cluster_oxa23_context.sh
#   AMRFinder protein ids are joined to Bakta locus tags to recover contig/start/stop.
#   Result: 848 blaOXA-23 copies in 650 genomes (181 genomes with >= 2 copies, up to 5) -> 58 clusters

# 8b. Cluster representatives vs the AB0057 chromosome (CP001182.2; confirmed Tn2006-in-AbaR4)
bash scripts/05_oxa23_context/02_cluster_reps_vs_AB0057.sh
#   Result: the three largest contexts (662/848 copies, 78.1 %) align over their full length
#           (97.6-99.7 % identity, 100 % query coverage)

# 8c. Build the curated reference database (v4)
bash scripts/05_oxa23_context/03_build_abar_reference_db.sh

# 8d. Classify ALL 848 copies (not only cluster representatives)
bash scripts/05_oxa23_context/04_classify_oxa23_context.sh
```

Key design decisions (each fixed a real bias found during the analysis):

- **Isolated regions instead of whole replicons.** AB0057 (`NC_011586.2`), Tn2008B (`CP085788.1`) and Tn2009 (`NZ_KM922672.1`) are whole chromosomes/plasmids; used whole, background similarity inflated the matches. The *bla*<sub>OXA-23</sub> block is located with the standalone Tn2006 probe (`JN129846.1`) and only that block ± 5 kb is retained (`CP085788.1:521382-534181`, `NZ_KM922672.1:64129-77323`, `NC_011586.2:577113-594181`). With the bias removed, **Tn2008B has no confirmed hits**, while Tn2009 is confirmed.
- **Normalised coverage** = alignment length / min(query length, reference length). Raw query coverage penalised complete matches to short references (e.g. Tn2008, 4,593 bp).
- **Thresholds:** < 99 % identity or < 90 % normalised coverage → *unassigned/atypical*. AbaR4, Tn2006 and the isolated AB0057 region overlap and are collapsed into **Tn2006/AbaR4-type**.
- Copy-level result on the 848 copies: Tn2006/AbaR4-type 538, unassigned/atypical 242, Tn2009 37, Tn2008 26, Tn2007 3, AbaR25 2, Tn2008B 0.
- **AbaR3/Tn6019, AbaR4, AbaR4a and AbaR25 calls** used in the association tests come from a **separate genome-wide screen** (Step 8e, Table S5), not from this flanking-region database.

#### Step 8e — Genome-wide AbaR island screen (part of Step 8; new script; run before Step 9)

Each island reference (AbaR3 = the resistance island of AB0057, `NC_011586.2:235738-326347`; AbaR4, AbaR4a, AbaR25 from the NCBI records listed in Step 8c) is aligned with BLASTN (BLAST+ 2.5.0+) against every genome. An island is called **present when ≥ 90 % identity alignments cover ≥ 80 % of its length** (union of aligned intervals). Calls for the AbaR4 variants overlap because they share most of their sequence.

```bash
bash scripts/05_oxa23_context/05_abar_genome_screen.sh   # -> results/abar_island_screen/abar_presence.tsv, abar_fraction_covered.tsv
```

**Validation of the screen (run on the study data):** the AB0057 chromosome (`GCF_000021245.2`), the source of the AbaR3 reference, covers 100 % of it (positive control). AbaR3 calls agree with the island's resistance cargo: *tet*(A) is present in 202/208 AbaR3-positive genomes versus 4/689 negatives, *catA1* in 188/208 versus 5/689; *tet*(A) is absent from all 643 genomes with < 0.7 AbaR3 coverage and present in 165/166 genomes with ≥ 0.9. Genomes between 0.5 and 0.7 coverage carry only the shared AbaR-type backbone (e.g. AB5075, which matches the two ends of the reference). Cross-check against the independent context classification (Table S6): AbaR4, AbaR4a and AbaR25 calls are almost absent from genomes without *bla*<sub>OXA-23</sub> (1/247, 1/247, 0/247) and from genomes with Tn2008 contexts (0/26), and are found in 50-58 % of genomes with a Tn2006/AbaR4-type context (the others carry Tn2006 outside a complete island or have fragmented assemblies). AbaR3 and AbaR4 are mutually exclusive (1 genome has both), as expected for alternative islands at the same site (*comM*). AbaR4, AbaR4a and AbaR25 overlap almost completely (237/245 AbaR4 calls are also AbaR25) and are therefore analysed as one **AbaR4-type** group. The association of AbaR3 with *bla*<sub>OXA-23</sub> keeps the same direction for cut-offs of 0.7, 0.8 and 0.9 (OR 0.14, 0.11, 0.13); AbaR4, AbaR4a and AbaR25 stay strongly positive at every cut-off tested. AbaR4b (= Tn2006, which contains *bla*<sub>OXA-23</sub>) and AbaR4c (a 4 kb *tniA* fragment found in 792/897 genomes) are not specific island calls.

### Step 9 — Association tests (verified script)

Fisher's exact test (SciPy v1.18.1; chosen over χ² because several cells are sparse or zero) of *bla*<sub>OXA-23</sub> carriage against ISAba1 and AbaR3/Tn6019, AbaR4, AbaR4a and AbaR25 (and the combined AbaR4-type group) across the 897 genomes:

```bash
python scripts/07_statistics/01_fisher_exact_oxa23_associations.py Supplementary_tables.xlsx results/abar_island_screen/abar_presence.tsv
```

### Step 10 — Sequence typing, capsule/OC loci, virulome and mobilome

**MLST** — Steps 1–2 (Pasteur scheme). **KL/OCL (Kaptive 2.0.6; prefixes match the study outputs, database paths and options are a template):**

```bash
kaptive.py -a ACB_ST1_pipeline/kept_ST1/*.fna -k <Acinetobacter_baumannii_k_locus_primary_reference.gbk> -o results/kaptive/acb_ST1_K_locus
kaptive.py -a ACB_ST1_pipeline/kept_ST1/*.fna -k <Acinetobacter_baumannii_OC_locus_primary_reference.gbk> -o results/kaptive/acb_ST1_OC_locus
# confidence tiers (perfect, very high, high, good, low, none) are reported; all tiers retained
```

**Virulome (verified script):** `scripts/06_typing_virulome_mobilome/01_build_acinetobacter_vfdb_and_screen.sh` rebuilds the **966-sequence *Acinetobacter* subset of VFDB Set B** and screens all genomes with ABRicate v1.0.1 at its **default thresholds (≥ 80 % identity, ≥ 80 % coverage)**. The per-genome tables are merged into the virulome matrix (Table S6; 897 genomes × 109 gene columns, checked cell by cell against the merged source table with 0 mismatches). Gene-family columns for Figure 3 combine genes of one system (e.g. *pil*, *csu*, *bfm*) and mark a family present only if **all** its genes are present.

**Mobilome (versions verified; ISEScan/IntegronFinder/IslandPath settings from logs, the rest templates):**

```bash
mob_recon -i ${ACC}.fna -o results/mobsuite/${ACC} -n 8                       # MOB-suite v3.1.9
isescan.py --seqfile ${ACC}.fna --output results/isescan/${ACC} --nthread 8    # ISEScan v1.7.3
integron_finder highquality_fastas/${ACC}.fna --outdir results/integrons --cpu 8           # IntegronFinder v2.0.6
# S5 rules (recovered by comparing Table S5 with the raw outputs; 897/897 genomes agree for integrons, IS and islands):
#   integrons = type reported anywhere in the genome; IS families = ISEScan complete copies (type c);
#   named IS = any ISfinder BLASTN hit (>=80% identity); N Plasmid Contigs = sum of num_contigs in mobtyper_results.txt;
#   replicons / MOB from mobtyper_results.txt (rep_cluster_ shortened to rep_); MPF = mate-pair-formation hit with qcovs >= 90 in
#   biomarkers.blast.txt (5 of 3,588 cells differ); island length = sum(abs(end - start)).
# ISfinder BLASTN search and PlasmidFinder: see scripts/06_typing_virulome_mobilome/03_isfinder_blastn_plasmidfinder.sh
# Replicon, relaxase and MPF types in Table S5 come from MOB-suite (Step 10). PlasmidFinder 2.1.6 was also run (Enterobacteriaceae and
# Gram-positive databases only; hits in 8 of 897 genomes, all non-Acinetobacter replicons) and is not used for the table.
# ResFinder (acquired genes, command read from its JSON): scripts/04_resistome/02_resfinder_batch.sh
# IslandPath-DIMOB v1.0.6: https://github.com/brinkmanlab/islandpath
```

### Step 11 — *Galleria mellonella* virulence and haemolymph proteome figure (wet lab + R scripts)

Groups of 10 larvae (230–280 mg) were injected with 10 µL of a 1 × 10⁸ CFU/mL suspension (PBS control; *A. baumannii* ATCC 19606 as the non-virulent comparator), incubated at 37 °C, three biological replicates (n = 30 per strain). Colonisation: haemolymph CFUs at 3, 6 and 9 h (three pooled larvae per aliquot, MacConkey agar; growth curves plotted in Origin). Haemolymph collected at 15 h post-infection was used for the proteomic comparison described in the manuscript (Methods 2.6–2.7).

```r
# Survival (R: survival, survminer, coin): Kaplan-Meier curves, global log-rank test, pairwise
# Fleming-Harrington (rho = 0, gamma = 1) tests with FDR correction, pairwise-p-value heatmap and workbook (Supplementary Table S8)
source("scripts/08_galleria_proteome/01_galleria_survival.R")
# Proteome heatmap (R: pheatmap): row z-scores of the intensities of the differentially abundant proteins
source("scripts/08_galleria_proteome/02_proteome_heatmap.R")
```

The scripts read `gm_a.baumannii_BR5.xlsx` (live larvae per time point and replicate) and `proteome/reguladas_pvalue.txt` (differentially abundant proteins); set the folder paths at the top of each script.
