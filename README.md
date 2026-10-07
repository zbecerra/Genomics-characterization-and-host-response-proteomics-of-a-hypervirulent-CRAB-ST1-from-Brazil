# Genomics and host-response proteomics of a hypervirulent carbapenem-resistant *Acinetobacter baumannii* ST1 isolate BR5 from Brazil

Reproducible pipeline, scripts and reference databases accompanying the manuscript:

> Becerra J., Merighi D.G.S., Vásquez-Ponce F., Martins-Gonçalves T., Pompeu G.C., Pascale C.B.A., Palmisano G., Lincopan N., Guzzo C.R. *Genomics characterization and host response proteomics of a hypervirulent carbapenem-resistant Acinetobacter baumannii ST1 isolate from Brazil.* [Journal, year — update on acceptance]

**Repository author and sole maintainer:** J. Becerra (zbecerra@usp.br), Department of Microbiology, Institute of Biomedical Sciences, University of São Paulo (USP), Brazil. The manuscript authors listed above are credited for the paper only; all code in this repository was written and curated by the repository author.

## Before release: open items

- [ ] Zenodo DOI for the archived release (see *Data availability*).
- [ ] Journal reference and year in the citation above (update on acceptance).
- [ ] Run `bash envs/export_envs.sh` on the cluster and commit the generated `envs/*.yml` files.
- [ ] Confirm the few conda environments still marked with an asterisk in [`docs/requirements.md`](docs/requirements.md) (Datasets CLI, RGI/CARD).
- [ ] Add `sessionInfo()` (R package versions) for the *G. mellonella* scripts.


## Overview

BR5 (GenBank assembly `GCA_044133445.1`; Brazil, 2021; sputum (tracheal secretion) from a fatal pneumonia case) is an ST1<sup>Pas</sup>, IC1, KL17:OCL1 carbapenem-resistant *A. baumannii* carrying *bla*<sub>OXA-23</sub>. This repository documents how BR5 was placed in the context of **896 additional ST1 genomes (897 in total)** downloaded from NCBI, and how the following were derived:

1. the ST1 genome collection (download, de-duplication, MLST filtering, quality filtering);
2. phylogenomics and recombination (Panaroo, Gubbins, ClonalFrameML, IQ-TREE);
3. the resistome and the genetic context of every *bla*<sub>OXA-23</sub> copy (AMRFinderPlus, RGI/CARD, Bakta, CD-HIT-EST, BLASTN against a curated transposon/island database);
4. KL/OCL typing, virulome (ABRicate + Acinetobacter-specific VFDB subset) and mobilome;
5. *Galleria mellonella* virulence assays and haemolymph proteomics (wet-lab; described in the manuscript).

## Pipeline overview

```
NCBI Datasets (taxid 470, dehydrated manifest)
 ├─ de-duplicate: GCF over GCA, latest version only
 ├─ download one genome at a time -> MLST (Pasteur, abaumannii_2) -> keep ST1      [scripts/01_dataset]
 ├─ re-verify every kept genome with MLST (1,256/1,256 = ST1)                      [scripts/01_dataset]
 └─ quality filtering (CheckM2, QUAST)  ->  897 genomes
        │
        ├─ Bakta v1.11.3 annotation
        │     └─ Panaroo v1.7.0 -> core-gene alignment (2,519 genes at >=98 %, 2,136,645 positions)
        │            ├─ Gubbins v3.4.3 -> recombination, per-branch r/m, filtered cgSNPs
        │            │       └─ SNP-sites -> SNP-dists -> IQ-TREE (GTR+G, -fconst, UFBoot2) -> iTOL
        │            └─ IQ-TREE (full alignment, GTR+G) -> ClonalFrameML (-emsim 100): R/θ, δ, ν
        ├─ AMRFinderPlus v4.2.7 (+ RGI v6.0.5 / CARD v4.0.1)
        │     └─ blaOXA-23 copies (n = 848) -> Bakta coordinates -> ±5 kb -> CD-HIT-EST
        │            -> BLASTN vs curated AbaR/Tn database (v4)  -> Fisher's exact tests   [scripts/05]
        ├─ Kaptive v2.0.6 (KL / OCL)
        ├─ ABRicate v1.0.1 + 966-sequence Acinetobacter VFDB subset (virulome)           [scripts/06]
        └─ MOB-suite, ISEScan, IntegronFinder, IslandPath-DIMOB (mobilome)
```

## Quickstart

1. Install the tools listed in [`docs/requirements.md`](docs/requirements.md) (one conda environment per tool; files in [`envs/`](envs/README.md)).
2. If you use `mlst` 2.17.6 with BLAST+ ≥ 2.10, apply the patch once: `bash scripts/00_setup/patch_mlst_blast_check.sh`.
3. Work from one project directory and run the scripts in numeric order (`scripts/01_dataset` … `scripts/08_galleria_proteome`). Each script starts with a header that states its inputs, outputs and status (see below).
4. The expected directory layout and every step are described in [`docs/pipeline_steps.md`](docs/pipeline_steps.md).


## Confidence per step

How each step was recovered. **Verified** = script recovered unchanged; **From log** = command read from a run log; **Reconstructed** = assembled from outputs and logs (flags not visible in a log are tool defaults); **Template** = parameters follow the Methods but the original command was not recorded. Where a template step feeds a supplementary table, its output was compared with the table (agreement is stated).

| Step | Content | Status |
|---|---|---|
| 1–2 | Download, de-duplication, ST1 filter, MLST re-verification | Verified |
| 3 | Quality filtering (CheckM2, QUAST) → 897 genomes | Verified scripts |
| 4 | Bakta annotation | Verified |
| 5 | Panaroo pangenome, core-gene alignment | Reconstructed (original command not kept) |
| 6a | Gubbins | Reconstructed from logs |
| 6b | IQ-TREE, ClonalFrameML | From log (ClonalFrameML thread option spelling not recorded) |
| 6c | cgSNP tree, SNP-dists | IQ-TREE from log; SNP-sites flags reconstructed; matrix verified |
| 7 | AMRFinderPlus; RGI; ResFinder | AMRFinderPlus and ResFinder from log; RGI reconstructed |
| 8 (8a–8d) | *bla*<sub>OXA-23</sub> flanks, CD-HIT-EST, reference database, classification | Verified scripts |
| 8e | Genome-wide AbaR island screen | New script, validated (positive control, gene content, threshold sensitivity) |
| 9 | Fisher's exact tests | Verified script |
| 10 | Virulome (ABRicate/VFDB) | Verified script |
| 10 | Kaptive, MOB-suite, ISEScan, IntegronFinder, IslandPath-DIMOB, ISfinder BLASTN | Template; Table S5 rules recovered and checked against raw outputs (897/897 genomes for integrons, IS and islands; replicons and relaxases agree; 5 of 3,588 MPF cells differ) |
| 11 | *G. mellonella* survival and proteome scripts (R) | Verified (original R scripts; paths made relative) |


## Documentation

| Page | Contents |
|---|---|
| [`docs/pipeline_steps.md`](docs/pipeline_steps.md) | Step-by-step commands (Steps 1–11) and how to read command labels |
| [`docs/requirements.md`](docs/requirements.md) | Tool versions with conda environment, databases, repository layout |
| [`docs/reproducibility_notes.md`](docs/reproducibility_notes.md) | Output files, troubleshooting (including the `mlst` patch), key references |


## Key results summary

| Analysis | Result |
|---|---|
| Genome collection | 897 ST1<sup>Pas</sup> genomes, 70 countries, 1982–2026 (1,256 ST1 candidates before quality filtering) |
| Capsule / outer-core loci | 19 KL types (KL17 n = 336, 37.5 %; KL1 n = 249, 27.8 %); 7 OCL types (OCL1 n = 653, 72.8 %) |
| Recombination (Gubbins) | 82.87 % of SNPs inside recombination blocks (83,487 of 100,744); pooled r/m 4.84 (83,487 / 17,257); mean per-branch r/m 2.93 over 1,793 branches; 4,347 blocks on 353 branches |
| Recombination (ClonalFrameML) | R/θ = 0.619 (95 % CI 0.599–0.637, 100 simulations); mean recombinant tract length δ = 495 bp (inferred imports average 666 bp); ν = 0.0237; r/m = 7.26 (95 % CI 7.03–7.51); 4,710 imports on 385 branches |
| BR5 sublineage | clusters with 4 Brazilian OXA-23-positive isolates (2010–2020) at 18–47 SNPs; ≥ 101 SNPs from all non-Brazilian genomes |
| *bla*<sub>OXA-23</sub> | 650/897 genomes (72.5 %); 848 copies; 58 contexts; Tn2006/AbaR4-type dominant, Tn2008 / Tn2009 / Tn2007 / AbaR25 minor |
| Associations (Fisher, 897 genomes) | ISAba1 OR 38.2 (p = 2.3 × 10⁻⁴⁹; 639/650 OXA-23⁺ vs 149/247 OXA-23⁻ genomes); AbaR3/Tn6019 in 208 genomes (74/650 vs 134/247), **inversely** associated with *bla*<sub>OXA-23</sub> (OR 0.11, p = 1.0 × 10⁻³⁸); AbaR4-type island (AbaR4, AbaR4a or AbaR25) in 284 genomes (283/650 vs 1/247; OR 189.7, p = 5.8 × 10⁻⁴⁸), partly by construction because these islands carry Tn2006 |
| Virulome | 966-sequence *Acinetobacter* VFDB subset; BR5: 111 predicted virulence genes |
| *G. mellonella* | BR5: 100 % mortality by 20 h (ATCC 19606: no mortality; log-rank p < 0.0001); ~10⁹ CFU/mL haemolymph by 9 h |
| Haemolymph proteome | 28 differentially abundant proteins (12 higher with ATCC 19606, 16 higher with BR5) |

## Supplementary tables

Supplementary tables (numbering as in the supplementary workbook):

| Table | Contents | Source in this repository |
|---|---|---|
| S1 | AST profile; resistome and virulome of BR5 | BR5 analysis (separate table) |
| S2 | Genomic and epidemiological information of the 897 genomes (accession, strain, date, country, source, QC metrics) | Steps 1–3 |
| S3 | Pairwise SNP-distance matrix (897 × 897; SNP-dists 1.2.0) | Step 6c, `acb_ST1_snp_distance_matrix.tsv` |
| S4 | AMR matrix: AMRFinderPlus, RGI/CARD (Perfect + Strict) and ResFinder | Step 7 |
| S5 | Mobile genetic elements matrix (MOB-suite, ISEScan + ISfinder BLASTN, IntegronFinder, IslandPath-DIMOB; AbaR islands from the genome-wide screen) | Steps 8e and 10 |
| S6 | Genetic-context classification of the 848 *bla*<sub>OXA-23</sub> copies (650 genomes) | Step 8d |
| S7 | Virulome matrix (ABRicate/VFDB) with K and OC loci (Kaptive) | Step 10 |
| S8 | *G. mellonella* survival: pairwise Fleming-Harrington p-values (FDR-adjusted) | Step 11, `01_galleria_survival.R` |
| Recombination and phylogeny workbook | Gubbins / ClonalFrameML summary statistics and per-branch tables | Step 6 |

## Data availability, citation and license

- BR5 reads and assembly: NCBI BioProject [PRJNA1148164](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1148164), BioSample [SAMN43188146](https://www.ncbi.nlm.nih.gov/biosample/SAMN43188146), assembly `GCA_044133445.1` (ASM4413344v1; released 2024-10-30; submitter University of Sao Paulo). Illumina NextSeq, 292x coverage, assembled with SPAdes v3.15 (the NCBI record lists "v. 0.4.8", a submitter typo); contig-level assembly, 42 contigs, 3,755,206 bp, N50 152,642 bp. Isolated in 2021 from sputum (tracheal secretion) in Presidente Prudente, Brazil.
- Comparator genomes: public NCBI accessions (Table S2).
- A versioned release of this repository will be archived on Zenodo: **[add DOI]**.
- License: MIT (see [`LICENSE`](LICENSE)); the scripts may be reused and modified with attribution. Third-party databases and tools keep their own licences.

Author and maintainer: J. Becerra — zbecerra@usp.br (issues and questions via GitHub Issues).
