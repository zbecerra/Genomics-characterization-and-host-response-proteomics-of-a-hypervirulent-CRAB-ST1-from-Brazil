# Requirements, databases and repository layout


## Software

| Tool | Version | Purpose | Conda env |
|---|---|---|---|
| Trimmomatic / SPAdes | 0.41 / 3.15 | BR5 read trimming and *de novo* assembly | not applicable (BR5 assembly) |
| NCBI Datasets CLI | 18.30.1 (conda base) | genome retrieval | (base) for v18.30.1; the `download` env holds v18.7.0 — confirm which ran Step 1* |
| mlst (Seemann) | 2.17.6 (BLAST+ patch, see [reproducibility notes](reproducibility_notes.md#troubleshooting)) | MLST, Pasteur scheme `abaumannii_2` | `mlst` |
| CheckM2 / QUAST | 1.1.0 / 5.2.0 | genome quality | QUAST: `quast_5.2.0`; CheckM2: `checkm2` |
| Bakta | 1.11.3 (DB v6.0) | annotation | `bakta` (system conda installation; export by path) |
| Panaroo | 1.7.0 | pangenome, core-gene alignment | `panaroo` |
| Gubbins | 3.4.3 (RAxML 8.2.12, GTRGAMMA) | recombination detection | `gubbins` (also holds RAxML, SNP-sites) |
| ClonalFrameML | 1.20 | R/θ, δ, ν | `clonalframe` |
| IQ-TREE | 3.1.2 | ML trees (GTR+G), UFBoot2 | `iqtree` (3.1.2; also in `alignment`) |
| SNP-sites / SNP-dists | 2.5.1 / 1.2.0 (conda env `snpsites`) | cgSNPs / distance matrix | `snpsites` |
| AMRFinderPlus | 4.2.7 (DB 2026-05-15.1) | AMR genes, point mutations | `amrfinder` |
| RGI / CARD | 6.0.5 / 4.0.1 | AMR cross-check (Perfect + Strict only) | `card_new` (RGI 6.0.5; the older `rgi` env holds RGI 4.0.3 and was not used); CARD database kept separately |
| ResFinder / database | 4.7.2 / 2.6.0 | acquired resistance genes (≥ 90 % identity, ≥ 60 % coverage) | `resfinder` |
| SeqKit / CD-HIT-EST | 2.13.0 / 4.8.1 | extraction / clustering of flanking regions | `alignment` |
| BLAST+ | 2.16.0 (flanking-region / context comparison, v2-v4 island databases are in BLAST v5 format) and 2.5.0+ (genome-wide AbaR island screen, ISfinder search) | similarity searches | `alignment` (2.16.0), `mge` / base (2.5.0)* |
| Entrez Direct (`efetch`) | any recent | reference downloads | any |
| ABRicate | 1.0.1 | virulence screening | `abricate` |
| Kaptive | 2.0.6 | KL / OCL typing | `kaptive` |
| MOB-suite / ISEScan / IntegronFinder / IslandPath-DIMOB | 3.1.9 / 1.7.3 / 2.0.6 / 1.0.6 | mobilome | `mobsuite` / `mge` / `integron`* / `islandpath` |
| PlasmidFinder | 2.1.6 (BLAST method; Enterobacteriaceae and Gram-positive databases) | screened, but not used for Table S5 (see Step 10) | `plasmidfinder` |
| ISfinder / BLAST+ | database downloaded 2026-07-28 (5,970 sequences) / BLAST+ 2.5.0+ (conda env `mge`) | IS search by BLASTN | `mge` |
| Python / SciPy | 3.13 / 1.18.1 | Fisher's exact test | user site (`pip install --user`) |
| iTOL | v6 (web) | tree visualisation | web |
| RStudio / R and packages (survival, survminer, coin, pheatmap, ggplot2) | RStudio 2025.09.2+418 (Build 418) / R 4.3.3; package versions not recorded | *G. mellonella* survival statistics, proteome heatmap | RStudio / R (not conda) |

## Databases

- **VFDB Set B nucleotide** (`VFDB_setB_nt.fas`, http://www.mgc.ac.cn/VFs/) → filtered to *Acinetobacter* entries (966 sequences).
- **Reference island/transposon sequences** (downloaded from NCBI by `scripts/05_oxa23_context/03_build_abar_reference_db.sh`): AbaR25 `JX481978.1`, TnAbaR4a `JN129845.1`, AbaR4c fragment `JN129847.1`, AbaR4 `HQ700358.2`, Tn2006 `JN129846.1`, Tn2007 `EF059914.1`, Tn2008 `KP780408.1`, Tn2008B `CP085788.1`, Tn2009 `NZ_KM922672.1`, AB0057 `NC_011586.2` (= `CP001182.2`).
- **Species-confirmation reference:** *A. baumannii* ATCC 19606, `CP045110.1`.
- AMRFinderPlus DB, CARD, Bakta DB, Kaptive *A. baumannii* K/OC locus databases.

\* Environment names marked with an asterisk were inferred from the environments found on the cluster (`conda list`), not from run logs. Environment files for the tools are collected in [`envs/`](../envs/README.md) (create them with `conda env create -f envs/<name>.yml`).

## Repository structure

```
.
├── README.md
├── docs/                        pipeline_steps.md, requirements.md, reproducibility_notes.md
├── envs/                        conda environment files (export_envs.sh generates them)
├── scripts/
│   ├── 00_setup/                patch_mlst_blast_check.sh (mlst + BLAST+ >= 2.10)
│   ├── 01_dataset/              01 download + ST1 filter, 02 re-verify MLST, 03a-c QC (CheckM2/QUAST) -> 897 genomes
│   ├── 02_annotation/           01_run_bakta.sh
│   ├── 03_pangenome_phylogeny/  01 Panaroo, 02 Gubbins, 03 ClonalFrameML, 04 cgSNP tree + SNP distances
│   ├── 04_resistome/            01_amrfinder_rgi_batch.sh
│   ├── 05_oxa23_context/        01-04: blaOXA-23 flanks -> clusters -> AbaR/Tn database -> classification; 05: genome-wide AbaR island screen
│   ├── 06_typing_virulome_mobilome/   01 virulome (VFDB, ABRicate), 02 Kaptive / MOB-suite / ISEScan / IntegronFinder
│   ├── 07_statistics/           01_fisher_exact_oxa23_associations.py
│   └── 08_galleria_proteome/    01 survival analysis (R), 02 proteome heatmap (R)
└── extras/                      not needed to reproduce the paper (QC cross-check, legacy recombination hotspots)
```

Script status is given at the top of every script: **verified** (recovered unchanged), **from log**, **reconstructed** or **template** (see "How to read this repository").

Expected working layout on the cluster (paths used by the scripts):

```
projects/assacinetobacter/
├── ACB_ST1_pipeline/{accessions,kept_ST1,logs}/    # Step 1-2 outputs
├── bakta_output/<accession>/<accession>.{fna,faa,gff3,tsv,...}
├── results/amr/amrfinder_results/<accession>.tsv
├── results/{pangenome/panaroo_results, gubbins, kaptive}/
├── databases/abar/                                 # Step 8
├── oxa23_context/                                  # Step 8
└── vfdb_rebuilt_results/<accession>_acb_vfdb_rebuilt.tsv
```
