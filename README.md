# sRNA Metatranscriptomics Pipeline for 20.440

A computational pipeline for the discovery, structural characterization, and functional analysis of small non-coding RNAs (sRNAs) in environmental microbiomes. The pipeline covers raw read download through differential expression, Random Forest classification, sRNA–mRNA interaction prediction, and figure generation

Two datasets are supported for comparison initially but only the AD dataset is used for this class (and CB is dropped from the later scripts but it should be easy to correct by pointing to the correct dataset directory) 

| ID | Environment | Timepoints |
|----|-------------|------------|
| **AD** | Atacama Desert halite | 2016, 2017 |
| **CB** | Chesapeake/Delaware Bay | Fall, Summer |

---

## Table of Contents

- [Overview](#overview)
- [Pipeline Steps](#pipeline-steps)
- [Requirements](#requirements)
- [Directory Structure](#directory-structure)
- [Metadata Format](#metadata-format)
- [Running the Pipeline](#running-the-pipeline)
- [Key Outputs](#key-outputs)
- [Script Reference](#script-reference)

---

## Overview

```
Raw reads (SRA)
    │
    ▼
01_download.sh          ← prefetch + fasterq-dump
    │
    ▼
02_kneaddata_assembly.sh ← QC (KneadData) → co-assembly (MEGAHIT) → annotation (Prodigal) → index (Bowtie2)
    │
    ├──▶ 03_humann.sh   ← functional profiling (HUMAnN + MetaPhlAn)
    │
    ├──▶ 04_snapt.sh    ← sRNA detection (SnapT)  [per timegroup]
    │         │
    │    05b_bam.sh     ← RNA-seq BAM alignment (Bowtie2 + SAMtools)
    │         │
    │    05_rnafold.sh  ← secondary structure (RNAfold)  [per timegroup]
    │         │
    │    06_intarna.sh  ← sRNA–mRNA interaction prediction (IntaRNA)  [per timegroup]
    │         │
    │    09_gene_to_uniref.sh  ← protein → UniRef90 (DIAMOND blastp)
    │
    ▼
07_feature_matrix.sh    ← build Random Forest feature matrix
    │
    ▼
08_run_analysis.sh      ← differential expression, RF classifier, target correlations
    │
08b_per_replicate_tpm.sh ← per-replicate TPM matrices (Pearson QC)
    │
    ▼
10_table_generation.sh  ← master summary tables
    │
    ▼
11_plot_figures.py      ← publication figures
```

---

## Pipeline Steps

### 01 — Download (`01_download.sh`)
Downloads raw FASTQ files from NCBI SRA using `prefetch` + `fasterq-dump`. Reads are renamed from SRR accession to sample name and organized by dataset. Accepts dataset (`AD`/`CB`/`ALL`) and data type (`Metagenome`/`Metatranscriptome`/`ALL`) filters.

### 02 — QC, Assembly & Annotation (`02_kneaddata_assembly.sh`)
1. **KneadData** — removes host contamination (human DNA/RNA and rRNA) from all samples.
2. **MEGAHIT** — co-assembles all DNA reads per dataset (min contig length 500 bp).
3. **Prodigal** — predicts protein-coding genes on the co-assembly (metagenome mode).
4. **Bowtie2-build** — indexes the co-assembly for RNA-seq alignment.

### 03 — Functional Profiling (`03_humann.sh`)
Runs **HUMAnN** (with MetaPhlAn) on both DNA and RNA samples to produce pathway abundance, gene family, and taxonomic profile tables. Tables are joined across samples at the end.

### 04 — sRNA Detection (`04_snapt.sh`)
Runs **SnapT** per timegroup on the combined metatranscriptomic reads, producing antisense (asRNA) and intergenic (itsRNA) small ncRNA annotations in GFF format. FASTA headers in the co-assembly are sanitized before the run.

### 05b — BAM Alignment (`05b_bam.sh`)
Aligns cleaned RNA-seq reads from each replicate back to the co-assembly using **Bowtie2**, then sorts and indexes with **SAMtools**. BAM files are required for per-replicate quantification.

### 05 — Secondary Structure (`05_rnafold.sh`)
Runs **RNAfold** (`--partfunc`) on every sRNA FASTA in parallel. Parses MFE, ensemble energy, GC content, and structural entropy into a summary TSV.

### 06 — sRNA–mRNA Interaction Prediction (`06_intarna.sh`)
1. Filters CDS regions to expressed contigs and caps target regions at 150 nt.
2. Runs **IntaRNA** (all-vs-all, no-seed mode) between sRNAs and mRNA 5′ regions.
3. Filters interactions by p-value (default `< 0.05`) and writes per-sRNA summary statistics.

### 07 — Feature Matrix (`07_feature_matrix.sh` / `07_build_feature_matrix.py`)
Assembles the Random Forest feature matrix from: SnapT FPKM log₂FC, RNAfold structural features, HUMAnN pathway abundances, MetaPhlAn taxonomic profiles, and (optionally) IntaRNA interaction features. Outputs both a full-sample and a species-level matrix.

### 08 — sRNA Analysis (`08_run_analysis.sh` / `08_srna_analysis.py`)
Core analysis script (five sections):
1. Differential expression between timepoints (permutation test).
2. Structural comparison: asRNA vs itsRNA.
3. Random Forest classifier — predicts sRNA up/downregulation from features.
4. sRNA–target expression correlation (log₂FC approximation).
5. Summary figures.

### 08b — Per-Replicate TPM (`08b_per_replicate_tpm.sh`)
Quantifies each replicate BAM individually with **StringTie** (`-e` mode against the SnapT GFF), then assembles a gene × sample TPM matrix. Used for Pearson correlation QC.

### 09 — Gene → UniRef90 (`09_gene_to_uniref.sh`)
Maps Prodigal protein predictions to UniRef90 clusters via **DIAMOND blastp** (e-value ≤ 1×10⁻⁵, top hit per query). Required for pathway-level target annotation in the table generation step.

### 10 — Table Generation (`10_table_generation.sh` / `10_table_generation.py`)
Integrates all upstream results into five master tables (see [Key Outputs](#key-outputs)).

### 11 — Figure Generation (`11_plot_figures.py`)
Produces all publication figures: alpha-diversity panels, pathway heatmaps, sRNA presence/absence volcano plots, Random Forest importance charts, and network diagrams. Run directly after step 10.

---

## Requirements

### Conda environments

| Environment | Key tools |
|-------------|-----------|
| `snapt_env` | SnapT, StringTie, Bowtie2, SAMtools, RNAfold, IntaRNA, MEGAHIT, Prodigal, Python ≥ 3.9 |
| `humann_env` | HUMAnN 4, MetaPhlAn 4, DIAMOND |
| `metagenome_env` | KneadData |

### Python packages (snapt_env)
```
numpy pandas scipy matplotlib seaborn scikit-learn statsmodels networkx adjustText
```

### Databases
| Database | Used by | Path (default) |
|----------|---------|---------------|
| Human genome (DNA) | KneadData | `~/human_dna_db` |
| Human transcriptome (RNA) | KneadData | `~/human_rna_db` |
| rRNA sequences | KneadData | `~/rrna_db` |
| MetaPhlAn SGB database | HUMAnN / MetaPhlAn | `mpa_vOct22_CHOCOPhlAnSGB_202403` | This is currently the database specified by biobakery, may change with updates from them
| UniRef90 (DIAMOND) | DIAMOND blastp | `~/humann_db/uniref/humann4_protein_database_filtered_v2019_06.dmnd` |
| NCBI NR (DIAMOND) | SnapT | `/path/to/nr.dmnd` |
| Rfam (covariance models) | SnapT | `/path/to/Rfam.cm` |

---

## Directory Structure

```
20_440/
├── Metadata.csv
├── AD/
│   ├── <sample>_knead_out/
│   ├── coassembly/
│   │   ├── final.contigs.fa
│   │   └── bowtie2_index.*
│   ├── annotation/
│   │   ├── prodigal.gff
│   │   └── prodigal.faa
│   ├── bam_files/
│   ├── snapt_2016/
│   ├── snapt_2017/
│   ├── rnafold_2016/
│   ├── rnafold_2017/
│   ├── intarna_2016/
│   ├── intarna_2017/
│   ├── per_replicate_tpm_2016/
│   ├── per_replicate_tpm_2017/
│   ├── Total_DNA_pathabundance_humann_table.tsv
│   ├── Total_RNA_pathabundance_humann_table.tsv
│   ├── Total_DNA_merged_abundance_table.tsv
│   ├── feature_matrix.tsv
│   └── analysis_output_2016_2017/
└── CB/
    └── (same structure, timegroups: Fall / Sum)
```

---

## Metadata Format

`Metadata.csv` is a **tab-separated** file with one row per sequencing run:

| Column | Description |
|--------|-------------|
| `BioProject` | NCBI BioProject / BioSample accession |
| `SRR` | SRA run accession |
| `Type` | `Metagenome` or `Metatranscriptome` |
| `Name` | Sample name (used as directory name) |
| `Title` | Free-text description |
| `Environment` | Habitat (e.g. `halite`) |
| `Latitude_N` | Latitude |
| `Longitude_W` | Longitude |
| `Dataset` | `AD` or `CB` |

---

## Running the Pipeline

All scripts are designed for SLURM (`sbatch`). Adjust the `#SBATCH` headers for your cluster. Replace `/home/apolonio/orcd/scratch/20_440` with your base directory throughout.

```bash
# 1. Download all reads
sbatch 01_download.sh AD ALL
sbatch 01_download.sh CB ALL

# 2. QC, assembly, annotation
sbatch 02_kneaddata_assembly.sh AD
sbatch 02_kneaddata_assembly.sh CB

# 3. Functional profiling
sbatch 03_humann.sh AD
sbatch 03_humann.sh CB

# 4. sRNA detection (both timegroups per dataset)
sbatch 04_snapt.sh AD 2016
sbatch 04_snapt.sh AD 2017
sbatch 04_snapt.sh CB Fall
sbatch 04_snapt.sh CB Sum

# 5. BAM alignment + secondary structure + interactions
sbatch 05b_bam.sh
sbatch 05_rnafold.sh AD 2016 && sbatch 05_rnafold.sh AD 2017
sbatch 06_intarna.sh AD 2016 && sbatch 06_intarna.sh AD 2017

# 6. Gene → UniRef90 mapping
sbatch 09_gene_to_uniref.sh

# 7. Feature matrix
sbatch 07_feature_matrix.sh AD
sbatch 07_feature_matrix.sh CB

# 8. Core analysis (per-replicate TPM first for full correlation)
sbatch 08b_per_replicate_tpm.sh AD 2016
sbatch 08b_per_replicate_tpm.sh AD 2017
sbatch 08_run_analysis.sh AD 2016 2017
sbatch 08_run_analysis.sh CB Fall Sum

# 9. Tables + figures
sbatch 10_table_generation.sh
python 11_plot_figures.py
```

---

## Key Outputs

| File | Description |
|------|-------------|
| `<DATASET>/feature_matrix.tsv` | Random Forest input features |
| `<DATASET>/analysis_output_*/` | DE results, RF model, correlation tables, summary plots |
| `<DATASET>/per_replicate_tpm_*/tpm_matrix.tsv` | Gene × replicate TPM matrix |
| `gene_to_uniref90.tsv` | Prodigal proteins → UniRef90 hits |
| `srna_master.tsv` | One row per sRNA × sample |
| `srna_target_pathway.tsv` | sRNA → target gene → pathway links |
| `pathway_master.tsv` | One row per pathway × sample |
| `srna_presence.tsv` | sRNA presence/absence per sample |
| `mechanism_table.tsv` | sRNA–pathway pairs with consistency flag |

---

## Script Reference

| # | Script | Language | Purpose |
|---|--------|----------|---------|
| 01 | `01_download.sh` | Bash | SRA download |
| 02 | `02_kneaddata_assembly.sh` | Bash | QC, assembly, annotation, indexing |
| 03 | `03_humann.sh` | Bash | HUMAnN functional profiling |
| 04 | `04_snapt.sh` | Bash | SnapT sRNA detection |
| 05 | `05_rnafold.sh` | Bash + Python | RNAfold secondary structure |
| 05b | `05b_bam.sh` | Bash | Bowtie2 alignment → BAM |
| 06 | `06_intarna.sh` | Bash + Python | IntaRNA interaction prediction |
| 07 | `07_feature_matrix.sh` | Bash | Feature matrix  |
| 07 | `07_build_feature_matrix.py` | Python | Feature matrix construction |
| 08 | `08_run_analysis.sh` | Bash | Random Forest snalysis  |
| 08 | `08_srna_analysis.py` | Python | DE, RF, correlation analysis |
| 08b | `08b_per_replicate_tpm.sh` | Bash + Python | Per-replicate StringTie TPM |
| 09 | `09_gene_to_uniref.sh` | Bash | DIAMOND blastp → UniRef90 |
| 10 | `10_table_generation.sh` | Bash | Table generation  |
| 10 | `10_table_generation.py` | Python | Master table generation |
| 11 | `11_plot_figures.py` | Python | Figure generation |
