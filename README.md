# sRNA_20440
Code used in the 20.440 class project: Computational Prediction of Community Abundance Dynamics Using Insights From Differential sRNA Expression  

# sRNA Regulatory Signals and Metabolic Pathway Shifts in Microbial Communities

## Repository Overview

This repository contains a complete computational pipeline for discovering small regulatory RNAs (sRNAs) in metagenomic assemblies, predicting their mRNA targets, and linking sRNA expression changes to metabolic pathway reorganization in response to environmental perturbations.

**Citation:** [To be added upon publication]

---

## Pipeline Architecture

```
Raw Reads (DNA + RNA)
    ↓
[02] Download from SRA
    ↓
[03] Quality Control (KneadData)
    ↓
[04] Assembly & Annotation (MEGAHIT, Prodigal, Bowtie2)
    ↓
[05] sRNA Discovery (SnapT) + Structure (RNAfold)
    ↓
[06] Target Prediction (IntaRNA)
    ↓
[07] Feature Matrix Construction
    ↓
[08] sRNA Analysis (Differential Expression, RF)
    ↓
[09] Metagenome Verification (Beta Diversity, Pathways)
    ↓
[10] sRNA-Pathway Integration
    ↓
[11] Validation (Mechanistic Hypotheses)
```

---

## Script Documentation

### **02_download.sh** – Download Sequencing Data from SRA

**Purpose:** Retrieve raw sequencing reads from NCBI Sequence Read Archive (SRA).

**Dependencies:** 
- SRA Toolkit v3.2.1 (`prefetch`, `fasterq-dump`)

**Usage:**
```bash
sbatch 02_download.sh [DATASET] [DATA_TYPE]

# Examples:
sbatch 02_download.sh AD ALL              # All Atacama samples
sbatch 02_download.sh CB Metatranscriptome # Only Chesapeake RNA
sbatch 02_download.sh ALL ALL             # Everything
```

**Parameters:**
- `DATASET`: AD (Atacama Desert) | CB (Chesapeake Bay) | ALL
- `DATA_TYPE`: Metagenome | Metatranscriptome | ALL

**Inputs:**
- `Metadata.csv`: Sample metadata with columns: bioproject, srr (SRA Run accession), type (Metagenome/Metatranscriptome), name (sample ID), dataset (AD/CB)

**Outputs:**
- `/{DATASET}/{SAMPLE_NAME}/*.fastq.gz`: Paired-end reads (R1, R2), gzip-compressed
- All raw data organized by dataset and sample name for downstream processing

**Key Features:**
- Parameterized to download subsets (single dataset, single data type, or all)
- Renames files from SRR accessions to human-readable sample names
- Skips already-downloaded samples (idempotent)
- Compresses output to save disk space

---

### **03_kneaddata_assembly.sh** – Quality Control and Read Filtering

**Purpose:** Remove host and ribosomal contamination; prepare reads for assembly.

**Dependencies:**
- KneadData v0.7.4
- Bowtie2 v2.4.2
- Reference databases:
  - Human genome (bowtie2 index): `human_dna_db/`
  - Human RNA (bowtie2 index): `human_rna_db/`
  - Ribosomal RNA (bowtie2 index): `rrna_db/`

**Usage:**
```bash
sbatch 03_kneaddata_assembly.sh [DATASET]

# Examples:
sbatch 03_kneaddata_assembly.sh AD
sbatch 03_kneaddata_assembly.sh CB
```

**Parameters:**
- `DATASET`: AD | CB
- `THREADS`: 8 (per-sample; total parallelization handled by SLURM)

**Inputs:**
- `/{DATASET}/{SAMPLE}/*_1.fastq.gz`: Raw R1 reads
- `/{DATASET}/{SAMPLE}/*_2.fastq.gz`: Raw R2 reads

**Outputs:**
- `/{DATASET}/{SAMPLE}_knead_out/*_kneaddata_paired_1.fastq`: Cleaned R1 reads
- `/{DATASET}/{SAMPLE}_knead_out/*_kneaddata_paired_2.fastq`: Cleaned R2 reads
- `/{DATASET}/{SAMPLE}_knead_out/*.log`: KneadData logs (contamination rates)

**Key Features:**
- Metagenomes: Filters human DNA contamination
- Metatranscriptomes: Removes both human RNA and ribosomal RNA (rRNA/tRNA)
- Removes intermediate bowtie2 alignments to save disk
- Outputs uncompressed (for faster downstream access)

**Output Statistics Expected:**
- Metagenomes: ~92% reads retained (8% human DNA)
- Metatranscriptomes: ~12% reads retained (88% rRNA)

---

### **04_assembly.sh** – Metagenome Co-assembly and Gene Annotation

**Purpose:** Assemble metagenomic reads into contigs and annotate genes.

**Dependencies:**
- MEGAHIT v1.2.9
- Prodigal v2.6.3
- Bowtie2 v2.4.2

**Usage:**
```bash
sbatch 04_assembly.sh [DATASET]

# Examples:
sbatch 04_assembly.sh AD
sbatch 04_assembly.sh CB
```

**Parameters:**
- `DATASET`: AD | CB
- `MEGAHIT` settings:
  - `--min-contig-len 500`: Minimum contig length (balances sensitivity/specificity)
  - `-t 48`: Threads
- `PRODIGAL` mode: `-p meta` (metagenomic mode)

**Inputs:**
- `/{DATASET}/{SAMPLE}_knead_out/*_kneaddata_paired_1.fastq`: All cleaned R1 reads (pooled)
- `/{DATASET}/{SAMPLE}_knead_out/*_kneaddata_paired_2.fastq`: All cleaned R2 reads (pooled)

**Outputs:**
- `/{DATASET}/coassembly/final.contigs.fa`: Co-assembled metagenome (FASTA)
- `/{DATASET}/coassembly/final.contigs.clean.fa`: Header-sanitized version (removes whitespace)
- `/{DATASET}/coassembly/bowtie2_index.*`: Bowtie2 index files (6 files)
- `/{DATASET}/annotation/prodigal.gff`: Gene annotations (GFF3 format)
- `/{DATASET}/annotation/prodigal.faa`: Protein sequences (FASTA)

**Key Features:**
- Co-assembly pools reads across replicates for improved contiguity
- MEGAHIT automatically determines optimal kmer sizes
- Prodigal uses metagenomic training models for gene prediction
- Outputs both standard and sanitized FASTA versions for compatibility with different tools

**Quality Metrics:**
- Assembly statistics logged in MEGAHIT output
- Expected N50: 5–50 kb (depending on community complexity)
- Expected gene count: 500k–1M+ genes (depending on assembly size)

---

### **04_snapt.sh** – sRNA Discovery and Expression Quantification

**Purpose:** Identify sRNAs in assembly and quantify their expression across samples.

**Dependencies:**
- SnapT v1.0 (Gelsinger et al. 2020)
- HISAT2 v2.2.1
- StringTie v2.1.7
- BLAST (nr database): `/home/apolonio/orcd/scratch/databases/nr/nr.dmnd`
- Rfam (covariance models): `/home/apolonio/orcd/scratch/databases/rfam/Rfam.cm`

**Usage:**
```bash
sbatch 04_snapt.sh [DATASET] [TIMEGROUP]

# Examples (Atacama):
sbatch 04_snapt.sh AD 2016
sbatch 04_snapt.sh AD 2017

# Examples (Chesapeake):
sbatch 04_snapt.sh CB Fall
sbatch 04_snapt.sh CB Sum
```

**Parameters:**
- `DATASET`: AD | CB
- `TIMEGROUP`: 
  - Atacama: 2016, 2017 (year of sampling)
  - Chesapeake: Fall, Sum (season)
- SnapT parameters:
  - `-r RF`: Use Rfam classification mode
  - `-t 48`: Threads
  - `-D nr_db`: NCBI NR BLAST database
  - `-R rfam_db`: Rfam covariance models

**Inputs:**
- `/{DATASET}/coassembly/final.contigs.clean.fa`: Assembly (from script 04)
- `/{DATASET}/annotation/prodigal.gff`: Gene annotations (from script 04)
- `/{DATASET}/{SAMPLE}_knead_out/*_kneaddata_paired_1.fastq`: All cleaned RNA R1 (for this timegroup)
- `/{DATASET}/{SAMPLE}_knead_out/*_kneaddata_paired_2.fastq`: All cleaned RNA R2 (for this timegroup)
- `/{DATASET}/combined_reads/`: Directory for pooled reads

**Outputs:**
- `/{DATASET}/snapt_{TIMEGROUP}/small_antisense_ncRNAs.gff`: Annotated asRNA (GFF3)
- `/{DATASET}/snapt_{TIMEGROUP}/small_intergenic_ncRNAs.gff`: Annotated itsRNA (GFF3)
- `/{DATASET}/snapt_{TIMEGROUP}/small_ncRNAs.gff`: Combined sRNA annotations
- `/{DATASET}/snapt_{TIMEGROUP}/small_antisense_ncRNAs.fa`: asRNA sequences (FASTA)
- `/{DATASET}/snapt_{TIMEGROUP}/small_intergenic_ncRNAs.fa`: itsRNA sequences (FASTA)
- `/{DATASET}/combined_reads/{DATASET}_{TIMEGROUP}_R1.fastq`: Pooled R1 reads (temporary)
- `/{DATASET}/combined_reads/{DATASET}_{TIMEGROUP}_R2.fastq`: Pooled R2 reads (temporary)

**Key Features:**
- Pools all RNA reads for a timegroup before sRNA discovery
- Identifies both asRNA (antisense) and itsRNA (intergenic) classes
- Filters out known non-regulatory RNAs (rRNA, tRNA) via Rfam
- Quantifies sRNA expression (FPKM, TPM) via StringTie
- Handles multiple timegroups within a dataset independently

**Output Format (GFF3):**
```
k141_251874	SnapT	sRNA	708	858	.	+	.	ID=srna_001;transcript_id=SNAPT_asRNA_1;FPKM=45.3;TPM=120.5
```

---

### **05_rnafold.sh** – RNA Secondary Structure Prediction

**Purpose:** Predict secondary structure and thermodynamic stability of sRNAs.

**Dependencies:**
- RNAfold v2.4.18 (ViennaRNA package)
- GNU Parallel

**Usage:**
```bash
sbatch 05_rnafold.sh [DATASET] [TIMEGROUP]

# Examples:
sbatch 05_rnafold.sh AD 2016
sbatch 05_rnafold.sh CB Fall
```

**Parameters:**
- `DATASET`: AD | CB
- `TIMEGROUP`: Year/season matching prior SnapT run
- `THREADS`: 36 (parallel RNAfold jobs)
- RNAfold options:
  - `--noPS`: Skip PostScript output (saves disk)
  - `--partfunc`: Compute partition function (ensemble thermodynamics)

**Inputs:**
- `/{DATASET}/snapt_{TIMEGROUP}/small_antisense_ncRNAs.fa`: asRNAs (from script 04)
- `/{DATASET}/snapt_{TIMEGROUP}/small_intergenic_ncRNAs.fa`: itsRNAs (from script 04)

**Outputs:**
- `/{DATASET}/rnafold_{TIMEGROUP}/rnafold_summary.tsv`: Summary table (TSV)
  - Columns: srna_id, length, gc_content, mfe, ensemble_energy, struct_entropy, structure
- `/{DATASET}/rnafold_{TIMEGROUP}/individual_seqs/`: Individual RNAfold output files (one per sRNA)
- `/{DATASET}/rnafold_{TIMEGROUP}/results/`: Parsed structure predictions

**Output Format (rnafold_summary.tsv):**
```
srna_id	length	gc_content	mfe	ensemble_energy	struct_entropy	structure
k141_251874:708-858(+)	150	0.520	-69.60	-66.40	0.342	((((((((((.(((((((((...))))))))))...)))))))).....))
```

**Metrics Computed:**
- **gc_content**: (G + C) / length, normalized to [0, 1]
- **mfe**: Minimum free energy (kcal/mol), negative = more stable
- **ensemble_energy**: Weighted average energy across all possible structures
- **struct_entropy**: Shannon entropy of base-pairing probability, [0, 1]
- **structure**: Secondary structure in dot-bracket notation

**Key Features:**
- Parallelizes across all sRNAs
- Handles both asRNA and itsRNA combined
- Computes both MFE and ensemble properties for robustness

---

### **06_intarna.sh** – sRNA-mRNA Target Prediction

**Purpose:** Identify putative mRNA targets for each sRNA via RNA-RNA binding prediction.

**Dependencies:**
- IntaRNA v2.2.2

**Usage:**
```bash
sbatch 06_intarna.sh [DATASET] [TIMEGROUP] [PVALUE_THRESHOLD]

# Examples:
sbatch 06_intarna.sh AD 2016 0.05
sbatch 06_intarna.sh CB Fall 0.05
```

**Parameters:**
- `DATASET`: AD | CB
- `TIMEGROUP`: Year/season matching prior runs
- `PVALUE_THRESHOLD`: Default 0.05 (adjust for sensitivity/specificity trade-off)
- IntaRNA options:
  - `--outMaxE 0`: Only report favorable interactions (ΔG < 0)
  - `--outDeltaE 100`: Report up to 100 kcal/mol less stable than optimal
  - `--noSeed`: Allow non-seed interactions
  - `-n 100`: Top 100 interactions per sRNA-target pair
  - `--threads 48`: Parallel threads

**Inputs:**
- `/{DATASET}/snapt_{TIMEGROUP}/small_antisense_ncRNAs.fa`: asRNAs
- `/{DATASET}/snapt_{TIMEGROUP}/small_intergenic_ncRNAs.fa`: itsRNAs
- `/{DATASET}/annotation/prodigal.gff`: Gene annotations (for CDS extraction)
- `/{DATASET}/coassembly/final.contigs.fa`: Assembly

**Outputs:**
- `/{DATASET}/intarna_{TIMEGROUP}/cds_regions.bed`: All CDS coordinates (BED format)
- `/{DATASET}/intarna_{TIMEGROUP}/cds_filtered.bed`: Expressed CDS only (first 150 bp)
- `/{DATASET}/intarna_{TIMEGROUP}/mrna_targets.fa`: mRNA target sequences (FASTA)
- `/{DATASET}/intarna_{TIMEGROUP}/interactions.tsv`: All IntaRNA predictions (raw)
- `/{DATASET}/intarna_{TIMEGROUP}/interactions_significant.tsv`: Filtered by p-value (TSV)
- `/{DATASET}/intarna_{TIMEGROUP}/intarna_summary.tsv`: Per-sRNA summary statistics (TSV)
- `/{DATASET}/intarna_{TIMEGROUP}/intarna_summary_fixed.tsv`: Header-corrected version (for downstream)

**Output Format (interactions_significant.tsv):**
```
srna_id	mrna_id	energy	ED1	ED2	Pu1	Pu2	hybridDP	start1	end1	start2	end2	pvalue
k141_527369:8116-8391(+)	ID=10359_1...::k141_251874:708-858(+)	-17.54	3.57	13.78	0.00305	1.95e-10	((((...))))	3	34	190	241	0.0133
```

**Output Format (intarna_summary.tsv):**
```
srna_id	n_targets	best_energy	mean_energy	std_energy	n_interactions	best_pvalue	mean_pvalue
k141_527369:8116-8391(+)	128	-17.54	-4.49	2.77	12017	1.57e-08	0.00168
```

**Key Features:**
- Filters targets to expressed genes (coverage > 0 in metatranscriptome)
- Limits targets to first 150 bp of CDS (regulatory region)
- Computes both raw and filtered interaction sets
- Provides per-sRNA statistics for quick reference

**Key Columns:**
- **energy**: Interaction free energy (kcal/mol); more negative = stronger binding
- **ED1, ED2**: Donor/acceptor interaction sites
- **Pu1, Pu2**: Target accessibility probabilities (0–1)
- **pvalue**: Statistical significance (logistic regression on sequence features)

---

### **07_build_feature_matrix.py** – Feature Integration and Aggregation

**Purpose:** Merge sRNA, structural, interaction, and community data into unified feature matrices for statistical analysis.

**Dependencies:**
- Python 3.6+
- pandas, numpy

**Usage:**
```bash
python3 07_build_feature_matrix.py \
    --dataset AD \
    --metadata Metadata.csv \
    --snapt AD/snapt_2016/small_antisense_ncRNAs.gff \
    --rnafold AD/rnafold_2016/rnafold_summary.tsv \
    --humann_rna AD/Total_RNA_pathabundance_humann_table.tsv \
    --humann_dna AD/Total_DNA_pathabundance_humann_table.tsv \
    --metaphlan AD/Total_merged_abundance_table.tsv \
    --output AD/feature_matrix.tsv \
    --intarna AD/intarna_2016/intarna_summary_fixed.tsv
```

**Required Arguments:**
- `--dataset`: AD | CB
- `--metadata`: Sample metadata CSV
- `--snapt`: Path to any SnapT GFF (used to locate dataset directory)
- `--rnafold`: RNAfold summary TSV (first timegroup)
- `--humann_rna`: HUMAnN RNA pathabundance table
- `--humann_dna`: HUMAnN DNA pathabundance table
- `--metaphlan`: MetaPhlAn species abundance table
- `--output`: Output TSV path

**Optional Arguments:**
- `--intarna`: IntaRNA summary TSV (optional; skipped if missing)

**Inputs:**
- All files from scripts 04, 05, 06, plus external HUMAnN/MetaPhlAn outputs

**Outputs:**
- `/{DATASET}/feature_matrix.tsv`: Timepoint-level features (1 row per timepoint pair)
- `/{DATASET}/feature_matrix_species_level.tsv`: Species-level features (N_species × N_timepoint_pairs)
- `/{DATASET}/feature_matrix_summary.json`: Summary metadata (features used, timepoints, data quality)

**Output Format (feature_matrix.tsv, excerpt):**
```
dataset	timepoint_T1	timepoint_T2	pair_label	srna_mean_length	asrna_fraction	mean_srna_log2fc_fpkm	n_srnas_upregulated	n_srnas_downregulated	mean_gc_content	mean_mfe	srna_mean_n_targets	rna_dna_ratio_T1	log2fc_rna_activity	target_mean_log2fc	n_species_total
AD	2016	2017	2016_vs_2017	218.5	0.612	-8.32	9	11	0.548	-65.3	2.8	0.042	-0.87	-2.14	7
```

**Output Format (feature_matrix_species_level.tsv, excerpt):**
```
dataset	pair_label	timepoint_T1	timepoint_T2	species	abundance_T1	abundance_T2	log2fc	mean_srna_log2fc_fpkm	asrna_fraction	n_srnas_upregulated	mean_mfe
AD	2016_vs_2017	2016	2017	Connexibacter disulfidoxidans	0.021	0.008	-1.39	-8.32	0.612	9	-65.3
```

**Key Features:**
- Normalizes all metrics (z-score transformation)
- Computes species-level residual abundance changes (removes baseline bias)
- Fills missing IntaRNA data gracefully
- Provides JSON summary for reproducibility tracking

---

### **08_srna_analysis.py** – sRNA Differential Expression and Classification

**Purpose:** Characterize sRNA expression changes, structural heterogeneity, and build predictive models.

**Dependencies:**
- Python 3.6+
- pandas, numpy, matplotlib, scipy, scikit-learn, shap

**Usage:**
```bash
python3 08_srna_analysis.py \
    --dataset AD \
    --snapt_t1 AD/snapt_2016 \
    --snapt_t2 AD/snapt_2017 \
    --rnafold_t1 AD/rnafold_2016/rnafold_summary.tsv \
    --rnafold_t2 AD/rnafold_2017/rnafold_summary.tsv \
    --feature_matrix AD/feature_matrix.tsv \
    --outdir AD/analysis_output_2016_2017 \
    --intarna_t1 AD/intarna_2016/intarna_summary_fixed.tsv \
    --intarna_t2 AD/intarna_2017/intarna_summary_fixed.tsv \
    --label_t1 2016 \
    --label_t2 2017
```

**Required Arguments:**
- `--dataset`: AD | CB
- `--snapt_t1`, `--snapt_t2`: Directories with SnapT output (timepoint 1, 2)
- `--rnafold_t1`, `--rnafold_t2`: RNAfold summary TSVs
- `--feature_matrix`: Feature matrix TSV (from script 07)
- `--outdir`: Output directory

**Optional Arguments:**
- `--intarna_t1`, `--intarna_t2`: IntaRNA summaries (optional)
- `--label_t1`, `--label_t2`: Custom timepoint labels

**Inputs:**
- sRNA annotations + expressions (SnapT GFF)
- Secondary structure metrics (RNAfold)
- Target predictions (IntaRNA)
- Feature matrix (script 07)

**Outputs:**
- `srna_expression_log2fc.tsv`: Per-sRNA expression changes (log2FC), detection status
- `asrna_vs_itsrna_structural.tsv`: Mann-Whitney U test results (structure comparison)
- `srna_complete_feature_table.tsv`: All features per sRNA
- `cv_scores_srna.tsv`: Cross-validation scores for sRNA-level RF
- `shap_importance_srna.tsv`: SHAP feature importance rankings
- `volcano_srna.pdf/png`: MA plot of sRNA expression
- `srna_turnover.pdf/png`: Pool composition and turnover visualization
- `struct_comparison.pdf/png`: Boxplots of structural features by sRNA class
- `srna_class_composition.pdf/png`: Pie charts of asRNA vs itsRNA composition
- `srna_log2fc_ranked.pdf/png`: Ranked barplot of sRNA log2FC
- `srna_fpkm_distributions.pdf/png`: Expression distributions by class and timepoint
- `shap_bar_srna.pdf/png`: SHAP importance barplot
- `shap_summary_srna.pdf/png`: SHAP summary plot

**Key Analyses:**

1. **Differential Expression:**
   - log2FC computed with pseudocount (1e-6)
   - Permutation test for mean FC significance

2. **Structural Comparison:**
   - Mann-Whitney U test (asRNA vs itsRNA)
   - Metrics: GC content, MFE, structure entropy

3. **sRNA-Level Random Forest:**
   - Predicts up/downregulation from structural + interaction features
   - 5-fold CV or LOO CV (data-dependent)
   - SHAP explainability included

---

### **09_metagenome_verification.py** – Independent Community-Level Validation

**Purpose:** Validate sRNA findings using independent metagenomic data (species composition, pathway shifts, beta diversity).

**Dependencies:**
- Python 3.6+
- pandas, numpy, matplotlib, scipy

**Usage:**
```bash
python3 09_metagenome_verification.py \
    --metaphlan AD/Total_merged_abundance_table.tsv \
    --humann_dna AD/Total_DNA_pathabundance_humann_table.tsv \
    --humann_rna AD/Total_RNA_pathabundance_humann_table.tsv \
    --srna_log2fc AD/analysis_output_2016_2017/srna_expression_log2fc.tsv \
    --outdir AD/metagenome_verification \
    --dataset AD \
    --t1 2016 \
    --t2 2017
```

**Required Arguments:**
- `--metaphlan`: MetaPhlAn abundance table
- `--humann_dna`: HUMAnN DNA pathway abundance
- `--humann_rna`: HUMAnN RNA pathway abundance
- `--srna_log2fc`: sRNA log2FC table (from script 08)
- `--outdir`: Output directory

**Optional Arguments:**
- `--dataset`: Dataset label (for plots)
- `--t1`, `--t2`: Timepoint labels

**Outputs:**
- `species_abundance_log2fc.tsv`: Per-species abundance changes
- `bray_curtis_matrix.tsv`: Pairwise Bray-Curtis distances (species)
- `pathway_log2fc.tsv`: Per-pathway RNA/DNA log2FC
- `beta_diversity.pdf/png`: Boxplot of within vs. between timepoint distances
- `pathway_rna_vs_dna_log2fc.pdf/png`: Scatter of RNA vs DNA pathway changes
- `convergent_evidence.pdf/png`: 3-panel summary figure

**Key Analyses:**

1. **Species Composition:**
   - Jaccard similarity (presence/absence overlap)
   - Comparison with sRNA Jaccard

2. **Beta Diversity:**
   - Bray-Curtis distances (within timepoint vs. between)
   - Mann-Whitney U test

3. **Pathway Activity:**
   - RNA vs DNA pathway log2FC comparison
   - Spearman correlation

---

### **10_srna_target_pathway_analysis.py** –  Integration of sRNA + Metabolic pathways

**Purpose:** Link sRNA targets to metabolic pathways; test mechanistic hypothesis.

**Dependencies:**
- Python 3.6+
- pandas, numpy, matplotlib, networkx, scipy

**Usage:**
```bash
python3 10_srna_target_pathway_analysis.py \
    --dataset AD \
    --srna_log2fc AD/analysis_output_2016_2017/srna_expression_log2fc.tsv \
    --intarna_summary AD/intarna_2016/intarna_summary_fixed.tsv \
    --intarna_targets AD/intarna_2016/interactions_significant_fixed.tsv \
    --prodigal_gff AD/annotation/prodigal.gff \
    --pathways AD/metagenome_verification/pathway_log2fc.tsv \
    --humann_dna AD/Total_DNA_pathabundance_humann_table.tsv \
    --humann_rna AD/Total_RNA_pathabundance_humann_table.tsv \
    --metadata Metadata.csv \
    --outdir AD/srna_target_pathway_analysis \
    --t1 2016 \
    --t2 2017
```

**Required Arguments:**
- `--dataset`: AD | CB
- `--srna_log2fc`: sRNA expression table
- `--intarna_summary`, `--intarna_targets`: IntaRNA outputs
- `--prodigal_gff`: Gene annotations
- `--pathways`: Pathway log2FC (from script 09)
- `--humann_dna`, `--humann_rna`: Pathway abundances
- `--metadata`: Sample metadata
- `--outdir`: Output directory

**Outputs:**
- `srna_target_summary.tsv`: Top N sRNAs + target counts + interaction energies
- `top_srna_target_interactions.tsv`: Detailed interactions (sRNA → gene → energy)
- `srna_target_network.pdf/png`: Network graph (sRNA → targets)
- `srna_target_pathway_summary.pdf/png`: sRNA changes vs. target counts, pathway shifts
- `top_srnas_target_counts.tsv`: Quick reference table

**Output Format (srna_target_summary.tsv):**
```
srna_id	log2fc_fpkm	srna_type	n_predicted_targets	interaction_energy_best	interaction_energy_mean
k141_728975:4614-4759(+)	-33.64	intergenic	100	-20.59	-16.56
```

---

### **11_validate_srna_pathway_links.py** –  Hypothesis Validation

**Purpose:** Test whether sRNA targeting patterns support predicted mechanistic model.

**Dependencies:**
- Python 3.6+
- pandas, numpy, matplotlib, scipy

**Usage:**
```bash
python3 11_validate_srna_pathway_links.py \
    --dataset AD \
    --srna_log2fc AD/analysis_output_2016_2017/srna_expression_log2fc.tsv \
    --intarna_targets AD/intarna_2016/interactions_significant_fixed.tsv \
    --pathways AD/metagenome_verification/pathway_log2fc.tsv \
    --humann_genefamilies AD/Total_DNA_genefamilies_humann_table.tsv \
    --outdir AD/validation_results
```

**Required Arguments:**
- `--dataset`: AD | CB
- `--srna_log2fc`: sRNA expression table
- `--intarna_targets`: IntaRNA predictions
- `--pathways`: Pathway log2FC
- `--humann_genefamilies`: Gene family abundances
- `--outdir`: Output directory

**Outputs:**
- `validation_results.tsv`: Validation metrics for N = 10, 20, 30, 50 sRNAs
- `validation_results.pdf/png`: 4-panel summary figure

**Output Format (validation_results.tsv):**
```
n_srnas	n_upregulated	n_downregulated	downreg_interactions	upreg_interactions	downreg_unique_genes	upreg_unique_genes	downreg_srnas_with_targets	upreg_srnas_with_targets	target_ratio	hypothesis_1_supported	hypothesis_2_supported
10	4	6	300	0	14	0	6	0	14.0	True	True
20	9	11	550	0	26	0	11	0	26.0	True	True
```

**Key Metrics:**

- **H1 (Downreg Coverage)**: % of downregulated sRNAs with ≥1 target (expected 100%)
- **H2 (Selectivity)**: Ratio of downregulated sRNA targets to upregulated sRNA targets (expected > 1.5)

---

## Installation & Setup

### **Dependencies**

Install all required tools via conda:

```bash
conda create -n srna_pipeline -c bioconda -c conda-forge \
    sratoolkit kneaddata megahit prodigal hisat2 stringtie snapt \
    rnafold intarna bowtie2 humann metaphlan matplotlib scipy scikit-learn shap networkx

conda activate srna_pipeline
```

### **Database Setup**

Download reference databases (requires ~100 GB disk space):

```bash
# Create database directory
mkdir -p /path/to/databases

# NCBI NR (for IntaRNA/SnapT)
cd /path/to/databases
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nr.*.tar.gz
tar -xzf nr.*.tar.gz
diamond makedb --in nr.faa -d nr.dmnd

# Rfam (for SnapT sRNA annotation)
wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
gunzip Rfam.cm.gz

# Human genome (for contamination filtering)
mkdir human_dna_db
bowtie2-build hg38.fa human_dna_db/hg_39

# Human RNA (for metatranscriptome QC)
mkdir human_rna_db
bowtie2-build human_mRNA.fa human_rna_db/human_hg38_refMrna

# Ribosomal RNA (SILVA)
mkdir rrna_db
bowtie2-build SILVA_128_SSUParc_LSUParc.fna rrna_db/SILVA_128
```

Update database paths in all scripts before running.

### **Directory Structure**

Expected input organization:

```
20_440/
├── Metadata.csv                    # Sample metadata
├── 02_download.sh
├── 03_kneaddata_assembly.sh
├── 04_assembly.sh
├── 04_snapt.sh
├── 05_rnafold.sh
├── 06_intarna.sh
├── 07_build_feature_matrix.py
├── 08_srna_analysis.py
├── 09_metagenome_verification.py
├── 10_srna_target_pathway_analysis.py
├── 11_validate_srna_pathway_links.py
│
├── AD/                             # Atacama dataset
│   ├── AD_S1_2016_02_R1_DNA/       # Downloaded SRA
│   ├── coassembly/
│   ├── annotation/
│   ├── snapt_2016/
│   ├── snapt_2017/
│   ├── rnafold_2016/
│   ├── rnafold_2017/
│   ├── intarna_2016/
│   ├── intarna_2017/
│   ├── analysis_output_2016_2017/
│   ├── metagenome_verification/
│   ├── srna_target_pathway_analysis/
│   └── validation_results/
│
└── CB/                             # Chesapeake Bay dataset
    ├── coassembly/
    ├── annotation/
    └── [similar subdirectories...]
```

---

## Running the Complete Pipeline

### **Option 1: Run Full Pipeline Sequentially**

```bash
# 1. Download data
sbatch 02_download.sh AD ALL

# 2. Quality control (run both DNA and RNA)
sbatch 03_kneaddata_assembly.sh AD

# 3. Assembly and gene annotation
sbatch 04_assembly.sh AD

# 4. sRNA discovery (run both timegroups)
sbatch 04_snapt.sh AD 2016
sbatch 04_snapt.sh AD 2017

# 5. Secondary structure
sbatch 05_rnafold.sh AD 2016
sbatch 05_rnafold.sh AD 2017

# 6. Target prediction
sbatch 06_intarna.sh AD 2016
sbatch 06_intarna.sh AD 2017

# 7. Feature integration (in Python)
python3 07_build_feature_matrix.py --dataset AD ...

# 8-11. Analyses (Python)
python3 08_srna_analysis.py --dataset AD ...
python3 09_metagenome_verification.py --dataset AD ...
python3 10_srna_target_pathway_analysis.py --dataset AD ...
python3 11_validate_srna_pathway_links.py --dataset AD ...
```


## Expected Outputs Summary

| Script | Key Outputs | File Size | Runtime |
|--------|------------|-----------|---------|
| 02 | Raw FASTQ files | 50–100 GB | 2–4 hrs |
| 03 | Cleaned FASTQ | 5–10 GB | 4–6 hrs |
| 04 | Assembly, genes, index | 1 GB | 8–12 hrs |
| 04_snapt | sRNA annotations | 100 MB | 12–24 hrs |
| 05 | Structures | 50 MB | 2–3 hrs |
| 06 | Interactions | 500 MB–1 GB | 4–8 hrs |
| 07 | Feature matrix | 10 MB | <1 min |
| 08 | DE analysis, figures | 50 MB | 5–10 min |
| 09 | Community validation | 30 MB | 2–5 min |
| 10 | Network, integration | 100 MB | 5–15 min |
| 11 | Validation metrics | 10 MB | <1 min |

---

## Troubleshooting

### **Common Issues**

**Out of memory errors (MEGAHIT, IntaRNA):**
- Reduce thread count (e.g., `-t 24` instead of 48)
- Increase allocated memory via `#SBATCH --mem=300GB`
- The CB samples were quite large so it's a good idea to regularly check storage space while things are running 

**Missing output files:**
- Check SLURM log files: `cat project_logs/slurm-[jobid].log`
- Verify input files exist before running downstream scripts
- Check for whitespace in FASTA headers (use `sed 's/>\(\S*\).*/>\1/' input.fa > output.fa`)

**IntaRNA not finding interactions:**
- Verify target sequences are >20 nt (minimum for meaningful binding prediction)
- Check that sRNA sequences are present in FASTA file
- Lower p-value threshold if too stringent (try 0.1 instead of 0.05)

**Script dependency issues:**
- Ensure Python scripts are executable: `chmod +x *.py`
- Activate conda environment before use: `conda activate srna_pipeline`
- Check Python version: `python3 --version` (should be 3.6+)

---
