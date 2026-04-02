#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH -n 48
#SBATCH -J kneaddata_assembly
#SBATCH -p mit_normal
#SBATCH --mem=300GB

set -euo pipefail

module load miniforge
conda activate snapt_env

# Directories
METADATA="/home/apolonio/orcd/scratch/Metadata.csv"
BASE_DIR="/home/apolonio/orcd/scratch"

HOST_REF="/home/apolonio/human_dna_db"
RNA_DB="/home/apolonio/human_rna_db"
RRNA_DB="/home/apolonio/rrna_db"
THREADS=8      # per-sample KneadData — DO NOT exceed 8

# Dataset flag — override at submission:
#   sbatch 03_kneaddata_assembly.sh AD
#   sbatch 03_kneaddata_assembly.sh CB
DATASET="${1:-AD}"

echo "===== KneadData + Assembly: $DATASET ====="

DATASET_DIR="${BASE_DIR}/${DATASET}"
ASSEMBLY_DIR="${DATASET_DIR}/coassembly"
ANNOTATION_DIR="${DATASET_DIR}/annotation"

#  Helper: list sample Names for a dataset 
get_samples() {
    local ds="$1"
    local tp="$2"
    tail -n +2 "$METADATA" \
        | awk -F'\t' -v ds="$ds" -v tp="$tp" \
            '$3==tp && $9==ds { gsub(/\r/,"",$4); print $4 }'
}

# STEP 1 — KneadData

echo " STEP 1: Clean-Up "

run_kneaddata() {
    local sample="$1"
    local type="$2"
    local sample_dir="${DATASET_DIR}/${sample}"
    local outdir="${DATASET_DIR}/${sample}_knead_out"

    if ls "${outdir}"/*_kneaddata_paired_1.fastq >/dev/null 2>&1; then
        rm -rf "$sample_dir"
    else
        echo "  WARNING: KneadData output missing for $sample"
    fi
}


# STEP 2 — MEGAHIT co-assembly (DNA only)

echo " STEP 2: MEGAHIT "

if [ -f "${ASSEMBLY_DIR}/final.contigs.fa" ]; then
    echo "  Assembly exists — skipping"
else
    R1_LIST=""; R2_LIST=""
    while IFS= read -r s; do
        kd="${DATASET_DIR}/${s}_knead_out"
        R1=$(ls "${kd}"/*_kneaddata_paired_1.fastq 2>/dev/null | head -n 1)
        R2=$(ls "${kd}"/*_kneaddata_paired_2.fastq 2>/dev/null | head -n 1)
        [[ -f "$R1" && -f "$R2" ]] && \
            R1_LIST="${R1_LIST},${R1}" && R2_LIST="${R2_LIST},${R2}" && \
            echo "  Added: $s"
    done < <(get_samples "$DATASET" "Metagenome")

    R1_LIST="${R1_LIST#,}"; R2_LIST="${R2_LIST#,}"
    [[ -z "$R1_LIST" ]] && echo "ERROR: No DNA reads for assembly" && exit 1

    rm -rf "$ASSEMBLY_DIR"
    megahit -1 "$R1_LIST" -2 "$R2_LIST" \
        -o "$ASSEMBLY_DIR" --min-contig-len 500 -t 48
fi

# STEP 3 — Prodigal annotation

echo " STEP 3: Prodigal "
mkdir -p "$ANNOTATION_DIR"

if [ ! -f "${ANNOTATION_DIR}/prodigal.gff" ]; then
    prodigal \
        -i "${ASSEMBLY_DIR}/final.contigs.fa" \
        -f gff -o "${ANNOTATION_DIR}/prodigal.gff" \
        -a "${ANNOTATION_DIR}/prodigal.faa" -p meta
fi

# STEP 4 — Bowtie2 index

echo " STEP 4: Bowtie2 index "
INDEX_PREFIX="${ASSEMBLY_DIR}/bowtie2_index"

if [ ! -f "${INDEX_PREFIX}.1.bt2" ]; then
    bowtie2-build --threads 48 \
        "${ASSEMBLY_DIR}/final.contigs.fa" "$INDEX_PREFIX"
fi

echo " Complete: $DATASET "
echo "Next -> sbatch 04_snapt.sh $DATASET"
