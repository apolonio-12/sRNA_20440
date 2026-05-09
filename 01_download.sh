#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH -n 30
#SBATCH -J download_reads
#SBATCH -p mit_normal
#SBATCH --mem=200GB
#SBATCH --output=/home/apolonio/orcd/scratch/20_440/project_logs/SRA_download2.log

set -euo pipefail

export PATH=$PATH:/home/apolonio/sratoolkit.3.2.1-alma_linux64/bin

# ── Configuration ─────────────────────────────────────────────────────────────
METADATA="/home/apolonio/orcd/scratch/20_440/Metadata.csv"
BASE_DIR="/home/apolonio/orcd/scratch/20_440"

# Filter flags — override at submission time:
#   sbatch 02_download.sh AD ALL
#   sbatch 02_download.sh CB Metatranscriptome
#   sbatch 02_download.sh ALL ALL
DATASET="${1:-ALL}"     # AD | CB | ALL
DATA_TYPE="${2:-ALL}"   # Metagenome | Metatranscriptome | ALL

echo "===== Download parameters ====="
echo "  Dataset filter:   $DATASET"
echo "  Data type filter: $DATA_TYPE"
echo "  Metadata file:    $METADATA"
echo "  Base directory:   $BASE_DIR"
echo "==============================="

mkdir -p "${BASE_DIR}/sra_tmp"

# ── Download function ─────────────────────────────────────────────────────────
download_run() {
    local srr="$1"
    local name="$2"
    local dataset="$3"

    local out_base="${BASE_DIR}/${dataset}"
    mkdir -p "$out_base"
    local sample_dir="${out_base}/${name}"

    if [ -d "$sample_dir" ]; then
        echo "  EXISTS: ${dataset}/${name} — skipping"
        return
    fi

    echo "  Downloading $srr → ${dataset}/${name}"
    mkdir -p "$sample_dir"

    prefetch "$srr" --output-directory "${BASE_DIR}/sra_tmp"

    fasterq-dump \
        --threads 20 \
        --split-files \
        --outdir "$sample_dir" \
        "${BASE_DIR}/sra_tmp/${srr}/${srr}.sra"

    # Rename SRR accession → sample name
    for suffix in "_1.fastq" "_2.fastq" ".fastq"; do
        [ -f "${sample_dir}/${srr}${suffix}" ] && \
            mv "${sample_dir}/${srr}${suffix}" \
               "${sample_dir}/${name}${suffix}"
    done

    # Compress
    #pigz -p 8 "${sample_dir}"/*.fastq 2>/dev/null || true

    # Clean prefetch cache
    rm -rf "${BASE_DIR}/sra_tmp/${srr}"

    echo "  Done: ${dataset}/${name}"
}

# ── Read metadata and loop ────────────────────────────────────────────────────
echo ""
echo "===== Starting downloads ====="

tail -n +2 "$METADATA" | while IFS=$'\t' read -r bioproject srr type name title env lat lon dataset; do

    # Strip Windows carriage returns
    srr=$(echo "$srr"         | tr -d '\r')
    type=$(echo "$type"       | tr -d '\r')
    name=$(echo "$name"       | tr -d '\r')
    dataset=$(echo "$dataset" | tr -d '\r')

    # Apply dataset filter
    if [[ "$DATASET" != "ALL" && "$dataset" != "$DATASET" ]]; then
        continue
    fi

    # Apply type filter
    if [[ "$DATA_TYPE" == "Metagenome" && "$type" != "Metagenome" ]]; then
        continue
    fi
    if [[ "$DATA_TYPE" == "Metatranscriptome" && "$type" != "Metatranscriptome" ]]; then
        continue
    fi

    download_run "$srr" "$name" "$dataset"

done

echo ""
echo "===== Download complete ====="
echo ""
echo "Submission examples for future reference:"
echo "  sbatch 02_download.sh AD  ALL               # all Atacama"
echo "  sbatch 02_download.sh CB  ALL               # all Chesapeake/Delaware"
echo "  sbatch 02_download.sh AD  Metagenome        # Atacama DNA only"
echo "  sbatch 02_download.sh CB  Metatranscriptome # CB RNA only"
echo "  sbatch 02_download.sh ALL ALL               # everything"
