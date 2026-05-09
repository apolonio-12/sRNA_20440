#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH -n 36
#SBATCH -J humann
#SBATCH -p mit_preemptable
#SBATCH --mem=250GB
#SBATCH --output=/home/apolonio/orcd/scratch/20_440/project_logs/humannAD.log

set -euo pipefail

source /home/apolonio/miniconda3/etc/profile.d/conda.sh
conda deactivate 
conda activate /home/apolonio/.conda/envs/humann_env

# ── Configuration ─────────────────────────────────────────────────────────────
METADATA="/home/apolonio/orcd/scratch/20_440/Metadata.csv"

MPA_DB="/home/apolonio/.conda/envs/humann_env/lib/python3.12/site-packages/metaphlan/metaphlan_databases"
MPA_DB_INDEX="mpa_vOct22_CHOCOPhlAnSGB_202403"

# Dataset flag — override at submission:
#   sbatch 04b_humann.sh AD
#   sbatch 04b_humann.sh CB
#   sbatch 04b_humann.sh AD Metagenome         (DNA only)
#   sbatch 04b_humann.sh AD Metatranscriptome  (RNA only)
#   sbatch 04b_humann.sh CB Both               (default)
DATASET="${1:-AD}"
TYPE_FILTER="${2:-Both}"    # Metagenome | Metatranscriptome | Both

BASE_DIR="/home/apolonio/orcd/scratch/20_440/${DATASET}"

echo "===== HUMAnN: dataset=$DATASET  type=$TYPE_FILTER ====="
echo "  Base: $BASE_DIR"

# ── Helper: get sample Names filtered by dataset + type ───────────────────────
get_samples() {
    local ds="$1"
    local tp="$2"   # Metagenome or Metatranscriptome
    tail -n +2 "$METADATA" \
        | awk -F'\t' -v ds="$ds" -v tp="$tp" \
            '$3==tp && $9==ds { gsub(/\r/,"",$4); print $4 }'
}

# ── Core HUMAnN function ──────────────────────────────────────────────────────
run_humann() {
    local sample="$1"
    local data_type="$2"    # Metagenome or Metatranscriptome

    local knead_dir="${BASE_DIR}/${sample}_knead_out"
    local out_dir="${BASE_DIR}/${sample}_humann_output"

    # Skip if already done
    if [ -d "$out_dir" ]; then
        echo "  EXISTS: $sample — skipping"
        return
    fi

    # Locate kneaddata paired outputs
    local R1
    local R2
    # Detect paired files
    R1=$(find "$knead_dir" -maxdepth 1 -type f -name "*_kneaddata_paired_1.fastq" | head -n 1)
    R2=$(find "$knead_dir" -maxdepth 1 -type f -name "*_kneaddata_paired_2.fastq" | head -n 1)

    # Detect single-end file (final cleaned output only)
    SINGLE=$(find "$knead_dir" -maxdepth 1 -type f \
        -name "*_kneaddata.fastq" | head -n 1)

    unset INPUT_FILE
    unset MODE

    if [[ -f "$R1" && -f "$R2" ]]; then
        MODE="paired"
        INPUT_FILE="${knead_dir}/${sample}_paired_combined.fastq"

        echo "  Concatenating paired reads: $sample ($data_type)"
        cat "$R1" "$R2" > "$INPUT_FILE"

    elif [[ -f "$SINGLE" ]]; then
        MODE="single"
        INPUT_FILE="$SINGLE"

        echo "  Using single-end reads: $sample ($data_type)"

    else
        echo "  WARNING: No valid kneaddata outputs for $sample — skipping"
        return
    fi

    mkdir -p "$out_dir"
    echo "  Running HUMAnN: $sample"

    humann \
        --input "$INPUT_FILE" \
        --output "$out_dir" \
        --threads 24 \
        --remove-temp-output \
        --metaphlan-options "--index $MPA_DB_INDEX --offline -t rel_ab_w_read_stats"

    # Delete combined to save space
    [[ "$MODE" == "paired" ]] && rm -f "$INPUT_FILE"
    echo "  Done: $sample → $out_dir"
}

# ═══════════════════════════════════════════════════════════════════════════════
# Run DNA samples
# ═══════════════════════════════════════════════════════════════════════════════
if [[ "$TYPE_FILTER" == "Metagenome" || "$TYPE_FILTER" == "Both" ]]; then
    echo ""
    echo "===== DNA (Metagenome) samples ====="
    while IFS= read -r sample; do
        run_humann "$sample" "Metagenome"
    done < <(get_samples "$DATASET" "Metagenome")
fi

# ═══════════════════════════════════════════════════════════════════════════════
# Run RNA samples
# ═══════════════════════════════════════════════════════════════════════════════
if [[ "$TYPE_FILTER" == "Metatranscriptome" || "$TYPE_FILTER" == "Both" ]]; then
    echo ""
    echo "===== RNA (Metatranscriptome) samples ====="
    while IFS= read -r sample; do
        run_humann "$sample" "Metatranscriptome"
    done < <(get_samples "$DATASET" "Metatranscriptome")
fi

# ═══════════════════════════════════════════════════════════════════════════════
# Join tables (mirrors your DE_Bay join commands exactly)
# ═══════════════════════════════════════════════════════════════════════════════
echo ""
echo "===== Joining HUMAnN tables ====="

# Collect DNA and RNA output dirs separately
DNA_OUT="${BASE_DIR}/Total_DNA_humann_output"
RNA_OUT="${BASE_DIR}/Total_RNA_humann_output"
mkdir -p "$DNA_OUT" "$RNA_OUT"

# Copy DNA outputs
if [[ "$TYPE_FILTER" == "Metagenome" || "$TYPE_FILTER" == "Both" ]]; then
    while IFS= read -r sample; do
        src="${BASE_DIR}/${sample}_humann_output"
        [ -d "$src" ] || continue
        find "$src" -name "*.tsv" -exec cp {} "$DNA_OUT/" \;
    done < <(get_samples "$DATASET" "Metagenome")
fi

# Copy RNA outputs
if [[ "$TYPE_FILTER" == "Metatranscriptome" || "$TYPE_FILTER" == "Both" ]]; then
    while IFS= read -r sample; do
        src="${BASE_DIR}/${sample}_humann_output"
        [ -d "$src" ] || continue
        find "$src" -name "*.tsv" -exec cp {} "$RNA_OUT/" \;
    done < <(get_samples "$DATASET" "Metatranscriptome")
fi

# ── Join DNA tables ───────────────────────────────────────────────────────────
if [ "$(ls -A $DNA_OUT 2>/dev/null)" ]; then
    echo "  Joining DNA tables..."
    humann_join_tables --input "$DNA_OUT" \
        --output "${BASE_DIR}/Total_DNA_genefamilies_humann_table.tsv" \
        --file_name gene

    humann_join_tables --input "$DNA_OUT" \
        --output "${BASE_DIR}/Total_DNA_pathabundance_humann_table.tsv" \
        --file_name pathabundance

    humann_join_tables --input "$DNA_OUT" \
        --output "${BASE_DIR}/Total_DNA_reactions_humann_table.tsv" \
        --file_name reactions
fi

# ── Join RNA tables ───────────────────────────────────────────────────────────
if [ "$(ls -A $RNA_OUT 2>/dev/null)" ]; then
    echo "  Joining RNA tables..."
    humann_join_tables --input "$RNA_OUT" \
        --output "${BASE_DIR}/Total_RNA_genefamilies_humann_table.tsv" \
        --file_name gene

    humann_join_tables --input "$RNA_OUT" \
        --output "${BASE_DIR}/Total_RNA_pathabundance_humann_table.tsv" \
        --file_name pathabundance

    humann_join_tables --input "$RNA_OUT" \
        --output "${BASE_DIR}/Total_RNA_reactions_humann_table.tsv" \
        --file_name reactions

    humann_join_tables --input "$RNA_OUT" \
        --output "${BASE_DIR}/Total_RNA_profile_humann_table.tsv" \
        --file_name profile
fi

echo ""
echo "===== HUMAnN complete: $DATASET ====="
echo "Key outputs:"
echo "  ${BASE_DIR}/Total_DNA_pathabundance_humann_table.tsv"
echo "  ${BASE_DIR}/Total_RNA_pathabundance_humann_table.tsv"
echo "Next: sbatch 07_feature_matrix.sh $DATASET"
