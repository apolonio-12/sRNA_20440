#!/bin/bash
#SBATCH -t 4:00:00
#SBATCH -n 4
#SBATCH -J feature_matrix
#SBATCH -p mit_normal
#SBATCH --mem=32GB
#SBATCH --output=/home/apolonio/orcd/scratch/20_440/project_logs/feature_matrix_%j.log

set -euo pipefail

#module load miniforge
#conda activate snapt_env
source /home/apolonio/miniconda3/etc/profile.d/conda.sh
conda activate snapt_env
# Usage:
#   sbatch 07_feature_matrix_final.sh AD
#   sbatch 07_feature_matrix_final.sh CB
DATASET="${1:-AD}"

BASE="/home/apolonio/orcd/scratch/20_440/${DATASET}"
SCRIPT="/home/apolonio/orcd/scratch/20_440/07_build_feature_matrix.py"
METADATA="/home/apolonio/orcd/scratch/20_440/Metadata.csv"

echo "Building feature matrix for: $DATASET"

# Check if IntaRNA output exists for first timegroup — include if available
if [[ "$DATASET" == "AD" ]]; then
    FIRST_TG="2016"
else
    FIRST_TG="Fall"
fi

INTARNA_ARG=""
INTARNA_FILE="${BASE}/intarna_${FIRST_TG}/intarna_summary.tsv"
if [ -f "$INTARNA_FILE" ]; then
    echo "  IntaRNA summary found — including interaction features"
    INTARNA_ARG="--intarna $INTARNA_FILE"
else
    echo "  No IntaRNA summary found — skipping interaction features"
    echo "  (Run 06_intarna_final.sh and rerun to add these features later)"
fi

python3 "$SCRIPT" \
    --dataset    "$DATASET" \
    --metadata   "$METADATA" \
    --snapt      "${BASE}/snapt_${FIRST_TG}/small_antisense_ncRNAs.gff" \
    --rnafold    "${BASE}/rnafold_${FIRST_TG}/rnafold_summary.tsv" \
    --humann_rna "${BASE}/Total_RNA_pathabundance_humann_table.tsv" \
    --humann_dna "${BASE}/Total_DNA_pathabundance_humann_table.tsv" \
    --metaphlan  "${BASE}/Total_DNA_merged_abundance_table.tsv" \
    --output     "${BASE}/feature_matrix.tsv" \
    $INTARNA_ARG

echo "===== Done: $DATASET ====="
echo "Outputs:"
echo "  ${BASE}/feature_matrix.tsv"
echo "  ${BASE}/feature_matrix_species_level.tsv"
echo "  ${BASE}/feature_matrix_summary.json"
