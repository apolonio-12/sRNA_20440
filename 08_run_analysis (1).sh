#!/bin/bash
# 08_run_analysis.sh
#
# Usage:
#   sbatch 08_run_analysis.sh AD 2016 2017
#   sbatch 08_run_analysis.sh CB Fall Sum

#SBATCH -t 12:00:00
#SBATCH -n 12
#SBATCH -J srna_analysis
#SBATCH -p mit_preemptable
#SBATCH --mem=250GB
#SBATCH --output=/home/apolonio/orcd/scratch/20_440/project_logs/srna_analysis_%1_%2_%3.log

source /home/apolonio/miniconda3/etc/profile.d/conda.sh
conda activate snapt_env


set -euo pipefail

DATASET="${1:-AD}"
T1="${2:-2016}"
T2="${3:-2017}"
BASE="/home/apolonio/orcd/scratch/20_440/${DATASET}"
SCRIPTS="/home/apolonio/orcd/scratch/20_440"

echo "===== sRNA Analysis: $DATASET  $T1 → $T2 ====="

# ── Verify inputs ─────────────────────────────────────────────────────────────
check() { [ -f "$1" ] || { echo "ERROR: Missing $1"; exit 1; }; }

check "${BASE}/snapt_${T1}/small_antisense_ncRNAs.gff"
check "${BASE}/snapt_${T1}/small_intergenic_ncRNAs.gff"
check "${BASE}/snapt_${T2}/small_antisense_ncRNAs.gff"
check "${BASE}/snapt_${T2}/small_intergenic_ncRNAs.gff"
check "${BASE}/rnafold_${T1}/rnafold_summary.tsv"
check "${BASE}/rnafold_${T2}/rnafold_summary.tsv"
check "${BASE}/feature_matrix.tsv"

# Confirm fixed IntaRNA summaries exist — run fix_intarna_summary.py first if not
for tg in "$T1" "$T2"; do
    fixed="${BASE}/intarna_${tg}/intarna_summary_fixed.tsv"
    if [ ! -f "$fixed" ]; then
        echo "ERROR: $fixed not found."
        echo "Run: python fix_intarna_summary.py --dataset $DATASET --timegroups $T1 $T2"
        exit 1
    fi
done

echo "All inputs verified."

OUTDIR="${BASE}/analysis_output_${T1}_${T2}"
mkdir -p "$OUTDIR"

python3 "${SCRIPTS}/08_srna_analysis.py" \
    --dataset        "$DATASET" \
    --snapt_t1       "${BASE}/snapt_${T1}" \
    --snapt_t2       "${BASE}/snapt_${T2}" \
    --rnafold_t1     "${BASE}/rnafold_${T1}/rnafold_summary.tsv" \
    --rnafold_t2     "${BASE}/rnafold_${T2}/rnafold_summary.tsv" \
    --feature_matrix "${BASE}/feature_matrix.tsv" \
    --outdir         "$OUTDIR" \
    --intarna_t1     "${BASE}/intarna_${T1}/intarna_summary_fixed.tsv" \
    --intarna_t2     "${BASE}/intarna_${T2}/intarna_summary_fixed.tsv" \
    --label_t1       "$T1" \
    --label_t2       "$T2" \
    --tpm_matrix_t1  "${BASE}/per_replicate_tpm_${T1}/tpm_matrix.tsv" \
    --tpm_matrix_t2  "${BASE}/per_replicate_tpm_${T2}/tpm_matrix.tsv"

echo ""
echo "===== Done: $OUTDIR ====="
ls -lh "$OUTDIR"
