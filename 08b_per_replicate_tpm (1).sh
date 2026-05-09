#!/bin/bash
# 08b_per_replicate_tpm.sh
# Quantifies each replicate BAM individually with StringTie
# Output enables per-replicate Pearson correlation (Section 5 full version)
#
# Usage:
#   sbatch 08b_per_replicate_tpm.sh AD 2016
#   sbatch 08b_per_replicate_tpm.sh AD 2017

#SBATCH -t 4:00:00
#SBATCH -n 16
#SBATCH -J per_rep_tpm
#SBATCH -p mit_normal
#SBATCH --mem=64GB
#SBATCH --output=/home/apolonio/orcd/scratch/20_440/project_logs/per_rep_tpm_%1_%2.log

source /home/apolonio/miniconda3/etc/profile.d/conda.sh
conda activate snapt_env
set -euo pipefail

DATASET="${1:-AD}"
TIMEGROUP="${2:-2016}"
BASE_DIR="/home/apolonio/orcd/scratch/20_440/${DATASET}"
BAM_DIR="${BASE_DIR}/bam_files"
SNAPT_OUT="${BASE_DIR}/snapt_${TIMEGROUP}"
TPM_DIR="${BASE_DIR}/per_replicate_tpm_${TIMEGROUP}"
mkdir -p "$TPM_DIR"

ANNOTATION="${SNAPT_OUT}/small_ncRNAs.gff"
[ -f "$ANNOTATION" ] || { echo "ERROR: $ANNOTATION not found"; exit 1; }

if [ "$DATASET" = "AD" ]; then
    BAM_PATTERN="${TIMEGROUP}"
else
    BAM_PATTERN="${TIMEGROUP}"
fi

echo "===== Per-replicate TPM: $DATASET $TIMEGROUP ====="
echo "  BAM pattern: *${BAM_PATTERN}*"

for bam in "${BAM_DIR}"/*${BAM_PATTERN}*.sorted.bam; do
    sample=$(basename "$bam" .sorted.bam)
    outdir="${TPM_DIR}/${sample}"
    mkdir -p "$outdir"

    if [ -f "${outdir}/${sample}_gene_abund.tsv" ]; then
        echo "  EXISTS: $sample — skipping"
        continue
    fi

    echo "  Quantifying: $sample"
    stringtie "$bam" \
        -G "$ANNOTATION" \
        -o "${outdir}/${sample}.gtf" \
        -e \
        -A "${outdir}/${sample}_gene_abund.tsv" \
        -p 16
    echo "  Done: $sample"
done

# ── Build TPM matrix across replicates ────────────────────────────────────────
echo ""
echo "Building TPM matrix..."

export TPM_DIR
python3 - <<'PYEOF'
import os
import pandas as pd

tpm_dir = os.environ["TPM_DIR"]
rows = []

for sample_dir in sorted(os.listdir(tpm_dir)):
    full_path = os.path.join(tpm_dir, sample_dir)
    if not os.path.isdir(full_path):
        continue
    abund = os.path.join(full_path, f"{sample_dir}_gene_abund.tsv")
    if not os.path.exists(abund):
        continue
    df = pd.read_csv(abund, sep="\t")
    df.columns = [c.strip() for c in df.columns]
    df["sample"] = sample_dir
    rows.append(df)

if not rows:
    print("ERROR: No abundance files found")
    exit(1)

all_abund = pd.concat(rows, ignore_index=True)
print(f"Samples: {all_abund['sample'].nunique()}")
print(f"Genes:   {all_abund['Gene ID'].nunique()}")

tpm_matrix = all_abund.pivot_table(
    index="Gene ID", columns="sample",
    values="TPM", aggfunc="first"
)
cov_matrix = all_abund.pivot_table(
    index="Gene ID", columns="sample",
    values="Coverage", aggfunc="first"
)

out_tpm = os.path.join(tpm_dir, "tpm_matrix.tsv")
out_cov = os.path.join(tpm_dir, "coverage_matrix.tsv")
tpm_matrix.to_csv(out_tpm, sep="\t")
cov_matrix.to_csv(out_cov, sep="\t")

print(f"Saved: {out_tpm}")
print(f"Saved: {out_cov}")
print(f"Matrix shape: {tpm_matrix.shape}  (genes × samples)")
PYEOF

echo ""
echo "===== Done: $DATASET $TIMEGROUP ====="
echo "  TPM matrix: ${TPM_DIR}/tpm_matrix.tsv"
echo "  Add --tpm_matrix to 08_srna_analysis.py when ready for full correlation"
