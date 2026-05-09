#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH -n 36
#SBATCH -J rnafold
#SBATCH -p mit_normal
#SBATCH --mem=200GB
#SBATCH --output=/home/apolonio/orcd/scratch/20_440/project_logs/rnafold%j.log

set -euo pipefail

#module load miniforge
#conda activate snapt_env

# Usage:
#   sbatch 05_rnafold.sh AD 2016
#   sbatch 05_rnafold.sh AD 2017
#   sbatch 05_rnafold.sh CB Fall
#   sbatch 05_rnafold.sh CB Sum
#sbatch --output=/home/apolonio/orcd/scratch/20_440/project_logs/rnafold_AD_2017.log     05_rnafold.sh AD 2016
DATASET="${1:-AD}"
TIMEGROUP="${2:-2016}"

BASE_DIR="/home/apolonio/orcd/scratch/20_440/${DATASET}"
SNAPT_OUT="${BASE_DIR}/snapt_${TIMEGROUP}"
RNAFOLD_DIR="${BASE_DIR}/rnafold_${TIMEGROUP}"
THREADS=36

mkdir -p "$RNAFOLD_DIR"

echo "===== RNAfold: dataset=$DATASET  timegroup=$TIMEGROUP ====="

# Find sRNA FASTA from SnapT output
SRNA_FA="${BASE_DIR}/snapt_${TIMEGROUP}/all_sRNAs_${TIMEGROUP}.fa"

if [ ! -f "$SRNA_FA" ]; then
    echo "Combining asRNA and itsRNA FASTA files..."
    cat "${SNAPT_OUT}/small_antisense_ncRNAs.fa" \
        "${SNAPT_OUT}/small_intergenic_ncRNAs.fa" \
        > "$SRNA_FA"
    echo "Combined: $(grep -c "^>" $SRNA_FA) total sRNAs"
fi

echo "sRNA FASTA: $SRNA_FA"
TOTAL=$(grep -c "^>" "$SRNA_FA")
echo "Total sRNAs: $TOTAL"

SPLIT_DIR="${RNAFOLD_DIR}/individual_seqs"
RESULTS_DIR="${RNAFOLD_DIR}/results"
mkdir -p "$SPLIT_DIR" "$RESULTS_DIR"

# ── Split FASTA into one file per sRNA ───────────────────────────────────────
if [ "$(ls -A $SPLIT_DIR 2>/dev/null | wc -l)" -eq 0 ]; then
    echo "Splitting FASTA..."
    awk '/^>/ {
        if (seqfile) close(seqfile)
        id = substr($0,2); gsub(/[^A-Za-z0-9_-]/,"_",id)
        seqfile = "'"$SPLIT_DIR"'/" id ".fa"
    } { print > seqfile }' "$SRNA_FA"
    echo "Split into $(ls $SPLIT_DIR | wc -l) files"
else
    echo "Split files already exist — skipping"
fi

# ── Parallel RNAfold ──────────────────────────────────────────────────────────
echo "Running RNAfold..."

fold_one() {
    local fa="$1"
    local out="${2}/$(basename $fa .fa).txt"
    [ -f "$out" ] && return
    RNAfold --noPS --partfunc < "$fa" > "$out" 2>/dev/null
}
export -f fold_one

ls "${SPLIT_DIR}"/*.fa | parallel -j "$THREADS" fold_one {} "$RESULTS_DIR"

DONE=$(ls "$RESULTS_DIR" | wc -l)
echo "RNAfold complete: $DONE / $TOTAL structures"

# ── Parse results into summary TSV ───────────────────────────────────────────
SUMMARY="${RNAFOLD_DIR}/rnafold_summary.tsv"
echo "$RESULTS_DIR" > /tmp/rnafold_results_dir.txt
echo "$SUMMARY"     > /tmp/rnafold_summary_path.txt

python3 - <<'PYEOF'
import os, re

results_dir = open("/tmp/rnafold_results_dir.txt").read().strip()
out_path    = open("/tmp/rnafold_summary_path.txt").read().strip()

rows = []
for fname in sorted(os.listdir(results_dir)):
    if not fname.endswith(".txt"):
        continue
    srna_id = fname.replace(".txt", "")
    lines = open(os.path.join(results_dir, fname)).readlines()
    if len(lines) < 3:
        continue
    seq      = lines[1].strip()
    mfe_line = lines[2].strip()
    mfe_m    = re.search(r"\((-?\d+\.\d+)\)", mfe_line)
    mfe      = float(mfe_m.group(1)) if mfe_m else float("nan")
    structure = mfe_line.split(" ")[0]
    ens = float("nan")
    for l in lines[3:]:
        m = re.search(r"ensemble energy:\s*(-?\d+\.\d+)", l)
        if m:
            ens = float(m.group(1))
            break
    gc = (seq.count("G") + seq.count("C")) / len(seq) if seq else float("nan")
    se = structure.count(".") / len(structure) if structure else float("nan")
    rows.append("\t".join([
        srna_id, str(len(seq)),
        f"{gc:.4f}", f"{mfe:.2f}", f"{ens:.2f}", f"{se:.4f}",
        structure
    ]))

with open(out_path, "w") as f:
    f.write("srna_id\tlength\tgc_content\tmfe\tensemble_energy\tstruct_entropy\tstructure\n")
    f.write("\n".join(rows) + "\n")

print(f"Written {len(rows)} rows to {out_path}")
PYEOF

echo "RNAfold summary: $SUMMARY"
echo "Next: sbatch 06_intarna.sh $DATASET $TIMEGROUP"
