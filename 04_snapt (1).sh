#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH -n 48
#SBATCH -J blast
#SBATCH -p mit_normal
#SBATCH --mem=250GB
#SBATCH --output=/home/apolonio/orcd/scratch/20_440/project_logs/snapt_%j.log

set -euo pipefail

#module load miniforge
#conda activate snapt_env

METADATA="/home/apolonio/orcd/scratch/20_440/Metadata.csv"
BASE_DIR="/home/apolonio/orcd/scratch/20_440"
NR_DB="/home/apolonio/orcd/scratch/databases/nr/nr.dmnd"
RFAM_DB="/home/apolonio/orcd/scratch/databases/rfam/Rfam.cm"

# Dataset and timepoint group flags:
#   AD samples: sbatch 04_snapt.sh AD 2016
#               sbatch 04_snapt.sh AD 2017
#   CB samples: sbatch 04_snapt.sh CB Fall
#               sbatch 04_snapt.sh CB Sum
DATASET="${1:-AD}"
TIMEGROUP="${2:-2016}"

DATASET_DIR="${BASE_DIR}/${DATASET}"
ASSEMBLY="${DATASET_DIR}/coassembly/final.contigs.fa"
ANNOTATION="${DATASET_DIR}/annotation/prodigal.gff"
SNAPT_OUT="${DATASET_DIR}/snapt_${TIMEGROUP}"
READS_DIR="${DATASET_DIR}/combined_reads"

mkdir -p "$SNAPT_OUT" "$READS_DIR"

echo "===== SnapT: dataset=$DATASET  timegroup=$TIMEGROUP ====="

# ── Helper: does sample name belong to this timegroup? ────────────────────────
matches_timegroup() {
    local name="$1"
    local tg="$2"
    local ds="$3"

    if [[ "$ds" == "AD" ]]; then
        # Match by year in sample name: 2016 or 2017
        echo "$name" | grep -q "_${tg}_" && return 0 || return 1
    else
        # CB: match by season keyword Fall or Sum
        echo "$name" | grep -q "_${tg}_" && return 0 || return 1
    fi
}

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 0 — Sanitize FASTA headers (strip everything after first whitespace)
# ═══════════════════════════════════════════════════════════════════════════════
ASSEMBLY_CLEAN="${DATASET_DIR}/coassembly/final.contigs.clean.fa"

if [[ -f "$ASSEMBLY_CLEAN" ]]; then
    echo "  Clean assembly exists — skipping header sanitization"
else
    echo "===== STEP 0: Sanitizing FASTA headers ====="
    
    # Check if headers actually need cleaning
    if grep -m1 "^>" "$ASSEMBLY" | grep -q " "; then
        echo "  Whitespace found in headers — stripping to bare IDs"
        sed 's/>\(\S*\).*/>\1/' "$ASSEMBLY" > "$ASSEMBLY_CLEAN"
        
        # Verify contig counts match
        ORIG=$(grep -c "^>" "$ASSEMBLY")
        CLEAN=$(grep -c "^>" "$ASSEMBLY_CLEAN")
        if [[ "$ORIG" != "$CLEAN" ]]; then
            echo "ERROR: Contig count mismatch after cleaning ($ORIG vs $CLEAN) — aborting"
            exit 1
        fi
        echo "  Verified: $CLEAN contigs in cleaned assembly"
    else
        echo "  Headers already clean — symlinking"
        ln -sf "$ASSEMBLY" "$ASSEMBLY_CLEAN"
    fi
fi

# Use the clean assembly from here on
ASSEMBLY="$ASSEMBLY_CLEAN"

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 1 — Concatenate RNA reads for this timegroup only
# ═══════════════════════════════════════════════════════════════════════════════
echo "===== STEP 1: Concatenating reads for $TIMEGROUP ====="

COMBINED_R1="${READS_DIR}/${DATASET}_${TIMEGROUP}_R1.fastq"
COMBINED_R2="${READS_DIR}/${DATASET}_${TIMEGROUP}_R2.fastq"

if [[ -f "$COMBINED_R1" && -f "$COMBINED_R2" ]]; then
    echo "  Combined reads exist — skipping"
else
    > "$COMBINED_R1"
    > "$COMBINED_R2"

    while IFS=$'\t' read -r bioproject srr type name title env lat lon dataset; do
        dataset=$(echo "$dataset" | tr -d '\r')
        type=$(echo "$type"       | tr -d '\r')
        name=$(echo "$name"       | tr -d '\r')

        [[ "$dataset" != "$DATASET" ]]        && continue
        [[ "$type"    != "Metatranscriptome" ]] && continue

        # Only include samples matching this timegroup
        matches_timegroup "$name" "$TIMEGROUP" "$DATASET" || continue

        knead_dir="${DATASET_DIR}/${name}_knead_out"
        R1=$(ls "${knead_dir}"/*_kneaddata_paired_1.fastq 2>/dev/null | head -n 1)
        R2=$(ls "${knead_dir}"/*_kneaddata_paired_2.fastq 2>/dev/null | head -n 1)

        if [[ -f "$R1" && -f "$R2" ]]; then
            echo "  Adding: $name"
            cat "$R1" >> "$COMBINED_R1"
            cat "$R2" >> "$COMBINED_R2"
        else
            echo "  WARNING: No paired reads for $name — skipping"
        fi

    done < <(tail -n +2 "$METADATA")

    LINES=$(wc -l < "$COMBINED_R1")
    echo "  Combined R1: $LINES lines (~$((LINES/4)) reads)"
fi

# Check we actually got reads
if [[ ! -s "$COMBINED_R1" ]]; then
    echo "ERROR: No reads collected for $DATASET $TIMEGROUP — check sample names in metadata"
    exit 1
fi

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 2 — Run SnapT for this timegroup
# ═══════════════════════════════════════════════════════════════════════════════
echo "===== STEP 2: SnapT ($TIMEGROUP) ====="

if [ -f "${SNAPT_OUT}/sRNA_annotations.tsv" ]; then
    echo "  SnapT output exists for $TIMEGROUP — skipping"
else
    snapt \
        -1 "$COMBINED_R1" \
        -2 "$COMBINED_R2" \
        -g "$ASSEMBLY" \
        -a "$ANNOTATION" \
        -D "$NR_DB" \
        -R "$RFAM_DB" \
        -o "$SNAPT_OUT" \
        -t 36 \
        -r RF

    echo "  SnapT complete for $TIMEGROUP"
fi

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 3 — Sanity check
# ═══════════════════════════════════════════════════════════════════════════════
echo "===== STEP 3: Output check ====="
ls -lh "$SNAPT_OUT"

if [ -f "${SNAPT_OUT}/sRNA_annotations.tsv" ]; then
    COUNT=$(tail -n +2 "${SNAPT_OUT}/sRNA_annotations.tsv" | wc -l)
    echo "sRNAs detected ($TIMEGROUP): $COUNT"
    head -n 3 "${SNAPT_OUT}/sRNA_annotations.tsv"
else
    echo "WARNING: sRNA_annotations.tsv not found — check SnapT logs in $SNAPT_OUT"
fi

echo "===== SnapT complete: $DATASET $TIMEGROUP ====="
echo ""
echo "Submit the other timegroup if not already running:"
if [[ "$DATASET" == "AD" ]]; then
    OTHER="2017"; [[ "$TIMEGROUP" == "2017" ]] && OTHER="2016"
    echo "  sbatch 04_snapt.sh AD $OTHER"
else
    OTHER="Sum"; [[ "$TIMEGROUP" == "Sum" ]] && OTHER="Fall"
    echo "  sbatch 04_snapt.sh CB $OTHER"
fi
echo ""
echo "Once both timegroups done:"
echo "  sbatch 05_rnafold.sh $DATASET 2016"
echo "  sbatch 05_rnafold.sh $DATASET 2017"
