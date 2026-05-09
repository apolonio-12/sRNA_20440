#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH -n 48
#SBATCH -J intarna
#SBATCH -p mit_normal
#SBATCH --mem=200GB
#SBATCH --output=/home/apolonio/orcd/scratch/20_440/project_logs/intaRNA_%1_%2.log

# Usage:
#   sbatch 06_intarna.sh AD 2016
#   sbatch 06_intarna.sh AD 2017
#   sbatch 06_intarna.sh CB Fall
#   sbatch 06_intarna.sh CB Sum

set -euo pipefail

DATASET="${1:-AD}"
TIMEGROUP="${2:-2016}"
PVAL_THRESH="${3:-0.05}"

BASE_DIR="/home/apolonio/orcd/scratch/20_440/${DATASET}"
SNAPT_OUT="${BASE_DIR}/snapt_${TIMEGROUP}"
ANNOTATION_DIR="${BASE_DIR}/annotation"
ASSEMBLY="${BASE_DIR}/coassembly/final.contigs.fa"
INTARNA_DIR="${BASE_DIR}/intarna_${TIMEGROUP}"
BAM_DIR="${BASE_DIR}/bam_files"
mkdir -p "$INTARNA_DIR"

echo "===== IntaRNA: dataset=$DATASET  timegroup=$TIMEGROUP  pval<${PVAL_THRESH} ====="

# ── Find sRNA FASTA ───────────────────────────────────────────────────────────
SRNA_FA="${SNAPT_OUT}/all_sRNAs_${TIMEGROUP}.fa"
if [ ! -f "$SRNA_FA" ]; then
    echo "Combining asRNA and itsRNA FASTA files..."
    cat "${SNAPT_OUT}/small_antisense_ncRNAs.fa" \
        "${SNAPT_OUT}/small_intergenic_ncRNAs.fa" \
        > "$SRNA_FA"
    echo "Combined: $(grep -c '^>' $SRNA_FA) total sRNAs"
fi
echo "sRNA FASTA: $SRNA_FA  ($(grep -c '^>' $SRNA_FA) seqs)"

# ── STEP 1: Get expressed contigs from SnapT output ──────────────────────────
EXPRESSED_CONTIGS="${INTARNA_DIR}/expressed_contigs.txt"
cut -f1 "${SNAPT_OUT}/small_ncRNAs.gff" | sort -u > "$EXPRESSED_CONTIGS"
echo "SnapT expressed contigs: $(wc -l < $EXPRESSED_CONTIGS)"

# Optionally expand with BAM coverage (AD: 2016/2017_T2, CB: Fall/Sum)
#if [ "$DATASET" = "AD" ]; then
#    BAM_PATTERN="${TIMEGROUP}_T2"
#else
#    BAM_PATTERN="${TIMEGROUP}"
#fi

#BAM_FILES=$(ls "${BAM_DIR}"/*${BAM_PATTERN}*.sorted.bam 2>/dev/null || true)
#if [ -n "$BAM_FILES" ]; then
#    N_BAMS=$(echo "$BAM_FILES" | wc -l)
#    echo "Adding BAM coverage filter ($N_BAMS replicates)..."
#    for bam in $BAM_FILES; do
#        samtools idxstats "$bam" | awk '$3 > 0 {print $1}'
#    done | sort -u >> "$EXPRESSED_CONTIGS"
#    sort -u "$EXPRESSED_CONTIGS" -o "$EXPRESSED_CONTIGS"
#    echo "Total expressed contigs after BAM filter: $(wc -l < $EXPRESSED_CONTIGS)"
#else
#    echo "No BAM files found for pattern *${BAM_PATTERN}* — using SnapT contigs only"
#fi

# ── STEP 2: Filter CDS to expressed contigs + cap at 150 nt ──────────────────
CDS_BED="${INTARNA_DIR}/cds_regions.bed"
CDS_FILTERED="${INTARNA_DIR}/cds_filtered.bed"
MRNA_FA="${INTARNA_DIR}/mrna_targets.fa"

if [ ! -f "$CDS_BED" ]; then
    echo "Extracting CDS from prodigal.gff..."
    awk '$3=="CDS"' "${ANNOTATION_DIR}/prodigal.gff" \
        | awk 'BEGIN{OFS="\t"}{print $1,$4-1,$5,$9,0,$7}' \
        > "$CDS_BED"
    echo "Total CDS: $(wc -l < $CDS_BED)"
fi

grep -Ff "$EXPRESSED_CONTIGS" "$CDS_BED" \
    | awk 'BEGIN{OFS="\t"}{
        if ($6=="+") print $1,$2,($2+150<$3?$2+150:$3),$4,$5,$6
        else         print $1,($3-150>$2?$3-150:$2),$3,$4,$5,$6
    }' > "$CDS_FILTERED"
echo "Filtered CDS: $(wc -l < $CDS_FILTERED)  (from $(wc -l < $CDS_BED) total)"

bedtools getfasta \
    -fi "$ASSEMBLY" \
    -bed "$CDS_FILTERED" \
    -s -name \
    > "$MRNA_FA"
echo "mRNA targets: $(grep -c '^>' $MRNA_FA)"

# ── STEP 3: Run IntaRNA ───────────────────────────────────────────────────────
INTERACTIONS="${INTARNA_DIR}/interactions.tsv"

if [ -f "$INTERACTIONS" ]; then
    echo "IntaRNA output exists — skipping (delete $INTERACTIONS to rerun)"
else
    N_QUERY=$(grep -c '^>' $SRNA_FA)
    N_TARGET=$(grep -c '^>' $MRNA_FA)
    echo "Running IntaRNA: ${N_QUERY} queries x ${N_TARGET} targets..."
    IntaRNA \
	    -q "$SRNA_FA" \
	    -t "$MRNA_FA" \
	    --outMode C \
	    --outCsvCols "id1,id2,E,ED1,ED2,Pu1,Pu2,hybridDP,start1,end1,start2,end2,P_E" \
	    --out "$INTERACTIONS" \
	    --outSep ";" \
	    --noSeed \
	    --outMaxE 0 \
	    --outDeltaE 100 \
	    -n 100 \
	    --threads 48
    echo "Raw interactions: $(( $(wc -l < $INTERACTIONS) - 1 )) rows"
fi

# ── STEP 4: Filter by p-value and summarize ───────────────────────────────────
SUMMARY="${INTARNA_DIR}/intarna_summary.tsv"
SIGNIFICANT="${INTARNA_DIR}/interactions_significant.tsv"

echo "$INTERACTIONS" > /tmp/intarna_in.txt
echo "$SUMMARY"      > /tmp/intarna_out.txt
echo "$SIGNIFICANT"  > /tmp/intarna_sig.txt
echo "$PVAL_THRESH"  > /tmp/intarna_pval.txt

python3 - <<'PYEOF'
import pandas as pd
import numpy as np

interactions = open("/tmp/intarna_in.txt").read().strip()
summary_path = open("/tmp/intarna_out.txt").read().strip()
sig_path     = open("/tmp/intarna_sig.txt").read().strip()
pval_thresh  = float(open("/tmp/intarna_pval.txt").read().strip())

df = pd.read_csv(interactions, sep=";")
df.columns = [c.strip() for c in df.columns]
df = df.rename(columns={"id1": "srna_id", "id2": "mrna_id",
                         "E": "energy", "P_E": "pvalue"})
df["energy"] = pd.to_numeric(df["energy"], errors="coerce")
df["pvalue"] = pd.to_numeric(df["pvalue"], errors="coerce")
df = df.dropna(subset=["energy"])

print(f"Total interactions: {len(df)}")

if "pvalue" in df.columns and df["pvalue"].notna().any():
    print(f"P_E range: {df['pvalue'].min():.4f} to {df['pvalue'].max():.4f}")
    df_sig = df[df["pvalue"] < pval_thresh].copy()
    print(f"Significant (p<{pval_thresh}): {len(df_sig)} interactions "
          f"across {df_sig['srna_id'].nunique()} sRNAs")
    if len(df_sig) == 0:
        print("WARNING: No significant hits — falling back to top 100 by energy")
        df_sig = df.sort_values("energy").groupby("srna_id").head(100)
else:
    print("WARNING: P_E column missing — falling back to top 100 by energy")
    df_sig = df.sort_values("energy").groupby("srna_id").head(100)
df_sig.sort_values(["srna_id", "energy"]).to_csv(sig_path, sep="\t", index=False)

summary = df_sig.groupby("srna_id").agg(
    n_targets      = ("mrna_id", "nunique"),
    best_energy    = ("energy",  "min"),
    mean_energy    = ("energy",  "mean"),
    std_energy     = ("energy",  "std"),
    n_interactions = ("mrna_id", "count"),
    **( {"best_pvalue": ("pvalue", "min"),
         "mean_pvalue": ("pvalue", "mean")}
        if "pvalue" in df_sig.columns else {} )
).reset_index()

summary.to_csv(summary_path, sep="\t", index=False)
print(f"\nSummary: {len(summary)} sRNAs → {summary_path}")
print(f"Sig hits → {sig_path}")
PYEOF

echo "===== Done: $DATASET $TIMEGROUP ====="
echo "  Summary:   $SUMMARY"
echo "  Sig hits:  $SIGNIFICANT"
