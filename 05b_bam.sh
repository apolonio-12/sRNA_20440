#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH -n 30
#SBATCH -J bam_files
#SBATCH -p mit_normal
#SBATCH --mem=200GB
#SBATCH --output=/home/apolonio/orcd/scratch/20_440/project_logs/bam.log

set -euo pipefail

BASE=~/orcd/scratch/20_440
IN_DIR=$BASE/AD
OUT_DIR=$BASE/AD/bam_files
INDEX=$BASE/AD/coassembly/bowtie2_index

mkdir -p $OUT_DIR

for d in $IN_DIR/AD_S1_*_RNA_knead_out; do
    sample=$(basename $d | sed 's/_knead_out//')

    R1=("$d"/*_kneaddata_paired_1.fastq)
    R2=("$d"/*_kneaddata_paired_2.fastq)

    echo "Processing $sample"

    bowtie2 -x "$INDEX" \
        -1 "${R1[0]}" \
        -2 "${R2[0]}" \
        -p 24 2> "$OUT_DIR/${sample}_bowtie2.log" | \
        samtools view -bS | \
        samtools sort -@ 8 -o "$OUT_DIR/${sample}.sorted.bam"

    samtools index "$OUT_DIR/${sample}.sorted.bam"

    echo "$sample complete."
done
