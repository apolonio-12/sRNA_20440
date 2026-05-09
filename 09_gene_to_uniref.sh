#!/bin/bash
#SBATCH -t 6:00:00
#SBATCH -n 48
#SBATCH -J download_reads
#SBATCH -p mit_normal
#SBATCH --mem=200GB
#SBATCH --output=/home/apolonio/orcd/scratch/20_440/project_logs/diamond.log

source /home/apolonio/miniconda3/etc/profile.d/conda.sh
conda activate /home/apolonio/.conda/envs/humann_env

diamond blastp \
  -q AD/annotation/prodigal.faa \
  -d /home/apolonio/humann_db/uniref/humann4_protein_database_filtered_v2019_06.dmnd \
  -o gene_to_uniref90.tsv \
  --outfmt 6 qseqid sseqid pident length evalue bitscore \
  --max-target-seqs 1 \
  --evalue 1e-5
