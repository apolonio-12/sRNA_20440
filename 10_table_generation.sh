#!/bin/bash
#SBATCH -t 6:00:00
#SBATCH -n 36
#SBATCH -J download_reads
#SBATCH -p mit_normal
#SBATCH --mem=200GB
#SBATCH --output=/home/apolonio/orcd/scratch/20_440/project_logs/table_generation.log
source /home/apolonio/miniconda3/etc/profile.d/conda.sh
conda activate snapt_env

python table_generation.py
