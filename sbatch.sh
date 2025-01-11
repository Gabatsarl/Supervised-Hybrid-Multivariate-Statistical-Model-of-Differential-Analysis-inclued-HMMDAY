#!/bin/bash
#SBATCH --account=def-amadou_cpu
#SBATCH --time=00-20:00
##SBATCH --gpus-per-node=2
#SBATCH --mem=20G
#SBATCH --job-name="name-of-job"
#SBATCH --mail-user=........@gmail.com
#SBATCH --mail-type=ALL

module load r/4.2.2
Rscript DESCRIPTIVE.R
