#!/bin/bash
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 02:00:00
#SBATCH --array=22,28

bash $PY/rnaseq_pipeline.sh "$SLURM_ARRAY_TASK_ID"_R1.fastq.gz "$SLURM_ARRAY_TASK_ID"_R2.fastq.gz GCF_000003745.3_12X_genomic.fna $SLURM_ARRAY_TASK_ID
