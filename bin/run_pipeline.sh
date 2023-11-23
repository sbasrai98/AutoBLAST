#!/bin/bash -l
#$ -cwd
#$ -S /bin/bash

date +%r

source /u/sbasrai/mambaforge/etc/profile.d/conda.sh
conda activate ./.venv/
module load java/11

snakemake -p --keep-incomplete --latency-wait 60 -j 100 \
--cluster "qsub -P abelsonlab -l h_vmem={resources.memory},h_rt={resources.time} -pe smp {threads} -V \
-e /.mounts/labs/abelsonlab/private/salman/viral-RNAseq-pipeline/pipetest/job_err/{wildcards.sample}_{rule}_err.txt \
-o /.mounts/labs/abelsonlab/private/salman/viral-RNAseq-pipeline/pipetest/job_out/{wildcards.sample}_{rule}_out.txt"

date +%r
