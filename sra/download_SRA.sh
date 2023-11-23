#!/bin/bash -l
#$ -cwd
#$ -S /bin/bash

# $1 is SRA ID

date +%r
if [ ! -d sra_dump/"$1" ]; then
	prefetch $1 -O sra_dump/

	date +%r
	fasterq-dump sra_dump/$1 -O sra_dump/
fi

# gzip the fastqs
if [ -f sra_dump/"$1"_1.fastq ]; then
	gzip sra_dump/"$1"_1.fastq
fi
if [ -f sra_dump/"$1"_2.fastq ]; then
	gzip sra_dump/"$1"_2.fastq
fi

date +%r

# while read line; do qsub -P abelsonlab -l h_vmem=2G,h_rt=0:70:0 -q u20build \
-e /.mounts/labs/abelsonlab/private/salman/viral-RNAseq-pipeline/sra/job_oe/"$line"_err.txt \
-o /.mounts/labs/abelsonlab/private/salman/viral-RNAseq-pipeline/sra/job_oe/"$line"_out.txt \
download_SRA.sh $line; done < sra_samples.txt