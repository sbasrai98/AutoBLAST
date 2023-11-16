#!/bin/bash -l
#$ -cwd
#$ -S /bin/bash

# $1 is SRA ID

date +%r
prefetch $1

date +%r
fasterq-dump $1

date +%r

# while read line; do qsub -P abelsonlab -l h_vmem=16G,h_rt=1:0:0 -q u20build getSRA.sh $line; done < sra_samples.txt