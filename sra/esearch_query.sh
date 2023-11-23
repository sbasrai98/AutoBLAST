#!/bin/bash -l
#$ -cwd
#$ -S /bin/bash

time esearch -db sra -query "vitis[Organism] AND \"rna seq\"[Strategy] AND illumina[Platform]" | efetch -format native -mode xml > all_vitis_RNAseq.xml
