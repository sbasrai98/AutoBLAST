# $1: contigs in FASTA format                                           
# $2: only 'spade' accepted                                             
# $3: sample name (used in output)                                      
# $4: parsed blast results file

python3 scripts/download_hits_align_contigs.py $1 $2 $3 $4

tar -czf "$3".tar.gz $3

