echo ------------------------
echo BLAST input: $1
echo BLAST output: $2
echo Sample ID: $3
echo ------------------------
echo Preparing input file...
time python3 qton.py $1
echo " "
echo Running BLAST...
time blastn -db nt -query $1 -max_target_seqs 1 -out $2 -evalue 1e-60 -num_threads 2 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue qlen slen staxids sscinames scomnames sblastnames sskingdoms stitle salltitles" 
#sed -i '1s/^/qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tqlen\tslen\tstaxids\tsscinames\tscomnames\tsblastnames\tsskingdoms\tstitle\tsalltitles\n/' $2
echo " "
echo Extracting viral sequences, summarizing results...
time python3 digest.py $2 $3
echo " "
echo Summary written to: "$3"_blast_summary.csv
