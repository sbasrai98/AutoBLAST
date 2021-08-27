echo ---------------------------
echo Left reads: $1
echo Right reads: $2
echo Host genome: $3
echo Sample ID: $4
echo ---------------------------

date +%r
echo Filtering host genome from reads...
module load bbmap
bbduk.sh in1=$1 in2=$2 out1="$4"_R1_clean.fq out2="$4"_R2_clean.fq ref=$3 -Xmx50g

date +%r
echo Assembling clean reads...
module load spades
spades.py -1 "$4"_R1_clean.fq -2 "$4"_R2_clean.fq --rnaviral -o $4

date +%r
module load StdEnv/2020
module load gcc/9.3.0
module load blast+/2.11.0

echo Running BLAST...
time blastn -db nt -query $4/contigs.fasta -max_target_seqs 1 -out "$4"_blasted.txt -evalue 1e-60 -num_threads 16 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue qlen slen staxids sscinames scomnames sblastnames sskingdoms stitle salltitles" 
sed -i '1s/^/qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tqlen\tslen\tstaxids\tsscinames\tscomnames\tsblastnames\tsskingdoms\tstitle\tsalltitles\n/' "$4"_blasted.txt

date +%r
module load scipy-stack/2020b
echo Extracting viral sequences, summarizing results...
python3 $PY/digestpd.py "$4"_blasted.txt $4
echo Summary written to: "$4"_summary.csv
echo All hits written to: "$4"_all_hits.csv

date +%r
# LOADING ENV WILL NEED TO BE PUT IN WRAPPER SCRIPT.
echo Loading Python virtual environment...
module load python/3.8.10
virtualenv --no-download tempenv
source tempenv/bin/activate
pip install --no-index --upgrade pip
pip install --no-index pillow
pip install --no-index pandas
pip install --no-index biopython

date +%r
echo Aligning contigs to NCBI sequences, generating scaffolds and diagrams...
python3 $PY/tie.py $4/contigs.fasta spade "$4"_virus_alignments "$4"_all_hits.csv > "$4"_tie_report.txt
deactivate
rm -rf tempenv

date +%r
echo Mapping reads to viral contigs...
module load bbmap
bbmap.sh in1="$4"_R1_clean.fq in2="$4"_R2_clean.fq ref=$4/contigs.fasta out="$4"_mapped_reads.sam

date +%r
echo Converting SAM file, counting mapped/unmapped reads...
module load samtools
# remove all unmapped reads:
sed -i '/\*[\t]0[\t]0[\t]\*[\t]\*[\t]0[\t]0/d' "$4"_mapped_reads.sam
# convert SAM to BAM. -S: input is a SAM file, -b: produce a BAM file as output
samtools view -S -b "$4"_mapped_reads.sam > "$4"_mapped_reads.bam
# sort the file (order positionally based upon their alignment coordinates on each chromosome)
samtools sort "$4"_mapped_reads.bam -o "$4"_mapped_reads.sort.bam
# index bam and then get stats; reads mapped for each contig
samtools index "$4"_mapped_reads.sort.bam
samtools idxstats "$4"_mapped_reads.sort.bam > "$4"_viral_read_counts.tsv
sed -i '1s/^/contig\tlength\tmapped\tunmapped\n/' "$4"_viral_read_counts.tsv
python3 $PY/addcounts.py "$4"
mkdir "$4"_read_mapping
mv "$4"_mapped_reads.* "$4"_read_mapping
mv "$4"_viral_read_counts.tsv "$4"_read_mapping

date +%r
echo Creating output directory
mkdir "$4"_processed
# mv $1 "$4"_processed
# mv $2 "$4"_processed
mv "$4"_R1_clean.fq "$4"_processed
mv "$4"_R2_clean.fq "$4"_processed
mv $4 "$4"_processed/"$4"_assembly
mv "$4"_blasted.txt "$4"_processed
mv "$4"_summary.csv "$4"_processed
mv "$4"_all_hits.csv "$4"_processed
mv "$4"_read_mapping "$4"_processed
date +%r