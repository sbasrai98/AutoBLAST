#configfile: "config.yaml"

with open('pipetest/sra_samples1.txt') as input_file:
	inputs = [x.rstrip() for x in input_file.readlines()]

rule all:
	input:
		expand("pipetest/parsed/{sample}_parsed_hits.txt", sample=inputs)

rule assemble_reads:
	input:
		r1="pipetest/fastq/{sample}_1.fastq",
		r2="pipetest/fastq/{sample}_2.fastq"
	output:
		folder=directory("pipetest/assembly/{sample}_assembly"),
		contigs="pipetest/assembly/{sample}_assembly/contigs.fasta"
	threads: 16
	resources:
		memory="8G",
		time="4:0:0"
	shell:
		"time spades.py -t {threads} -1 {input.r1} -2 {input.r2} --rnaviral -o {output.folder}"

rule BLAST_contigs:
	input:
		"pipetest/assembly/{sample}_assembly/contigs.fasta"
	output:
		"pipetest/blasted/{sample}_contigs_blasted.txt"
	threads: 24
	resources:
		memory="16G",
		time="12:0:0"
	shell:
		'time blastn -db nt -query {input} -max_target_seqs 1 -out {output} -evalue 1e-60 -num_threads {threads} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue qlen slen staxids sscinames scomnames sblastnames sskingdoms stitle salltitles"; '
		"sed -i '1s/^/qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tqlen\tslen\tstaxids\tsscinames\tscomnames\tsblastnames\tsskingdoms\tstitle\tsalltitles\
		/' {output}"

rule parse_BLAST_results:
	input:
		"pipetest/blasted/{sample}_contigs_blasted.txt"
	output:
		"pipetest/parsed/{sample}_parsed_hits.txt"
	threads: 1
	resources:
		memory="8G",
		time="0:30:0"
	shell:
		"time python3 scripts/parse_BLAST_results.py {input} {output}"

rule download_hits_align_contigs:
	input:
		contigs="pipetest/inputs/{sample}_contigs.fa",
		parsed_hits="pipetest/parsed/{sample}_parsed_hits.txt"
	output:
		directory("pipetest/results/{sample}")
	shell:
#		"bash scripts/dl_wrapper.sh {input.contigs} spade pipetest/zipped_results/{wildcards.sample}/ {input.parsed_hits}"
		"python3 scripts/download_hits_align_contigs.py {input.contigs} spade {output} {input.parsed_hits}"


# sample_22

#bash scripts/dl_wrapper.sh pipetest/inputs/sample_22_contigs.fa spade pipetest/zipped_results/sample_22 pipetest/parsed/sample_22_parsed_hits.txt