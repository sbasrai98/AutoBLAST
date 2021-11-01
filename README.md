# viral-RNAseq-pipeline
A pipeline to automate viral load characterization from RNAseq data (rnaseq_pipeline.sh). It is 
configured to run on Compute Canada's HPC clusters.

# Workflow

 * Filter out the host genome (BBDuk)
 * Assemble contigs from clean reads (SPAdes)
 * BLAST contigs
 * Extract viral sequences, generate summary (digestpd.py)
 * Download BLAST hits, align contigs, build scaffold/consensus sequences, generate coverage diagrams (tie.py, contigmap.py)
 * Map reads to viral contigs, get read counts (BBMap, SAMtools)
	  - Post-assembly mapped read counts can be useful as a rough proxy for relative abundance, but this is not always guaranteed to be accurate!
 * Predict protein domains (HMMER)

