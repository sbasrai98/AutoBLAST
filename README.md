# viral-RNAseq-pipeline
A pipeline to automate viral load characterization from RNAseq data.

# Dependencies
 * Python
	* pandas
	* biopython
	* pillow
 * External
	* BBmap
	* SPAdes
	* BLAST+

# Workflow
 * Filter out the host genome (BBDuk)
 * Assemble contigs from clean reads (SPAdes)
 * BLAST contigs
 * Extract viral sequences, generate summary
 * Download BLAST hits, align contigs, build scaffold/consensus sequences, generate coverage diagrams
