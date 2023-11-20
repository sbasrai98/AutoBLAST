# %%
import argparse
import pandas as pd
from Bio import SeqIO

def extract_contigs(parsed_hits: str, contigs: str) -> list:
    """Extract all contigs of interest from an assembly.

    parsed_hits: path to tab-delimited BLAST hits file
    contigs: path to contig assembly FASTA file

    Returns a list of SeqRecord objects."""

    assembly = list(SeqIO.parse(contigs, 'fasta'))
    hits = pd.read_table(parsed_hits)
    extracted = [x for x in assembly if x.id.split('_')[1] in list(map(str, hits['Contig']))]
    return extracted

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("parsed_hits", help="tab-delimited BLAST hits file")
    parser.add_argument("contigs", help="contig assembly FASTA file")
    parser.add_argument("output", help="name of output file to write")
    args = parser.parse_args()

    extracted_contigs = extract_contigs(args.parsed_hits, args.contigs)
    SeqIO.write(extracted_contigs, args.output, 'fasta')
