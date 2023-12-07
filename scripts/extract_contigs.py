# %%
import os
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
    parser.add_argument("output", help="name of output folder to write")
    args = parser.parse_args()

    os.makedirs(args.output)
    extracted_contigs = extract_contigs(args.parsed_hits, args.contigs)
    for c in extracted_contigs:
        folder = args.output.rstrip('/').split('/')[-1] 
        c.id = folder + c.id[4:]
        c.description = ''
        filename = '_'.join(c.id.split('_')[:2]) + '.fa'
        SeqIO.write(c, f'{args.output}/{filename}', 'fasta')
