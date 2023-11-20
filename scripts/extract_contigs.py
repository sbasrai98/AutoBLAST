import re
import argparse
from Bio import SeqIO

def extract_contigs(parsed_hits: str, contigs: str) ->


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Enter: python3 "+sys.argv[0]+" <contigs.fa> <contig_type> <output.fa> <list_of_contigs>")
        quit()
    
    fin = open(sys.argv[1])
    contig_data = fin.read()
    fin.close()
    contig_list = sys.argv[4].split(',')

    if sys.argv[2] == 'nuv': pull_contigs = nuv
    if sys.argv[2] == 'prime': pull_contigs = prime
    if sys.argv[2] == 'spade': pull_contigs = spade
    pull_contigs(sys.argv[3], contig_data, contig_list)
