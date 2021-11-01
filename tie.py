import pandas as pd
import os, re, sys
from Bio import SeqIO, Align, Entrez
Entrez.email = 'sbasrai@uwaterloo.ca'
from Bio.Seq import Seq
from PIL import Image, ImageDraw
from contigselection import write, nuv, prime, spade, nuvselect, primeselect, spadeselect
from contigmap import make_aligner, map_seqs, build_msa, consensus, contig_diagram

if len(sys.argv) != 5:
    print('Enter: python3 '+sys.argv[0]+' <contigs.fa> <contig_type> <sample_name> <sample_blast_hits.csv>') 
    quit()

contig_file = sys.argv[1]
cont_file = open(contig_file)
contig_data = cont_file.read() # reassigned this variable below..
cont_file.close()

if sys.argv[2] == 'nuv': pull_contigs = nuv
if sys.argv[2] == 'prime': pull_contigs = prime
if sys.argv[2] == 'spade': pull_contigs = spade

cfa = '_contigs.fa'
smp = sys.argv[3]
os.makedirs(smp)

hits = pd.read_csv(sys.argv[4])
all_conts = list(hits['Contig'])
if sys.argv[2] == 'nuv': contig_data = nuvselect(contig_data, all_conts)
if sys.argv[2] == 'prime': contig_data =  primeselect(contig_data, all_conts)
if sys.argv[2] == 'spade': contig_data =  spadeselect(contig_data, all_conts)
write(smp+'/'+smp+cfa, contig_data)
viruses = list(dict.fromkeys(hits['Virus']))

for v in viruses: 
    os.makedirs(smp+'/'+v)
    cur_virus = hits[hits['Virus'] == v]
    vir_conts = list(cur_virus['Contig'])
    pull_contigs(smp+'/'+v+'/'+v+cfa, contig_data, vir_conts)
    access = list(dict.fromkeys(cur_virus['Accession']))

    for a in access:
        cur_acc = cur_virus[cur_virus['Accession'] == a]
        acc_conts = list(cur_acc['Contig'])
        os.makedirs(smp+'/'+v+'/'+a)
        pull_contigs(smp+'/'+v+'/'+a+'/'+a+cfa, contig_data, acc_conts)
        handle = Entrez.efetch(db="nucleotide", id=a, rettype="fasta", retmode="text")
        write(smp+'/'+v+'/'+a+'/'+a+'.fa', handle.read())

        # contig map
        contigs = list( SeqIO.parse(smp+'/'+v+'/'+a+'/'+a+cfa, 'fasta') )
        ref = SeqIO.read(smp+'/'+v+'/'+a+'/'+a+'.fa', 'fasta')
        name = a+'_contigs_to_'+str(ref.id)
        posns = map_seqs(ref, contigs)
        if len(posns) > 0:
            msa = build_msa(posns, contigs)
            SeqIO.write(msa, smp+'/'+v+'/'+a+'/'+name+' msa.fa', 'fasta')
            if len(posns) > 1:
                msacollapse = consensus(msa)
                SeqIO.write(msacollapse, smp+'/'+v+'/'+a+'/'+name+'_consensus.fa', 'fasta')
            im = contig_diagram(posns, ref)
            im.save(smp+'/'+v+'/'+a+'/'+name+'_diagram.png', quality=100)
