import pandas as pd
import os, re, sys
from Bio import SeqIO, Align, Entrez
Entrez.email = 'sbasrai@uwaterloo.ca'
from Bio.Seq import Seq
from PIL import Image, ImageDraw
from contigselection import write, nuv, prime
from contigmap import make_aligner, map_seqs, build_msa, consensus, contig_diagram

if len(sys.argv) != 5:
    print('Enter: python3 '+sys.argv[0]+' <contigs.fa> <contig_type> <sample_name> <sample_blast_hits.csv>') 
    quit()

contig_file = sys.argv[1] #'sample16contigs.fasta'
cont_file = open(contig_file)
contig_data = cont_file.read()
cont_file.close()

if sys.argv[2] == 'nuv': pull_contigs = nuv
if sys.argv[2] == 'prime': pull_contigs = prime

cfa = ' contigs.fa'
smp = sys.argv[3]  #'sample16'
os.makedirs(smp) # will use sample name/id

hits = pd.read_csv(sys.argv[4]) #'16blasthits.csv')
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

        # contig map  ...could clean up...
        contigs = list( SeqIO.parse(smp+'/'+v+'/'+a+'/'+a+cfa, 'fasta') )
        ref = SeqIO.read(smp+'/'+v+'/'+a+'/'+a+'.fa', 'fasta')
        name = a+' contigs to '+str(ref.id)
        posns = map_seqs(ref, contigs)
        if len(posns) > 0:
            msa = build_msa(posns, contigs)
            SeqIO.write(msa, smp+'/'+v+'/'+a+'/'+name+' msa.fa', 'fasta')
            msacollapse = consensus(msa)
            SeqIO.write(msacollapse, smp+'/'+v+'/'+a+'/'+name+' consensus.fa', 'fasta')
            im = contig_diagram(posns, ref)
            im.save(smp+'/'+v+'/'+a+'/'+name+' diagram.png', quality=100)
