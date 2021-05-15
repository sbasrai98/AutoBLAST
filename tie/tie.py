###########################
from Bio import Entrez
Entrez.email = 'sbasrai@uwaterloo.ca'
###########################
import pandas as pd
import os
import re
from nuvselfunc import nuvselect

os.chdir('tie')

contig_file = 'NuVs-22-Contigs.fa'
cont_file = open(contig_file)
contig_data = cont_file.read()
cont_file.close()

cfa = ' contigs.fa'
smp = 'hitsort'
os.mkdir(smp) # will use sample name/id
#os.chdir('hitsort')

# contig file for all sample hits
fout = open(smp+'/'+smp+cfa, 'w')
fout.write(contig_data)
fout.close()

hits = pd.read_csv('22blasthits.csv')

viruses = list(dict.fromkeys(hits['Virus']))
for v in viruses: 
    os.makedirs(smp+'/'+v)
    cur_virus = hits[hits['Virus'] == v]
    vir_conts = list(cur_virus['Contig'])
    fout = open(smp+'/'+v+'/'+v+cfa, 'w')
    fout.write(nuvselect(contig_data, vir_conts))
    fout.close() 
    access = list(dict.fromkeys(cur_virus['Accession']))
    for a in access:
        cur_acc = cur_virus[cur_virus['Accession'] == a]
        acc_conts = list(cur_acc['Contig'])
        os.makedirs(smp+'/'+v+'/'+a)
        fout = open(smp+'/'+v+'/'+a+'/'+a+cfa, 'w')
        fout.write(nuvselect(contig_data, acc_conts))
        fout.close() 
        handle = Entrez.efetch(db="nucleotide", id=a, rettype="fasta", retmode="text")
        fout = open(smp+'/'+v+'/'+a+'/'+a+'.fa', 'w')
        fout.write(handle.read())
        fout.close()
        



        
