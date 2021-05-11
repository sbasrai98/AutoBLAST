import pandas as pd 
import re

tbl = pd.read_table('22test.txt')

keepcol = tbl[ ['qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen', 'sscinames', 'scomnames', 'sblastnames', 'sskingdoms', 'stitle'] ]
keepcol = keepcol.drop_duplicates(subset=['qseqid'])
filter = 'Viruses'
keepcol = keepcol[(keepcol['sscinames'] == filter) |
               (keepcol['scomnames'] == filter) |
               (keepcol['sblastnames'] == filter) |
               (keepcol['sskingdoms'] == filter)]

keepcol.sort_values(['sscinames', 'pident', 'length'], ascending= [True, False, False], inplace=True)
keepcol = keepcol[ ['qseqid', 'sscinames', 'pident', 'length', 'qlen', 'slen', 'sseqid', 'stitle'] ]

contig_num = re.compile('_[0-9]+')
keepcol['qseqid'] = keepcol['qseqid'].map(lambda x: contig_num.findall(x)[0][1:])
access = re.compile('\|[^a-z]+\.[0-9]+')
keepcol['sseqid'] = keepcol['sseqid'].map(lambda x: access.findall(x)[0][1:])

summary = {'': [], 'Virus': [], 'Contigs': [], 'Length': [], 'Identity': []}

for n in list(dict.fromkeys(keepcol['sscinames'])):
    onevirus = keepcol[keepcol['sscinames'] == n]
    summary['Virus'].append(n)
    summary['Contigs'].append(onevirus.shape[0])
    summary['Length'].append(str(onevirus['qlen'].min())+' - '+str(onevirus['qlen'].max()))
    summary['Identity'].append(str(onevirus['pident'].min())+' - '+str(onevirus['pident'].max()))
    summary[''].append('')
summary = pd.DataFrame(summary, columns= list(summary))
summary.sort_values(['Contigs', 'Length', 'Identity'], ascending= [False, True, False], inplace=True)

fout = open('sample 22 virus summary.csv', 'w')
fout.write('BLAST Summary:\n')
summary.to_csv(fout, index=False)
fout.write('\nAll hits:\n')
header = ['Contig','Virus','Identity','Hit Length','Query Length','Subject Length','Accession','Title']
fout.write(','.join(header)+'\n')
keepcol.to_csv(fout, index = False, header=False) 
fout.close()