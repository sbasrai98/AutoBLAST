import sys

data = open(sys.argv[1]).readlines()
data = list(filter(lambda x: not(x.startswith('#')), data))

newtsv = []
for line in data:
    values = []
    line = line.split(' ')
    line = list(filter(lambda x: x != '', line))
    line = line[:22] + [' '.join(line[22:])]
    newtsv.append('\t'.join(line))

outname = sys.argv[1][::-1]
ext = outname.find('.')+1
outname = outname[ext:][::-1]+'_parsed.tsv'

fout = open(outname, 'w')
headings = ['tname', 'taccess', 'tlen', 'qname', 'qaccess', 'qlen', 'full_E-value', 'full_score', 'full_bias', '#', 'of', 'c-Evalue', 
            'i-Evalue', 'dom_score', 'dom_bias', 'hmmfrom', 'hmmto', 'alifrom', 'alito', 'envfrom', 'envto', 'acc', 'description']
fout.write('\t'.join(headings)+'\n')
for line in newtsv: fout.write(line)
fout.close()
