import pandas as pd
import sys

readcounts = open(sys.argv[1]+'_viral_read_counts.tsv').readlines()
allhits = pd.read_csv(sys.argv[1]+'_all_hits.csv')
original_cols = list(allhits)

readcol = []
for cont in allhits['Contig']:
    line = list(filter(lambda x: x.startswith('NODE_'+str(cont)+'_'), readcounts))[0]
    readcol.append(int(line.split('\t')[2]))
allhits['Reads'] = readcol
allhits['Avg Read Depth'] = allhits['Reads'] / allhits['Query Length']

neworder = original_cols[:-1] + ['Reads', 'Avg Read Depth'] + \
           [original_cols[-1]]
allhits = allhits[neworder]
allhits.to_csv(sys.argv[1]+'_all_hits.csv', index=False, header=list(allhits))
