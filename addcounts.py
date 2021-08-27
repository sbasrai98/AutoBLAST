import pandas as pd
import sys

readcounts = open(sys.argv[1]+'_viral_read_counts.tsv').readlines()
allhits = pd.read_csv(sys.argv[1]+'_all_hits.csv')

readcol = []
for cont in allhits['Contig']:
    line = list(filter(lambda x: x.startswith('NODE_'+str(cont)+'_'), readcounts))[0]
    readcol.append(int(line.split('\t')[2]))
rfcol = []
for r in readcol:
    rfcol.append(r/sum(readcol)*100)
rdcol = []
for i in range(len(allhits)):
    rdcol.append(readcol[i]/allhits['Query Length'][i])
allhits['Reads'] = readcol
allhits['Relative Fraction'] = rfcol
allhits['Avg Read Depth'] = rdcol

allhits.to_csv(sys.argv[1]+'_all_hits.csv', index=False, header=list(allhits))


allhits = pd.read_csv(sys.argv[1]+'_all_hits.csv')

sumhits = pd.read_csv(sys.argv[1]+'_summary.csv')

sumreads = []
sumrdcol = []
for virus in sumhits['Virus']:
    selection = allhits[allhits['Virus']==virus]
    sr = sum(selection['Reads'])
    sumreads.append(sr)
    totlen = sum(selection['Query Length'])
    sumrdcol.append(sr/totlen)
sumrfcol = []
for sr in sumreads:
    sumrfcol.append(sr/sum(sumreads)*100)

sumhits['Total Reads'] = sumreads
sumhits['Relative Fraction'] = sumrfcol
sumhits['Avg Read Depth'] = sumrdcol

sumhits.to_csv(sys.argv[1]+'_summary.csv', index=False, header=['']+list(sumhits)[1:])
