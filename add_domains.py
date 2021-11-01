import pandas as pd
import sys
import re

hmmer_hits = pd.read_table(sys.argv[1])
significant = hmmer_hits[hmmer_hits['full_E-value'] < float(1e-09)]
contigs = list(map(lambda x: x[:-7], significant['qname']))
contigs = list(dict.fromkeys(contigs))

annotations = {}
for id in contigs:
    selection = significant[significant['qname'].str.startswith(id)]
    selection = selection.sort_values(by='full_E-value')
    domains = list(dict.fromkeys(selection['tname']))
    contig_num = int(re.search('NODE_[0-9]+_', id).group(0)[5:-1])
    annotations[contig_num] = '|'.join(domains)

# add annotation column to all_hits.csv
all_hits = pd.read_csv(sys.argv[2])
cols = list(all_hits)

predicted_domains = []
for num in all_hits['Contig']:
    if num in annotations:
        predicted_domains.append(annotations[num])
    else:
        predicted_domains.append(' ')
all_hits['Predicted Domains'] = predicted_domains


neworder = cols[:-1] + ['Predicted Domains', cols[-1]]
all_hits = all_hits[neworder]
all_hits.to_csv(sys.argv[2], index=False, header=list(all_hits))
