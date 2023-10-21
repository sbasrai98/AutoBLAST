import pandas as pd 
import re
import sys

# blast results file + filter(string) -> DataFrame trimmed/formatted 
def filter_hits(blast_results_file, filter):
    tbl = pd.read_table(blast_results_file)
    keepcol = tbl[ ['qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen', 'sscinames', 'scomnames', 'sblastnames', 'sskingdoms', 'stitle'] ]
    keepcol = keepcol.drop_duplicates(subset=['qseqid'])
    keepcol = keepcol[(keepcol['sscinames'] == filter) |
                (keepcol['scomnames'] == filter) |
                (keepcol['sblastnames'] == filter) |
                (keepcol['sskingdoms'] == filter)]
    keepcol.sort_values(['sscinames', 'pident', 'length'], ascending= [True, False, False], inplace=True)
    keepcol = keepcol[ ['qseqid', 'sscinames', 'pident', 'length', 'qlen', 'slen', 'sseqid', 'stitle'] ]
    contig_num = re.compile('_[0-9]+') # get contig numbers
    keepcol['qseqid'] = keepcol['qseqid'].map(lambda x: contig_num.findall(x)[0][1:])
    access = re.compile('\|[^a-z]+\.[0-9]+') # get accession IDs
    keepcol['sseqid'] = keepcol['sseqid'].map(lambda x: access.findall(x)[0][1:])
    return keepcol

# DataFrame(all hits) -> DataFrame(summary)
def summarize_hits(hits):
    summary = {'': [], 'Virus': [], 'Contigs': [], 'Length': [], 'Identity': [], 'mxlen': [], 'mxid': []}
    for n in list(dict.fromkeys(hits['sscinames'])):
        onevirus = hits[hits['sscinames'] == n] # DataFrame for each virus
        summary['Virus'].append(n)
        summary['Contigs'].append(onevirus.shape[0])
        summary['Length'].append(str(onevirus['qlen'].min())+' - '+str(onevirus['qlen'].max()))
        summary['Identity'].append(str(onevirus['pident'].min())+' - '+str(onevirus['pident'].max()))
        summary[''].append('')
        summary['mxlen'].append(onevirus['qlen'].max())
        summary['mxid'].append(onevirus['pident'].max())
    summary = pd.DataFrame(summary, columns=list(summary))
    summary.sort_values(['Contigs', 'mxlen', 'mxid'], ascending= [False, False, False], inplace=True)
    summary.drop(labels=['mxlen', 'mxid'], axis=1, inplace=True)
    return summary

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Enter: python3 "+sys.argv[0]+" <blastresults.txt> <sample_name>")
        quit()

    all_hits = filter_hits(sys.argv[1], 'Viruses') #'test files/22test.txt'
    header = ['Contig','Virus','Identity','Hit Length','Query Length','Subject Length','Accession','Title']
    all_hits.to_csv(sys.argv[2]+'_all_hits.csv', index=False, header=header)

    sum_hits = summarize_hits(all_hits)
    sum_hits.to_csv(sys.argv[2]+'_summary.csv', index=False, header=list(sum_hits))
