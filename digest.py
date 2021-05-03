import re
import sys
import argparse

parser = argparse.ArgumentParser(description='Extract and summarize BLAST results based on some filter')
parser.add_argument('blast_results', help='file containing BLAST results')
parser.add_argument('sample_id', help='identifier for sample')
parser.add_argument('-filter', default='Viruses', help="string used to filter BLAST hits (default = 'Viruses')" )
args = parser.parse_args()

# load blast hits for sample n
def loadhits(): 
    fin = open(args.blast_results)
    blasthits = list(map(lambda x: x.split('\t'), fin.readlines()))
    fin.close()
    return blasthits

# trim to top hits only
def tophits(blasthits):
    top = []
    logged = []
    for h in blasthits:
        if not(h[0] in logged):
            top.append(h)
            logged.append(h[0])
    return top

# trim to viruses, only keep relevant info, sort by virus into dictionary
def virustrim(tophits):
    virushits = list(filter(lambda x: args.filter in str(x), tophits))
    short = list(map(lambda x: [x[0], x[1], x[2], x[3], x[14], x[11], x[12]], virushits))
    virkeys = list(map(lambda x: x[4], short))
    vhitsdict = dict.fromkeys(list(dict.fromkeys(virkeys)))
    for key in virkeys:
        vhitsdict[key] = list(filter(lambda x: x[4] == key, short))
    return vhitsdict

# sorts hits by contig count, then min length (for summary)
def summaryorder(hit):
    hit = hit.split(',')
    lower = re.compile('[0-9]+')
    hlow = int(lower.findall(hit[4])[0])
    return int(hit[3])*-1000000 - hlow

# summarize contigs data
def summarize(hitsdict):
    summary = [] 
    summary.append('\n')
    summary.append('Sample '+args.sample_id+' BLAST summary\n')
    summary.append(',,'+args.filter+',Contigs,Hit Length,Identity\n')
    hitlist = []
    for k in hitsdict:
        #virus = k
        contigs = str(len(hitsdict[k]))
        hitlens = list(map(lambda x: int(x[3]), hitsdict[k]))
        hlrange = str(min(hitlens))+' - '+str(max(hitlens))
        pidents = list(map(lambda x: float(x[2]), hitsdict[k]))
        idrange = str(min(pidents))+' - '+str(max(pidents)) 
        hitlist.append(',,'+','.join([k, contigs, hlrange, idrange])+'\n')
    hitlist.sort(key=summaryorder)
    summary.extend(hitlist)
    summary.append('\n')   
    return summary

# create output file
def output(hitsdict):
    data = []
    lst = []
    for k in hitsdict:
        lst.extend(hitsdict[k])

    contig_num = re.compile('_[0-9]+')
    #annotate = re.compile(f'\|{sample}\|.*')
    access = re.compile('\|[^a-z]+\.[0-9]+')
    for hit in lst:
        contig = contig_num.findall(hit[0])[0][1:]
        virus = hit[4]
        pid = hit[2]
        alnlen = hit[3]
        annotation = 'N/A'
        qlen = hit[5]
        slen = hit[6]
        accession = access.findall(hit[1])[0][1:]
        data.append(','.join([args.sample_id, contig, virus, pid, alnlen, 
                              annotation, qlen, slen, accession])+'\n')
    lowpid = list(filter(lambda x: float(x.split(',')[3]) < 95, data))
    lowpid.sort(key = lambda x: float(x.split(',')[3]))
    lowpid = lowpid[:15]
    headings = ','.join(['Sample #','Contig #',args.filter,'Identity',
                          'Hit Length','Annotation','Query Length',
                          'Subject Length','Accession'])+'\n'
    summary = summarize(hitsdict)

    fout = open('Sample '+args.sample_id+' '+args.filter+' BLAST summary.csv', 'w')
    for s in summary:
        fout.write(s)
    fout.write('Lowest identity top hits:\n')
    fout.write(headings)
    for l in lowpid:
        fout.write(l)
    fout.write('\n')
    fout.write('All hits:\n')
    fout.write(headings)
    for d in data:
        fout.write(d)
    fout.close()

output(virustrim(tophits(loadhits())))



