import re
import sys

## Call this script in the command line with:
## python3 digest.py <blastresultsfile> <sampleID>

if len(sys.argv) < 3:
    print("Not enough arguments.\nEnter: python3 digest.py <blastresultsfile> <sampleID>")
    quit()

# load blast hits for sample n
def loadhits(): 
    fin = open(sys.argv[1])
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
    virushits = list(filter(lambda x: 'Viruses' in x, tophits))
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
    summary.append('Sample '+sys.argv[2]+' BLAST summary\n')
    summary.append(',,Virus,Contigs,Hit Length,Identity\n')
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
        data.append(','.join([sys.argv[2], contig, virus, pid, alnlen, 
                              annotation, qlen, slen, accession])+'\n')
    lowpid = list(filter(lambda x: float(x.split(',')[3]) < 95, data))
    lowpid.sort(key = lambda x: float(x.split(',')[3]))
    lowpid = lowpid[:15]
    headings = ','.join(['Sample #','Contig #','Virus','Identity',
                          'Hit Length','Annotation','Query Length',
                          'Subject Length','Accession'])+'\n'
    summary = summarize(hitsdict)

    fout = open(sys.argv[2]+'_blast_summary.csv', 'w')
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



