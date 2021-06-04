import re
import sys

def nuvselect(contig_data, contig_list):
    '''
    contig_data: string of all contigs
    contig_list: list of int(or ints as strings)
    returns string of selected contigs
    '''
    contigs = ''
    contig_list = list(map(str, contig_list))
    for num in contig_list:
        hit = re.compile(f'>sequence_{num}\|.+\|.*\n.*\n')
        findhit = hit.findall(contig_data)
        try:
            contigs+=findhit[0]
        except:
            print('Contig '+num+' not found')
    return contigs

def primeselect(contig_data, contig_list):
    '''
    contig_data: string of all contigs
    contig_list: list of int(or ints as strings)
    returns string of selected contigs
    '''
    contigs = ''
    contig_list = list(map(str, contig_list))
    for num in contig_list:
        hit = re.compile(f'(>Contig_{num}.+\n([A-Z]*\n)+)')
        findhit = hit.findall(contig_data)
        try:
            contigs+=findhit[0][0]
        except:
            print('Contig '+num+' not found')
    return contigs

# function to split .fa into header and sequence of all seqs
def readfa(filein):
    fin = open(filein)
    filein = fin.read()
    fin.close()
    filein = filein.split('>')[1:] #first char needs to be >
    seqs = []
    for f in filein:
        f = '>'+f
        f = f.split('\n')
        header = f[0]
        sequence = ''.join(f[1:])
        seqs.append([header, sequence])
    return seqs

def splitup(file):
    snivs = readfa(file)
    splitsnivs = []
    for s in snivs:
        if 'N' in s[1]:
            while s[1].count('NN') != 0:
                s[1] = s[1].replace('NN', 'N')
            
            conts = s[1].split('N')
            for n in range(len(conts)):
                newcont = [ s[0]+' '+str(n) , conts[n] ]
                splitsnivs.append(newcont)
        else:
            splitsnivs.append(s)

        fout = open(file, 'w')
        for s in splitsnivs:
            fout.write(s[0]+'\n')
            fout.write(s[1]+'\n')
        fout.close()

def write(filename, data):
    fout = open(filename, 'w')
    fout.write(data)
    fout.close()

def nuv(filename, contig_data, contig_list):
    write(filename, nuvselect(contig_data, contig_list))

def prime(filename, contig_data, contig_list):
    write(filename, primeselect(contig_data, contig_list))
    splitup(filename)

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Enter: python3 "+sys.argv[0]+" <contigs.fa> <contig_type> <output.fa> <list_of_contigs>")
        quit()
    
    fin = open(sys.argv[1])
    contig_data = fin.read()
    fin.close()
    contig_list = sys.argv[4].split(',')

    if sys.argv[2] == 'nuv': pull_contigs = nuv
    if sys.argv[2] == 'prime': pull_contigs = prime
    pull_contigs(sys.argv[3], contig_data, contig_list)
    