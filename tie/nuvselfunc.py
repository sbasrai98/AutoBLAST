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





























