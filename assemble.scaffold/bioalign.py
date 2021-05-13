from Bio import SeqIO
from Bio.Seq import Seq
import sys

if len(sys.argv) != 3:
    print('enter: python3 '+sys.argv[0]+' <reference.fa> <contigs.fa>') 
    quit()

contig_file = sys.argv[2]
reference = sys.argv[1]
contigs = list( SeqIO.parse(contig_file, 'fasta') )
ref = SeqIO.read(reference, 'fasta')

from Bio import Align
aligner = Align.PairwiseAligner()
aligner.mismatch_score = -2
aligner.gap_score = -3
aligner.left_gap_score = 0
aligner.right_gap_score = 0

# ref(Seq) + contigs(list of Seq) -> list of mappings
def map_seqs(ref, seqs):
    positions = []
    omitted = []
    for s in seqs:
        foraln = aligner.align(s.seq, ref.seq)[0]
        revaln = aligner.align(s.reverse_complement().seq, ref.seq)[0]
        if foraln.score > revaln.score:
            usealn = foraln
            direction = 'F'
        else:
            usealn = revaln
            direction = 'R'
        if len(usealn.aligned[1]) != 1: # should only be one aligned segment
            # FAULTY INDEL DETECTION...triggers for some other case too, not just indel
            omitted.append(s.description)
            continue
        ref_startidx = usealn.aligned[1][0][0]
        #stopidx = usealn.aligned[1][0][1] # never used..
        con_startidx = usealn.aligned[0][0][0]
        if usealn.score < 50: # score too low, should not be aligned
            # maybe include parameter: min score relative to contig size?
            omitted.append(s.description)
            continue
        if ref_startidx < con_startidx:
            print('%i nucleotides were trimmed from 5\' end of:' % con_startidx)
            print(s.description)
            s.description += ' (%i nt trimmed from 5\' end)' % con_startidx
        positions.append([s.description, ref_startidx, con_startidx, direction])
    
    print('%i contigs mapped to:' % len(positions))
    print(ref.description)
    print('%i contigs not mapped' % len(omitted))
    positions.sort(key=lambda x: x[1])
    return positions

# list of mappings + list of Seqs -> fasta MSA as list of Seqs
def build_msa(positions, seqs):
    msa_seqs = []
    for h in positions:
        sequence = list(filter(lambda x: x.description == h[0], seqs))[0]
        if h[3] == 'F':
            sequence = sequence.seq
        else:
            sequence = sequence.reverse_complement().seq
        seqid = h[0]
        line = '-'*h[1]
        line += sequence[h[2]:]
        msa_seqs.append( SeqIO.SeqRecord(id=seqid, seq=line, description=''))
    return msa_seqs

# list of Seqs(msa) -> Seq (consensus)
def consensus(msa):
    seqs = msa
    maxlen = max( list(map(lambda x: len(x.seq), seqs)) )
    # pad ends with gaps
    for s in seqs:
        pad = '-'*(maxlen - len(s.seq))
        s.seq += pad
    collapsed = ''
    for n in range(maxlen):
        posns = list(map(lambda x: x.seq[n], seqs))
        add = list(filter(lambda x: x != '-', posns))
        if add == []:
            collapsed += '-'
        else:
            # get list of highest count letters, choose first one
            highestcount = max( list(map(lambda x: add.count(x), add)) )
            highestchars = list(filter(lambda x: add.count(x) == highestcount, add))
            collapsed += highestchars[0]
    return SeqIO.SeqRecord(id='consensus', seq=Seq(collapsed), description='')

# Mapping same set of contigs to multiple refs? Will overwrite files..
posns = map_seqs(ref, contigs)
msa = build_msa(posns, contigs)
SeqIO.write(msa, sys.argv[2]+' msa.fa', 'fasta')

msacollapse = consensus(msa)
SeqIO.write(msacollapse, sys.argv[2]+' consensus.fa', 'fasta')
