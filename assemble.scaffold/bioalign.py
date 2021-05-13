from Bio import SeqIO
from Bio.Seq import Seq
contig_file = 'split farmaTo contigs.fa'
reference = 'MT002973 ToBRFV.fasta'
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
    for s in seqs:
        foraln = aligner.align(s.seq, ref.seq)[0]
        revaln = aligner.align(s.reverse_complement().seq, ref.seq)[0]
        if foraln.score > revaln.score:
            usealn = foraln
            direction = 'F'
        else:
            usealn = revaln
            direction = 'R'
        if len(usealn.aligned) != 2: # problem...
            print('The following sequence was not included due to a possible indel:')
            print(s.description)
            continue
        startidx = usealn.aligned[1][0][0]
        stopidx = usealn.aligned[1][0][1]
        positions.append([s.description, startidx, stopidx, direction])
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
        line += sequence
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

posns = map_seqs(ref, contigs)
msa = build_msa(posns, contigs)
SeqIO.write(msa, 'tobrfv msa.fa', 'fasta')

msacollapse = consensus(msa)
SeqIO.write(msacollapse, 'tobrfv consensus.fa', 'fasta')
