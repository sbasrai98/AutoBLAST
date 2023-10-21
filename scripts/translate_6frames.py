from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys

seqs = list(SeqIO.parse(sys.argv[1], 'fasta'))

translated = []
for s in seqs:
    for i in range(3):
        frame = s.seq[i:]
        if len(frame) % 3 != 0:
            frame = frame[:-(len(frame) % 3)]
        translated.append(SeqRecord(seq=frame.translate(), id=s.id+'_frame'+str(i+1), description=''))

        # reverse complement, frames 4, 5, 6
        frame = s.seq.reverse_complement()[i:]
        if len(frame) % 3 != 0:
            frame = frame[:-(len(frame) % 3)]
        translated.append(SeqRecord(seq=frame.translate(), id=s.id+'_frame'+str(i+4), description=''))

outname = sys.argv[1][::-1]
ext = outname.find('.')+1
outname = outname[ext:][::-1]+'_translated.fa'
SeqIO.write(translated, outname, 'fasta')
