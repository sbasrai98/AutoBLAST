import re
import sys

## Call this script in the command line with:
## python3 select.py <sequencefile> <list> <outputfile>

if len(sys.argv) < 4:
    print("Not enough arguments.\nEnter: python3 nuvselect.py <sequencefile> <list> <outputfile>")
    quit()

fin = open(sys.argv[1])
data = fin.read()
fin.close()

writethese = sys.argv[2].split(',')
try:
    writethese = list(map(int, writethese))
except:
    print('Incorrect list format. Separate each number with ",".\nExample: 3,48,94,213,4089')
writethese = list(map(str, writethese))

contigs = ''
for num in writethese:
    hit = re.compile(f'(>Contig_{num}.+\n([A-Z]*\n)+)')
    findhit = hit.findall(data)
    try:
        contigs+=findhit[0][0]
    except:
        print('Contig '+num+' not found')

 #   hit = re.compile(f'>sequence_{num}\|.+\|.*\n.*\n')
 #   findhit = hit.findall(data)[0]
 #   contigs+=findhit

fout = open(sys.argv[3], 'w')
fout.write(contigs)
fout.close()





























