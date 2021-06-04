import sys

## Replace all "?" characters with "N"

## Call this script in the command line with:
## python3 qton.py <sequencefile>

if len(sys.argv) != 2:
    print("qton.py takes exactly one argument.\nEnter: python3 qton.py <sequencefile>")
    quit()

fin = open(sys.argv[1])
data = fin.read()
fin.close()

data = data.replace('?', 'N')

fout = open(sys.argv[1], 'w')
fout.write(data)
fout.close()

