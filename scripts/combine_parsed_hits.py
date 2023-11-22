import os
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("hits_folder", help="directory containing .tsv files to combine")
#parser.add_argument("output", help="file name for output")
args = parser.parse_args()
args.hits_folder = args.hits_folder.rstrip('/')

samples = [file for file in os.listdir(args.hits_folder) if file.endswith('_parsed_hits.txt')]
combined = pd.read_table(f'{args.hits_folder}/{samples[0]}')
for s in samples[1:]:
    combined = pd.concat([combined, pd.read_table(f'{args.hits_folder}/{s}')])
combined.to_csv(f'{args.hits_folder}/combined_hits.txt', index=False, sep='\t')
