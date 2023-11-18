# %%
import os
import argparse
import pandas as pd
import xml.etree.ElementTree as ET

# %%
parser = argparse.ArgumentParser()
parser.add_argument("sample_list", help="file containing SRA accessions")

#with open(args.sample_list) as file_in: # 
with open('sra_samples1.txt') as file_in:
    samples = [x.rstrip() for x in file_in.readlines()]

os.makedirs('xml', exist_ok=True) # make xml directory if it doesn't exist

# for s in samples:
    # os.system(f"esearch -db sra -query '{s}[accession]' | efetch -format native -mode xml > xml/{s}_meta.xml")

# %%

# from root, get desired fields

metadata_paths = {}

# root is EXPERIMENT_PACKAGE_SET. path starts from next node.
metadata_paths['Experiment Title'] = ['EXPERIMENT_PACKAGE','EXPERIMENT','TITLE']
metadata_paths['Design Description'] = ['EXPERIMENT_PACKAGE','EXPERIMENT','DESIGN','DESIGN_DESCRIPTION']
metadata_paths['Library Name'] = ['EXPERIMENT_PACKAGE','EXPERIMENT','DESIGN','LIBRARY_DESCRIPTOR','LIBRARY_NAME']
metadata_paths['Library Strategy'] = ['EXPERIMENT_PACKAGE','EXPERIMENT','DESIGN','LIBRARY_DESCRIPTOR','LIBRARY_STRATEGY']
metadata_paths['Library Source'] = ['EXPERIMENT_PACKAGE','EXPERIMENT','DESIGN','LIBRARY_DESCRIPTOR','LIBRARY_SOURCE']
metadata_paths['Library Selection'] = ['EXPERIMENT_PACKAGE','EXPERIMENT','DESIGN','LIBRARY_DESCRIPTOR','LIBRARY_SELECTION']
metadata_paths['Sequencing Platform'] = ['EXPERIMENT_PACKAGE','EXPERIMENT','PLATFORM']
metadata_paths['Organization Name'] = ['EXPERIMENT_PACKAGE','Organization','Name']
metadata_paths['Organization Country'] = ['EXPERIMENT_PACKAGE','Organization','Address','Country']
metadata_paths['Organization City'] = ['EXPERIMENT_PACKAGE','Organization','Address','City']
metadata_paths['Contact Country'] = ['EXPERIMENT_PACKAGE','Organization','Contact','Address','Country']
metadata_paths['Contact City'] = ['EXPERIMENT_PACKAGE','Organization','Contact','Address','City']
metadata_paths['Study Title'] = ['EXPERIMENT_PACKAGE','STUDY','DESCRIPTOR','STUDY_TITLE']
metadata_paths['Study Abstract'] = ['EXPERIMENT_PACKAGE','STUDY','DESCRIPTOR','STUDY_ABSTRACT']
metadata_paths['Scientific Name'] = ['EXPERIMENT_PACKAGE','SAMPLE','SAMPLE_NAME','SCIENTIFIC_NAME']
metadata_paths['Taxon ID'] = ['EXPERIMENT_PACKAGE','SAMPLE','SAMPLE_NAME','TAXON_ID']
metadata_paths['Sample Location'] = []
metadata_paths['Sample Attributes'] = ['EXPERIMENT_PACKAGE','SAMPLE','SAMPLE_ATTRIBUTES'] # TODO: sample attributes, will need to use findall and create key/value pairs for each

sra_metadata = pd.DataFrame(columns=metadata_paths.keys())

for s in samples: #  SRR17589118  SRR25318402

    tree = ET.parse(f'xml/{s}_meta.xml')
    root = tree.getroot()

    metadata = dict.fromkeys(metadata_paths.keys())
    for k in metadata:
        try:
            current = root
            for tag in metadata_paths[k]:
                current = current.find(tag)
            if k == 'Sequencing Platform':
                metadata[k] = current[0][0].text
            elif k == 'Sample Attributes':
                attributes = ''
                for child in current:
                    if child[0].text == 'geo_loc_name':
                        metadata['Sample Location'] = child[1].text
                    else:
                        attributes += f'{child[0].text}:{child[1].text}, '
                attributes = attributes.rstrip(', ')
                metadata[k] = attributes
            else:
                metadata[k] = current.text
        except AttributeError:
            pass
    sra_metadata.loc[s] = metadata.values()

sra_metadata.index.rename('SRA Accession', inplace=True)
sra_metadata.to_csv('sra_metadata3.csv')