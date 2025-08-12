#!/usr/bin/env python3

import os
import sys
import datetime
import hashlib

class defaults:
    accession_file = 'Webin-accessions-2025-08-11T08_07_29.940+01_00.txt'
    study = 'PRJEB87316'
    instrument_model = ['NextSeq 2000', 'Illumina MiSeq']
    platform = 'ILLUMINA'
    library_name = ['']#needs to be extracted from the file names
    library_source = 'GENOMIC'
    library_selection = 'PCR'
    library_strategy = 'AMPLICON'
    library_layout = 'PAIRED'
    insert_size = '183'
    run_folder = ['/array-1/data/raw/240412_M00316_0204_000000000-L8G28/raw/', '/array-1/data/raw/230913_SN435_19_AAF37FLM5/raw/']
    destination = 'ENA_upload'

def calc_md5(datei):
    """function to calculate the md5 of a given file"""
    md5_hash = hashlib.md5()
    with open(datei,"rb") as f:
        for byte_block in iter(lambda: f.read(4096),b""):
            md5_hash.update(byte_block)
    return md5_hash.hexdigest()

dataset = {}
if not os.path.exists(defaults.accession_file):
    print(f"Input file {defaults.accession_file} does not exist.")
    sys.exit(1)
else:
    with open(defaults.accession_file, 'r') as infile:
        header = True
        for line in infile:
            if header:
                header = False
                continue
            else:
                line = line.strip().split('\t')
                if line[0] == 'SAMPLE':
                    dataset[line[2]] = [line[1]]
                    
if not os.path.exists(defaults.destination):
    os.makedirs(defaults.destination)

for f_enum, folder in enumerate(defaults.run_folder):
    if not os.path.exists(folder):
        print(f"Run folder {folder} does not exist.")
        sys.exit(1)
    for datei in os.listdir(folder):
        if datei.endswith('.fastq.gz'):
            ID = datei.split('_',2)[-1].rsplit('_',2)[0]
            if ID in dataset:
                print(F"Processing file {datei} for ID {ID}                                           ", end='\r', flush=True)
                dataset[ID].append((datei, calc_md5(os.path.join(folder, datei)), os.path.join(folder, datei), f_enum))

with open(f'{datetime.datetime.now().strftime("%y%m%d")}_ENA_fastq_submission.tsv','w') as outfile:
    outfile.write('\t'.join(['FileType', 'fastq', 'Read submission file type']+['']*24)+'\n')
    outfile.write('\t'.join(['study', 'sample', 'library_name', 'library_strategy', 'library_source', 'library_selection', 'library_layout', 'insert_size', 'library_construction_protocol', 'platform', 'instrument_model', 'forward_file_name', 'forward_file_md5', 'reverse_file_name', 'reverse_file_md5'] + ['']*8) + '\n')
    for ID in dataset:
        temp = sorted(dataset[ID][1:])
        entries = [defaults.study, dataset[ID][0], ID, defaults.library_strategy, defaults.library_source, defaults.library_selection, defaults.library_layout, defaults.insert_size, '', defaults.platform, defaults.instrument_model[temp[0][3]], temp[0][0], temp[0][1], temp[1][0], temp[1][1]]
        outfile.write('\t'.join(entries)+'\n')
        for i in range(2):
            print(F"Syncing file {temp[i][0]}                                                          ", end='\r', flush=True)
            os.system(f"rsync -a --checksum {temp[i][2]} {defaults.destination}/{temp[i][0]}")
            with open(f"{temp[i][0]}.md5", 'w') as md5_out:
                md5_out.write(temp[i][1])
        
        with open(F"{defaults.destination}/{ID}.manifest.txt", 'w') as f:
            f.write(F"STUDY\t{defaults.study}\n")
            f.write(F"SAMPLE\t{ dataset[ID][0]}\n")
            f.write(F"NAME\t{ID}\n")
            f.write(F"PLATFORM\t{defaults.platform}\n")
            f.write(F"INSTRUMENT\t{defaults.instrument_model[temp[0][3]]}\n")
            f.write(F"INSERT_SIZE\t{defaults.insert_size}\n")
            f.write(F"LIBRARY_NAME\t{ID}\n")
            f.write(F"LIBRARY_SOURCE\t{defaults.library_source}\n")
            f.write(F"LIBRARY_SELECTION\t{defaults.library_selection}\n")
            f.write(F"LIBRARY_STRATEGY\t{defaults.library_strategy}\n")
            f.write(F"FASTQ\t{temp[0][0]}\n")
            f.write(F"FASTQ\t{temp[1][0]}\n")
        
        

print()
print(f'lftp -u Webin-66927 ftp://webin.ebi.ac.uk -e "mirror -R {os.getcwd()}/ENA_upload; quit"')