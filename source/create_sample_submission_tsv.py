#!/usr/bin/env python3

import os
import sys


folders = ['/array-1/data/raw/240412_M00316_0204_000000000-L8G28/raw/', '/array-1/data/raw/230913_SN435_19_AAF37FLM5/raw/']
output_filename = '250808_ENA_Samples.tsv'

header = ['tax_id', 'scientific_name', 'sample_alias', 'sample_title', 'sample_description', 'collection date', 'geographic location (country and/or sea)', 'collected_by', 'serotype (required for a seropositive sample)', 'strain', 'sub_strain', 'genotype', 'culture_medium']

species = [('511693', 'BL21(AI)'), ('469008' ,'BL21(DE3)')]
species_name = 'Escherichia coli'

with open(output_filename, 'w') as outfile:
    outfile.write('\t'.join(header) + '\n')
    outfile.write('\t'.join(['#units']+['']*(len(header)-1)) + '\n')
    for f_enum, folder in enumerate(folders):
        for filename in os.listdir(folder):
            if filename.endswith('R1_001.fastq.gz'):
                sample_alias = filename.split('_', 2)[2].rsplit('_', 2)[0]
                
                outfile.write('\t'.join([species[f_enum][0], species_name, sample_alias] + ['']*5 + ['NA', 'BL21' , species[f_enum][1], '', '']) + '\n')