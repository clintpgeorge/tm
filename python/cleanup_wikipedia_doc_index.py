#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import os 
from shutil import copy2 

'''
Checks whether the files in INPUT_FILE exist in 
the data folder and saves the existing files details 
into OUTPUT_FILE.  
'''

input_file = '/home/cgeorge/workspace/whales_tires/whales_tires.csv'
output_file = '/home/cgeorge/workspace/whales_tires/whales-tires-cleaned.csv'
data_folder = '/home/cgeorge/workspace/whales_tires/docs'  
output_data_folder = '/home/cgeorge/workspace/whales_tires/docs-cleaned'

if not os.path.exists(output_data_folder):
    os.makedirs(output_data_folder)
    
doc_id = 0
subjects = []
with open(input_file) as fp: 
    with open(output_file, 'w') as fw:
        print >>fw, 'docid;category;subject' # headers 
        for line_num, line in enumerate(fp):
            if line_num > 0:
                _, category, subject, _, _ = line.strip().split(';')
                doc_path = os.path.join(data_folder, subject + '.txt')
                if os.path.exists(doc_path) and (subject not in subjects):  
                    print >>fw, '%d;%s;%s' % (doc_id, category, subject)
                    dest_file_name = os.path.join(output_data_folder, subject)
                    try:
                        copy2(doc_path, dest_file_name)
                    except (IOError, os.error) as why:
                        print doc_path, dest_file_name, str(why)
                    doc_id += 1
                    subjects.append(subject)

print doc_id, 'documents found.' 



