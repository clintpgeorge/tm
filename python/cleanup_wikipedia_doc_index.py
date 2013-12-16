#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import os 


'''
Checks whether the files in INPUT_FILE exist in 
the data folder and saves the existing files deatils 
into OUTPUT_FILE.  
'''

input_file = '/home/cgeorge/workspace/whales_tires/whales_tires.csv'
output_file = '/home/cgeorge/workspace/whales_tires/whales-tires.txt'
data_folder = '/home/cgeorge/workspace/whales_tires/docs'  

doc_id = 0
with open(input_file) as fp: 
    with open(output_file, 'w') as fw:
        print >>fw, 'docid,category,subject,docpath' # headers 
        for line_num, line in enumerate(fp):
            if line_num > 0:
                _, category, subject, _, _ = line.strip().split(';')
                doc_path = os.path.join(data_folder, subject + '.txt')
                if os.path.exists(doc_path):
                    print >>fw, '%d,%s,%s,%s' % (doc_id, category, subject, doc_path)
                    doc_id += 1

print doc_id, 'documents found.' 



