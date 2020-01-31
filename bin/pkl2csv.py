#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 17 13:49:29 2019

@author: megbert
"""

## this file converts a pickle file to a csv file
## 1 input: pickle_file.pkl
## the csv file is save in the same location as the pkl file,
## with the same name - but the csv extension pickle_file.csv

import os
import argparse
import pickle
import csv

parser = argparse.ArgumentParser()
parser.add_argument('pkl_file')
args = parser.parse_args()


if __name__ == "__main__":
    
    # Open pickle file
    with open(args.pkl_file, 'rb') as pkl_file:
        pkl_file = pickle.load(pkl_file)
        
    # Grab all possible dict keys in the file
    header = []
    for line in pkl_file:
        line_h = line.keys()
        for h in line_h:
            if h not in header:
                header.append(h)
    
    # Write to a csv file
    with open('{}.csv'.format(args.pkl_file[:-4]), 'w') as csv_file:
        csv_out = csv.writer(csv_file, delimiter=',', quotechar='"')
        csv_out.writerow(header)
        for line in pkl_file:
            temp = []
            for h in header:
                if h in line:
                    value = str(line[h])
                    if ',' in value:
                        value = '"{}"'.format(value)
                    temp.append(value)
                else:
                    temp.append('NA')
            csv_out.writerow(temp)
            
            
        
        
        
        
