#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 13:52:33 2020

@author: megbert
"""
"""The purpose of this script is to find all the ligands that the PDB of interest,
    and all PDBs within a sequence similarilty cuttoff, load them in a pymol file 
    for easy visualization. This script will also output an excel file with smiles
    information, etc.

"""

import argparse
import os
import sys
import json
import time
import pickle
import pandas as pd
import numpy as np
import multiprocessing as mp
from multiprocessing import Pool
from subprocess import call
from subprocess import check_call
from subprocess import check_output
import urllib.request
import xml.etree.ElementTree as ET
sys.path.insert(0, sys.path[0])
from get_ligands import get_LIGANDS_queue
from get_smiles import get_SMILES_queue


# SCRIPTS
MAKE_PSE = os.path.join(sys.path[0], 'align_ligs_pse.py')
PKL2CSV = os.path.join(sys.path[0], 'pkl2csv.py')

class SimilarLigands:
    
    def __init__(self, pdb, chain, save_dir, mw_filter, similarity):
        self.pdb = pdb
        self.chain = chain
        self.save_dir = save_dir
        self.similarity = similarity
        self.mw_filter = mw_filter
        
        # Get Similar PDBs
        self.sim_list = self.pdb_get_similars()
        self.pdb_list = list(set([x.split('.')[0] for x in self.sim_list]))
        print(self.pdb_list)
                
        # Get Ligands for similar PDBs
        self.lig_dir = os.path.join(self.save_dir, 'ligands')
        try:
            os.stat(self.lig_dir)
        except:
            os.mkdir(self.lig_dir)
        self.lig_dict = self.get_ligands_thread()
        self.ligand_list = list(set([x for y in [z['ligand_ids'] for z in self.lig_dict] for x in y]))
        print(self.ligand_list)

        # Get Smiles for ligands
        self.smiles_dir = os.path.join(self.save_dir, 'smiles')
        try:
            os.stat(self.smiles_dir)
        except:
            os.mkdir(self.smiles_dir)
        self.smiles_dict = self.get_smiles_thread()    
        
        # Combine info to get relevant ligands for alignment
        self.sim_ligs, self.pdb_ligs_dict = self.combine()
        print(pd.DataFrame.from_dict(self.pdb_ligs_dict))
        
        # Filter table to get a unique set of ligands
        self.unique_dict = self.filter_table()
        

    # Get similar PDBS
    def pdb_get_similars(self):
        pdb_chain = self.pdb + '.' + self.chain
        url = 'http://www.rcsb.org/pdb/rest/sequenceCluster?cluster=%d&structureId=%s' % (self.similarity, pdb_chain)
        req = urllib.request.Request(url)
        f = urllib.request.urlopen(req)
        result = f.read()
        root = ET.fromstring(result)
        sim_list = []    
        for pdbchain in root.iter('pdbChain'):
            sim_list.append(pdbchain.attrib['name'])    
        return sim_list
    
    # Get ligands for similar pdbs
    def get_ligands_thread(self):
        lig_dict = get_LIGANDS_queue(self.pdb_list, self.lig_dir)
        return lig_dict
    
    # Get smiles for ligands
    def get_smiles_thread(self):
        smiles_dict = get_SMILES_queue(self.ligand_list, self.smiles_dir)
        return smiles_dict

    # Combine ligands and smiles information        
    def combine(self):
        
        sim_ligs = []
        pdb_ligs_dict = []
        
        for entry in self.sim_list:
            pdb = entry.split('.')[0]
            chain = entry.split('.')[1]
            pdb_lig_info = [x for x in self.lig_dict if x['pdb'].upper() == pdb.upper()][0]
            for l in pdb_lig_info['ligand_ids']:
                if chain in pdb_lig_info['ligand_chains'][l]:
                    lig_info = [x for x in self.smiles_dict if x['lig'] == l]
                    if len(lig_info) > 0:
                        lig_info = lig_info[0]
                        if float(lig_info['mw']) > self.mw_filter:
                            sim_ligs.append('{}.{}.{}'.format(pdb, chain, l))                            
                            temp = {'pdb':pdb, 'chain':chain}
                            temp.update(lig_info)
                            pdb_ligs_dict.append(temp)
                            
        # sort alpha-numerically for easy searching
        sim_ligs = sorted(sim_ligs)
        pdb_ligs_dict = sorted(pdb_ligs_dict, key = lambda i: i['pdb']) 
                            
        return sim_ligs, pdb_ligs_dict       

    def filter_table(self):
#        unique_list = ['{}_{}'.format(x['lig'], x['mw']) for x in self.pdb_ligs_dict]
        unique_dict = []; unique_list = []
        for entry in self.pdb_ligs_dict:
            uid = '{}_{}'.format(entry['lig'], entry['mw'])
            pdb_chain = '{}_{}'.format(entry['pdb'], entry['chain'])
            
            if uid not in unique_list:
                # add entry to unique_dict
                unique_dict.append(entry)
                unique_list.append(uid)

            elif uid in unique_list:
                # update the existing entry in unique_dict
                for x in unique_dict:
                    if '{}_{}'.format(x['lig'], x['mw']) == uid:
                        if 'additional_pdbs' in x:
                            x['additional_pdbs'].append(pdb_chain)
                        elif 'additional_pdbs' not in x:
                            x.update({'additional_pdbs':[pdb_chain]})
                                            
        return unique_dict                
                




if __name__ == "__main__":
        
    parser = argparse.ArgumentParser()
    parser.add_argument('pdb')
    parser.add_argument('chain')
    parser.add_argument('--save_dir', required=False, default = os.getcwd())
    parser.add_argument('--mw_filter', required=False, default=100, help="Ligands with MW < XX Da will be filtered out")
    parser.add_argument('--sequence_similarity', required=False, default=90, help="PDBs with SS > XX % will be considered")
    args = parser.parse_args()
    
    ## Make a directory to deposit results
    save_dir = os.path.join(args.save_dir, '{}_{}'.format(args.pdb, args.chain))
    try:
        os.stat(save_dir)
    except:
        os.mkdir(save_dir)
    
    # Find the similar ligands
    s = SimilarLigands(args.pdb, args.chain, save_dir, args.mw_filter, args.sequence_similarity)

    sl_file = os.path.join(save_dir, 'sim_pdb_chain_ligs.txt')
    with open(sl_file, 'w') as sl_out:
        for entry in s.sim_ligs:
            sl_out.write(entry)
            sl_out.write('\n')
            
    # Output pickle table (ALL CASES) & convert to CSV
    pkl_file = os.path.join(save_dir, 'pdb_ligs_dict.pkl')
    with open(pkl_file, 'wb') as pkl_out:
        pickle.dump(s.pdb_ligs_dict, pkl_out)
    cmd = ['python', PKL2CSV, '--', pkl_file]
    print(cmd)
    check_call(cmd)

    # Output pickle table (FILTERED CASES) & convert to CSV
    uni_pkl_file = os.path.join(save_dir, 'pdb_ligs_dict_UNIQUE.pkl')
    with open(uni_pkl_file, 'wb') as pkl_out:
        pickle.dump(s.unique_dict, pkl_out)
    cmd = ['python', PKL2CSV, '--', uni_pkl_file]
    print(cmd)
    check_call(cmd)

        
    # Make pymol session    
    cmd = ['pymol', '-qrc', MAKE_PSE, '--', '{}.{}'.format(args.pdb, args.chain), sl_file, save_dir]
    print(" ".join(cmd))
    check_call(cmd)



