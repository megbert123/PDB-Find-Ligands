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
import wget
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

def create_dir(dir_name):
    try:
        os.stat(dir_name)
        exists = True
    except:
        os.mkdir(dir_name)
        exists = False
    return exists



class SimilarLigands:
    
    def __init__(self, pdb, chain, save_dir, mw_filter, similarity):
        self.pdb = pdb.upper()
        self.chain = chain
        self.save_dir = save_dir
        self.similarity = similarity
        self.mw_filter = mw_filter
        
        # Get Similar PDBs
        self.fasta_file = self.get_fasta()
        self.mmseq_file = self.run_mmseq()
        self.sim_list = self.pdb_get_similars()
        self.pdb_list = list(set([x.split('_')[0] for x in self.sim_list]))
        print(self.pdb_list)
                
        # Get Ligands for similar PDBs
        self.lig_dir = os.path.join(self.save_dir, 'ligands')
        create_dir(self.lig_dir)
        self.lig_dict = self.get_ligands_thread()
        self.ligand_list = list(set([x for y in [z['ligand_ids'] for z in self.lig_dict] for x in y]))

        # Get Smiles for ligands
        self.smiles_dir = os.path.join(self.save_dir, 'smiles')
        create_dir(self.smiles_dir)
        self.smiles_dict = self.get_smiles_thread()    

        # Combine info to get relevant ligands for alignment
        self.sim_ligs, self.pdb_ligs_dict = self.combine()
        print(pd.DataFrame.from_dict(self.pdb_ligs_dict))
        
        # Filter table to get a unique set of ligands
        self.unique_dict = self.filter_table()
        

    # Download FASTA seq from PDB
    def get_fasta(self):
        url = 'https://www.rcsb.org/fasta/entry/{}/download'.format(self.pdb)
        fasta_file = os.path.join(self.save_dir, '{}.fasta'.format(self.pdb))
        urllib.request.urlretrieve(url, fasta_file)

        with open(fasta_file, 'r') as f:
            f = f.readlines()
    
        for idx, line in enumerate(f):
            if line.startswith('>'):
                subject = line.split('|')[1]
                subject = subject.split(' ')[1]
                print(subject)
                if self.chain in subject:
                    header_idx = idx
                    break

        chain_fasta = f[header_idx]
        for line in f[header_idx+1:]:
            if line.startswith('>') == False:
                chain_fasta += line
            elif line.startswith('>') == True:
                break

        chain_fasta_file = os.path.join(self.save_dir, '{}.{}.fasta'.format(self.pdb, self.chain))
        with open(chain_fasta_file, 'w') as f_out:
            f_out.write(chain_fasta)
        return chain_fasta_file

    # Run MMseqs2 to find PDBs with similar sequences
    def run_mmseq(self):
       
        # Create mmseq PDB database (if not already in existence)
        pdb_database = os.path.join(os.path.split(self.save_dir)[0], 'PDB_database', 'pdb')
        db_exists = create_dir(os.path.split(pdb_database)[0])
        if db_exists == False:
            os.chdir(os.path.split(self.save_dir)[0])
            cmd = ['mmseqs', 'databases', 'PDB', pdb_database, 'tmp']
            check_call(cmd)

        # Search FASTA seq against PDB database
        os.chdir(self.save_dir)
        aln_res_file = '{}.{}.alnRes.m8'.format(self.pdb, self.chain)
        cmd = ['mmseqs', 'easy-search', self.fasta_file, pdb_database, aln_res_file, 'tmp']
        check_call(cmd)
        return aln_res_file


    # Extract similar PDB IDs from mmseqs alnRsn.m8 file
    def pdb_get_similars(self):
         
        # MMseqs2 alnRsn.m8 file: Query ID, PDB ID, Sequence Identity, Alignment Lengthm # Mismatches, 
        # Gaps, Query Domain Start, Query Domain End, Target Domain Start, Target Domain End, E-value, Bit Score
        with open(self.mmseq_file, 'r') as f:
            f = f.readlines()
        pdb_list = []
        for line in f:
            line = line.split('\t')
            if float(line[2]) >= self.similarity:
                pdb_list.append(line[1].upper())
        return pdb_list

    
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
            pdb, chain = entry.split('_')
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
    parser.add_argument('--sequence_similarity', required=False, default=0.9, help="PDBs with SS > XX*100 % will be considered, Enter # between 0 and 1")
    args = parser.parse_args()
    
    ## Make a directory to deposit results
    save_dir = os.path.join(args.save_dir, '{}_{}'.format(args.pdb, args.chain))
    create_dir(save_dir)

    # Find the similar ligands
    s = SimilarLigands(args.pdb, args.chain, save_dir, args.mw_filter, float(args.sequence_similarity))

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
    check_call(cmd)

    # Output pickle table (FILTERED CASES) & convert to CSV
    uni_pkl_file = os.path.join(save_dir, 'pdb_ligs_dict_UNIQUE.pkl')
    with open(uni_pkl_file, 'wb') as pkl_out:
        pickle.dump(s.unique_dict, pkl_out)
    cmd = ['python', PKL2CSV, '--', uni_pkl_file]
    check_call(cmd)

    # Make pymol session    
    cmd = ['pymol', '-qrc', MAKE_PSE, '--', '{}.{}'.format(args.pdb, args.chain), sl_file, save_dir]
    check_call(cmd)



