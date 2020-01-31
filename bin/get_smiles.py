#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 18 15:55:16 2018

@author: megbert
"""
import sys
import csv
import argparse
import time
import pickle
import os
import wget
import logging
from threading import Thread
from queue import Queue
import pandas as pd
from rdkit import Chem


def get_SMILES(ligand, directory=os.getcwd()):       
    ligand = ligand.upper()

    # Download the LIGAND cif with wget
    url = 'https://files.rcsb.org/ligands/view/{}.cif'.format(ligand)
    file_path = os.path.join(directory, '{}.cif'.format(ligand))
    try:
        os.stat(file_path)
    except:
        try:
            wget.download(url, file_path)
        except:
            print('No infomation for lig :', ligand)            
    
    # Extract SMILES from the LIGAND cif file
    if os.path.exists(file_path):
        # Open the page & contents
        with open(file_path, 'r') as cif_file:
            cif_file = cif_file.readlines()
        
        # Grab the (OpenEye, Isomeric smiles for the ligand)
        smiles_options = []
        for line in cif_file:
            # SMILES
            if ligand in line and 'SMILES' in line:
                smiles_line = line.split() # smiles string is always the last entry in the list
                smiles_options.append(smiles_line[-1])
                smiles = smiles_line[-1].strip('"') # set smiles (gets overridden to last entry, OK)
            # MOLECULAR WEIGHT
            if line.startswith('_chem_comp.formula_weight'):
                if len(line.split()) > 1:
                    mw = line.split()[1]
                else:
                    mw = None
            # NAME
            if line.startswith('_chem_comp.name'):
                if len(line.split()) > 1:
                    name = line.split()[1].strip('"')
                else:
                    name = None

        
        # If the smiles selected is not isomeric, try to update it
        mol = Chem.MolFromSmiles(smiles); chiral = []        
        if mol != None:
            chiral = Chem.FindMolChiralCenters(mol)
        if len(chiral) == 0:
            for s in smiles_options:
                s = s.strip('"')
                mol = Chem.MolFromSmiles(s)
                if mol != None:
                    if mol.GetNumAtoms() > 4:
                        chiral = Chem.FindMolChiralCenters(mol)
                        if len(chiral) > 0:
                            smiles = s
                            break

        try:
            smiles
        except:
            smiles = 'n/a'
                
    # Return empty vectors if not able to download LIGAND cif file
    elif os.path.exists(file_path) == False:
        smiles = 'n/a'; mw = None; name = None
                
    return({'lig':ligand, 'name':name, 'smiles':smiles, 'mw':mw})


def get_SMILES_4_thread_queue2(q, result):
    while not q.empty():
        q_info = q.get()
        q_index = q_info[0]
        ligand = q_info[1].upper()
        directory = q_info[2]
        result_dict = get_SMILES(ligand, directory)        
        result[q_index].update(result_dict)
        q.task_done()        
    return True


def get_SMILES_queue(ligand_list, directory=os.getcwd()):
        
    # Setup thread queue
    q = Queue(maxsize=0)
    num_threads = min(50, len(ligand_list))
    results = [{} for x in ligand_list]
    for i in range(len(ligand_list)):
        q.put((i,ligand_list[i], directory))    
        
    # Starting worker threads on queue processing
    for i in range(num_threads):
        worker = Thread(target=get_SMILES_4_thread_queue2, args=(q,results))
        worker.setDaemon(True)
        worker.start()
    q.join()        

    return results
    
if __name__ == "__main__":
    
#    results = get_SMILES('0QH')    
#    print(results)
  
#    ligand_list = ['0Z4', '616', 'CIR', 'FRW', 'GIS', 'PEA', '2BL', 'TFK', 'XE', 'BAF', 'GLJ', '0QN', 'ETF', 'ACY', 'J54', 'BDK', 'TYJ', '0Z1', 'BAA', '4E4', 'SEI', 'IL0', '0QH', 'IBR', 'ASP', 'EOH', 'TPX', 'ACN', 'ACT', 'FPA', 'DBU', 'FR1', 'TRS', 'PHE', 'YNM', 'CNT', '0Z2', 'TPY', 'I3C', '0Z3', 'ICL', 'TFI', 'TSU', 'ICU', 'ARG', 'LYS', '0Z0', 'FLC', 'CL5', 'ALQ', 'SE4', 'BBL', 'VAI', '2Z5', '6NA', '0P2', 'CD', '681', '217', 'IOD', 'SUJ']
    ligand_list = ['LBM', 'XCP', 'LIU', 'F3Q', '1XJ', 'TPO', 'HOH', '1E9', '1XV', 'AJE', '398', '43B', 'CL', '1Y1', '2PE', 'ACE', 'MH8', 'PEG', 'SO4']        
    
    results = get_SMILES_queue(ligand_list)
    print(results)
    
