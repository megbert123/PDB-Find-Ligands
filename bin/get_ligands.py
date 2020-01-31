#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 18 15:55:16 2018

@author: megbert
"""

import csv
import argparse
import time
import pickle
import os
import sys
import pandas as pd
import wget
from queue import Queue
from threading import Thread

   
def get_LIGANDS(pdb, directory=os.getcwd()):       
    pdb = pdb[:4].lower()

    # Download the PDB file with wget
    url = 'https://files.rcsb.org/view/{}.pdb'.format(pdb)
    file_path = os.path.join(directory, '{}.pdb'.format(pdb))
    try:
        os.stat(file_path)
    except:
        try:
            wget.download(url, file_path)
        except:
            print('No ligands for pdb id :', pdb)            
    
    # Extract ligand IDs & chains from the PDB file
    if os.path.exists(file_path):
        # Open the page & contents
        with open(file_path, 'r') as pdb1_file:
            pdb1_file = pdb1_file.readlines()
        
        # Grab the ligand names from the file
        ligands = []
        for line in pdb1_file:
            if line.startswith('FORMUL'):
                lig = line.split()[2]
                ligands.append(lig)
                
        # Grab the corresponding chains for each ligand
        ligand_chains = {}
        for line in pdb1_file:
            if line.startswith('HETATM'):
                ligand = line[17:20].strip()
                chain = line[21:22]
                if ligand not in ligand_chains:
                    ligand_chains.update({ligand:[chain]})
                elif ligand in ligand_chains and chain not in ligand_chains[ligand]:
                    ligand_chains[ligand].append(chain)
                    
    # If not able to download PDB file, return empty vectors
    elif os.path.exists(file_path) == False:
        ligands  = []
        ligand_chains = [] 
        
    return({'pdb':pdb, 'ligand_ids':ligands, 'ligand_chains':ligand_chains})


def get_LIGANDS_4_thread_queue(q, result):  
    while not q.empty():
        q_info = q.get()
        q_index = q_info[0]
        pdb = q_info[1][:4].lower()
        directory = q_info[2]
        result_dict = get_LIGANDS(pdb, directory)
        result[q_index].update(result_dict)
        q.task_done()        
    return True    
    
def get_LIGANDS_queue(pdb_list, directory=os.getcwd()):
        
    # Setup thread queue
    q = Queue(maxsize=0)
    num_threads = min(50, len(pdb_list))
    results = [{} for x in pdb_list]
    for i in range(len(pdb_list)):
        q.put((i,pdb_list[i], directory))    
        
    # Starting worker threads on queue processing
    for i in range(num_threads):
        worker = Thread(target=get_LIGANDS_4_thread_queue, args=(q,results))
        worker.setDaemon(True)    
        worker.start()
    q.join()
        
    return results    
    
    
    
if __name__ == "__main__":
    
    results = get_LIGANDS('1ELA')    
    print(results)

#    pdb_list = ['1MCV', '1JIM', '2EST', '3ODF', '2BDB', '1ELC', '2BD5', '1QR3', '1FLE', '1L1G', '1EAI', '1LKB', '1HAY', '2V0B', '1C1M', '1EST', '1OKX', '1UVP', '1EAT', '3MO3', '1NES', '6Q8S', '3MU4', '3HGP', '2G4U', '1ELA', '4GVU', '2D26', '7EST', '1ELE', '2DE8', '2FOE', '1HAZ', '1E37', '5AVD', '3MU0', '1H9L', '1FZZ', '2FOB', '2BLQ', '1GWA', '1EAU', '6EST', '2BD2', '1MMJ', '2BD4', '2BB4', '2BD7', '1HAX', '3EST', '2FOH', '9EST', '4YM9', '3HGN', '4EST', '1QNJ', '1BMA', '1E38', '1HV7', '2A7C', '3MU1', '3MNX', '1E35', '2IOT', '2FO9', '1UVO', '1EAS', '2FOD', '2A7J', '2DE9', '2BLO', '3MU8', '1ELD', '2V35', '3MTY', '5EST', '2CV3', '2G4T', '1ELB', '1LKA', '3MNS', '2BD9', '1QGF', '2OQU', '2FOA', '3MOC', '1BTU', '1E34', '1QIX', '1UO6', '1ESA', '3MNC', '2BDA', '3E3T', '1HB0', '1GVK', '2BD3', '3MNB', '8EST', '3MU5', '1ELG', '3MO6', '3ODD', '3UOU', '2H1U', '1LVY', '2BD8', '2BDC', '2FOF', '1L0Z', '1ELF', '1INC', '3MO9', '2FOG', '2FOC', '1ESB', '1B0E', '1E36']#    results = get_LIGANDS_queue(pdb_list)
#    print(results)

