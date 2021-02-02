#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 16:24:11 2020

@author: megbert
"""

import argparse
from pymol import cmd
import os
import json


parser = argparse.ArgumentParser()
parser.add_argument('reference_pdb_chain')
parser.add_argument('pdb_chain_lig_file')
parser.add_argument('save_dir')
args = parser.parse_args()

aligned_lig_dir = os.path.join(args.save_dir, 'aligned_ligands')
try:
    os.stat(aligned_lig_dir)
except:
    os.mkdir(aligned_lig_dir)


cmd.fetch(args.reference_pdb_chain)
cmd.show_as('cartoon')

with open(args.pdb_chain_lig_file, 'r') as align_ligs:
    align_ligs = align_ligs.readlines()
    align_ligs = sorted(align_ligs)

print(align_ligs)

ligand_list = []
for pcl in align_ligs:
    print(pcl)
    pcl = pcl.strip().split('.')
    pdb_chain = '{}.{}'.format(pcl[0], pcl[1])
    ligand = pcl[2]    
    ligand_list.append(ligand)
    
    # load PDB into pymol, align to reference protein and save ligand
    cmd.fetch(pdb_chain)
    cmd.align(pdb_chain, args.reference_pdb_chain)
    pcl_name = '{}_{}_{}'.format(pcl[0], pcl[1], pcl[2])
    cmd.create(pcl_name, '{} and het and resn {}'.format(pdb_chain, ligand))
    cmd.save(os.path.join(aligned_lig_dir, pcl_name+'.pdb'), pcl_name)
    cmd.show_as('sticks', pcl_name)

    # delete protein if not reference protein
    if pdb_chain != args.reference_pdb_chain:
        cmd.delete(pdb_chain)

# Group by ligands for easy analysis
ligand_list = sorted(list(set(ligand_list)))
for lig in ligand_list:
    cmd.group(lig, '*_{}'.format(lig))
        
        
cmd.save(os.path.join(args.save_dir, 'all_ligands.pse'))
cmd.quit()
    
    
