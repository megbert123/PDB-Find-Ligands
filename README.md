# PDB-Find-Ligands
This project provides code to search for all ligands in the PDB that bind your protein of interest.



* isomeric smiles note

This script will return isomeric smiles for the ligand, whenever possible. If no isomeric smiles are in the .cif file, then canoncial smiles will be returned instead. Interestingly, sometimes (as in the case of 1BMA_A ligand 0QH) the isomeric smiles displayed on the PDB webpage (https://www.rcsb.org/ligand/0QH) is differnet from any of the smiles in the cif file (https://files.rcsb.org/ligands/view/0QH.cif'). There is one set of isomeric smiles in the cif file, it is labeled SMILES_CANONICAL CACTVS, but the stereochemistry is different from the isomeric smiles on the webpage. [i.e. the cif isomeric smiles show the rings in the cis formation, whereas the isomeric smiles on the web page show the rings in the trans formation - if you check the actual structure of the ligand in 1BMA_A_0QH, the rings are actually probably somewhere between cis and trans.]

