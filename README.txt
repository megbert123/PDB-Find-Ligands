# PDB-Find-Ligands
This project provides code to search for all ligands in the PDB that bind a given protein of interest.

$ python find_ligands_wrapper.py -- <PDB_ID> <CHAIN> <--save_dir> <--mw_filter> <--sequence_similarity>
    INPUTS:
      PDB_ID    == PDB_ID of the protein of interest (i.e. 1ELA)
      CHAIN     == Chain of PDB_ID that you want considered (similar proteins are found by searching the sequences of the PDB_ID / Chain)
      SAVE_DIR  == Directory where you want the results saved (optional)
      MW_FILTER == Minimum MW of returned ligands (default = 100 Da), because most ions / crystal additives have MW < 100 Da
      SEQ_SIM   == Sequence similarity allowed for searching for similar PDB_ID (i.e. 30 for 30% SS)
    OUTPUTS:
      - makes a folder <PDB_CHAIN> on the save directory to store all output files
      - makes a <ligands> and <smiles> folder to store ligand / smiles information download from the PDB
      - makes an <aligned_ligands> folder to store ligand pdb files with coordinates reflecting alignment to reference PDB_ID
      - PYMOL SESSION - all_ligands.pse, contains all ligands identified, aligned to reference PDB
      - pdb_ligs_dict.pkl (.csv) table containing PDB_ID / CHAIN / LIGAND_ID / LIGAND_NAME / SMILES / MW / ADDITIONAL_PDBs
      - pdb_ligs_dict_UNIQUE.pkl (.csv) table containing same as above, but with only one row per ligand ID
                                                                                                           
How it works
    $ find_ligands_wrapper.py
        1.) Find Similar PDBs (pdb_get_similars) using the built-in PDB function, searching by PDB_ID/CHAIN
        2.) Find Ligands associated with PDB_ID/CHAIN returned from (1) using wget from PDB website <PDB_ID>.pdb [get_ligamds.py]
        3.) Find SMILES associated with ligands pulled in (2) using wget from PDB website <LIG_ID>.cif, *see note* [get_smiles.py]
        4.) Combine all information (PDB_ID / CHAIN / LIG_ID / SMILES) and filter out ligands with MW < MW_FILTER
        5.) Output data into pdb_ligs_dict_(UNIQUE).pkl / .csv
        6.) Create pymol session 'all_ligands.pse' by aligning PDB_ID (and ligand) to reference PDB. [align_ligs_pse.py]
            
 
Isomeric Smiles Note: This script will return isomeric smiles for the ligand, whenever possible. If no isomeric smiles are in the .cif file, then canoncial smiles will be returned instead. Interestingly, sometimes (as in the case of 1BMA_A ligand 0QH) the isomeric smiles displayed on the PDB webpage (https://www.rcsb.org/ligand/0QH) is differnet from any of the smiles in the cif file (https://files.rcsb.org/ligands/view/0QH.cif'). There is one set of isomeric smiles in the cif file, it is labeled SMILES_CANONICAL CACTVS, but the stereochemistry is different from the isomeric smiles on the webpage. [i.e. the cif isomeric smiles show the rings in the cis formation, whereas the isomeric smiles on the web page show the rings in the trans formation - if you check the actual structure of the ligand in 1BMA_A_0QH, the rings are actually probably somewhere between cis and trans.]

