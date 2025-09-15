#! /usr/bin/python3
import json
from itertools import islice

lib = 'lib-ccpnmrV3_to_cyana.lib'
'''
seq = 'SilB41.seq'
prot = 'name.prot'
'''
def nomenclature(lib=lib, prot=None,seq=None, output_file= None):
    
    # creat a dictionaire with link amino acid and number in the sequence
    seq_dic = {} # empty dict
    
    with open(seq, 'r') as f:
        for line in f:
            if line.strip():
                amino_acid, res_num = line.strip().split() #splite number and aa
                seq_dic[int(res_num)] = amino_acid.upper() #Put the aa in upper letter
    
    #creat a dictionary with the file lib writ in json with {} --> go to dict
    with open (lib, 'r') as file:
        lib_dic = json.load(file) # read only str
        
   # Initialiser la liste pour accumuler les traductions
    atom_translet = []
    
    #open and read the .prot file
    with open (prot, 'r') as ff:
        items = list(islice(ff, 3))# Skip header line *3 == next (ff) *3
        for line in ff:
            if line.strip():
                parts = line.strip().split()
                atom_name = parts[-2] # take the atmos name
                seq_code = int(parts[-1]) # take the number of aa
                aa=[]
                aa = seq_dic[seq_code]
                if aa in ['ALA', 'ARG', 'GLY', 'HIS', 'TRP', 'VAL', 'ILE', 'LEU']:
                    try:
                        # Trying to access amino acid-specific translation
                        atom_translet.append(lib_dic[aa][atom_name])
                    except KeyError:
                        try:
                            # If unsuccessful, try general translation
                            atom_translet.append(lib_dic["GENERAL"][atom_name])
                        except KeyError:
                            # If unsuccessful, use atom name as default value
                            atom_translet.append(atom_name)
                else :
                    try:
                        #try the translation in general
                        atom_translet.append(lib_dic["GENERAL"][atom_name])
                    except KeyError:
                        # If unsuccessful, use atom name as default value
                        atom_translet.append(atom_name)
    return atom_translet
            
def file_transforme(atom_translet, prot=None, output_file= None):
    with open(prot, 'r') as f:
        # Conserver l'en-tête
        header = list(islice(f, 3))
        # Lire les lignes restantes
        lines = list(f)
    
    # Écrire dans le fichier de sortie
    with open(output_file, 'w') as f_out:
        # Écrire l'en-tête
        for line in header:
            f_out.write(line)
        # Modifier la troisième colonne
        for i, line in enumerate(lines):
            if line.strip():
                parts = line.strip().split()
                parts[3] = atom_translet[i]  # Remplacer la troisième colonne
                f_out.write("\t".join(parts) + "\n")
            else:
                f_out.write(line)
    
     

        
        
        