#! /usr/bin/python3

import pandas as pd
import numpy as np
import monentlature
import json

'''
#information
attribution_prot = '/home/biosys/Bureau/python_cy/for_cyana_ccpnm/Bform/attrib_Bform_full.csv' # need to be change

sequence_code_dict = 'proline.txt'
'''

def process_protein_data(attribution_prot, sequence_code_dict,output_file):
    # read the file CSV
    data = pd.read_csv(attribution_prot)
    
    with open(sequence_code_dict, 'r') as fichier:
        proline = json.load(fichier) 

    #Verify if the colone 'Total\nPeak Count'exist
    if 'Total\nPeak Count' in data.columns:
        #Creat a new file without any zero or NaN in colomne Total\nPeak Count
        data_zero = data[(data['Total\nPeak Count'] !=0) & (pd.notna(data['Total\nPeak Count']))]
        
    else :
        print ('There is no "Total Peak Count" colomns in the file')
        
    # Delete several columns
    columns_to_delete = ['uniqueId','NmrAtom', 'ResidueType', 'Total\nPeak Count']
    #Creates a new list containing only the column names from columns_to_delete
    data_columns= data_zero.drop(columns=[col for col in columns_to_delete if col in data_zero.columns])
    #change the numerotation of the proline
    data_columns['SequenceCode']=data_columns['SequenceCode'].replace(proline)
    #remove all the attribution -1(duplicate)
    data_columns= data_columns[~data_columns['SequenceCode'].astype(str).str.contains('-1')]

    # Creat a new column with crescent number, begining by 1
    data_columns.insert(0, 'Index', range(1, len(data_columns) + 1))

    # Re order the columns
    colonnes_ordre = ['Index', 'Value\n(ppm)', 'Value Error\n(ppm)', 'AtomName', 'SequenceCode']
    data_columns = data_columns[colonnes_ordre]

    # Exchange all the None in 0
    data_columns.replace("None", 0, inplace=True)

    #Save the file
    data_columns.to_csv(output_file +'/name.prot', sep='\t', index=False)

    #converte file
    trad = monentlature.nomenclature(seq=output_file + '/prot.seq',prot= output_file +'/name.prot',output_file=output_file)
    monentlature.file_transforme(trad,  prot= output_file +'/name.prot',output_file = output_file+ '/attrib_cya.prot')
           
