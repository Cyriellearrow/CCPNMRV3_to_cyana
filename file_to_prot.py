#!python3

import pandas as pd
import numpy as np

#information
attribution_prot = '/home/biosys/Bureau/python_cy/for_cyana_ccpnm/Bform/attrib_Bform_full.csv' # need to be change
sequence = '/home/biosys/Bureau/python_cy/for_cyana_ccpnm/Bform/sequencecode_clean.csv' #  need to be change
sequence_code_dict = {
    "45-1": 44,
    "55-1": 54,
    "58-1": 57,
    "65-1": 64,
    "69-1": 68,
    "77-1": 76
}  # here put the proline changement 


# read the file CSV
data = pd.read_csv(attribution_prot)

#Verify if the colone 'Total\nPeak Count'exist

if 'Total\nPeak Count' in data.columns:
    #Creat a nex file without any zero or NaN in colomne Total\nPeak Count
    data_zero = data[(data['Total\nPeak Count'] !=0) & (pd.notna(data['Total\nPeak Count']))]
    
else :
    print ('There is no "Total Peak Count" colomns in the file')

# Delete several columns
columns_to_delete = ['uniqueId','NmrAtom', 'ResidueType', 'Total\nPeak Count']
data_columns= data_zero.drop(columns=[col for col in columns_to_delete if col in data_zero.columns])

data_columns['SequenceCode']=data_columns['SequenceCode'].replace(sequence_code_dict)

data_columns= data_columns[~data_columns['SequenceCode'].astype(str).str.contains('-1')]

# Creat a new column with crescent number, begining by 1
data_columns.insert(0, 'Index', range(1, len(data_columns) + 1))

# Re order the columns
colonnes_ordre = ['Index', 'Value\n(ppm)', 'Value Error\n(ppm)', 'AtomName', 'SequenceCode']
data_columns = data_columns[colonnes_ordre]

# Exchange all the None in 0
data_columns.replace("None", 0, inplace=True)

#Save the file
data_columns.to_csv('name.prot', sep='\t', index=False)