#! /usr/bin/python3
import sys 
import pandas as pd
'''
file_15 = '/home/biosys/Bureau/python_cy/for_cyana_ccpnm/Bform/15NNOE.csv'
file_13 = '/home/biosys/Bureau/python_cy/for_cyana_ccpnm/Bform/C13NOE.csv'
path_save = '/home/biosys/Bureau/python_cy/python_clean/'
version= '3' #input('cyana 2 vs cyana 3.98 (2/3)')
'''
def process_files(file_13, file_15, version, path_save):
    def change_format(file, path_save=path_save, variable='N'):
        file = file.drop(columns=['_object']) # remove the first colomn _object
        file.insert(0, 'Index', range(1, len(file) + 1)) # add number increasing
        file.insert(4,'1',1)
        file.insert(5,'U','U')
        file = file.replace('None', 0)
        file['Column_0.00e+00'] = '0.00e+00'
        for i in range(1, 6):
            file[f'Column_0_{i}'] = 0
        file = file.to_csv(index=False, header=False, sep='\t')
        
        if version == 3:
            text_n = "# Number of dimensions 3\n#FORMAT xeasy3D\n#INAME 1 HN\n#INAME 2 H\n#INAME 3 N\n#SPECTRUM N15NOESY  HN H N\n"
            text_c = "# Number of dimensions 3\n#FORMAT xeasy3D\n#INAME 1 HC\n#INAME 2 H\n#INAME 3 C\n#SPECTRUM C13NOESY  HC H C\n"
        elif version == 2:
            text_n = "# Number of dimensions 3\n#FORMAT xeasy3D\n#INAME 1 HN\n#INAME 2 H\n#INAME 3 N\n#CYANAFORMAT HhN\n"
            text_c = "# Number of dimensions 3\n#FORMAT xeasy3D\n#INAME 1 HC\n#INAME 2 H\n#INAME 3 C\n#CYANAFORMAT HhC\n"
            
        if variable == 'N':
                with open(path_save +'/15N.peaks', 'w') as text:
                    text.write(text_n)
                    text.write(file)
            
        else:
            with open(path_save +'/13C.peaks', 'w') as text:
                text.write(text_c)
                text.write(file)
                    
    
    if version == 3 or 2:
        data_13 = pd.read_csv(file_13)
        data_15 = pd.read_csv(file_15)
        
        #change columns for 13C
        new_order = [ '_object', 'Pos F1', 'Pos F3', 'Pos F2', 'Volume'] # new order for data in carbon
        data_13 = data_13[new_order] #change the order
        
        #formatation
        change_format(file=data_13, variable='C')
        change_format(file=data_15)
    
    else:
        print('out')

def main():
    file_13C = sys.argv[1]
    File_15N = sys.argv[2]
    version = sys.argv[3]
    save_path = sys.argv[4]

    
    version = int(version)
    process_files(file_13C,File_15N,version,save_path)

if __name__ == "__main__":
    main()