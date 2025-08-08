import fasta_to_seq
import file_to_prot
import to_xeasy

#run seq first
converted_sequence = fasta_to_seq.convert_sequence()
fasta_to_seq.write_to_file(seq_list=converted_sequence)



'''

# run file to prot
file_to_prot.process_protein_data()


#run to xeasy
to_xeasy.process_files()

print ('finsh')
'''