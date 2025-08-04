#!python3
# parameter
sequence = "APEKKVLFWYDPMKPDTKFDKPGKSPFMDMDLVPKYADESG" # put the sequence
output_file = "/home/biosys/Bureau/python_cy/SilB41.seq" # Path output
start_number = 43  # number beginning


# Dico aa conversion 1 letter to 3 letters
one_to_three = {
    'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
    'E': 'Glu', 'Q': 'Gln', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
    'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
    'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val'
}

def convert_sequence(seq):
    converted_seq = [one_to_three[aa] for aa in seq]
    return converted_seq

def write_to_file(seq_list, output_file, start_number):
    with open(output_file, 'w') as file:
        for i, aa in enumerate(seq_list):
            file.write(f"{aa}\t{start_number + i}\n")

converted_sequence = convert_sequence(sequence)

write_to_file(converted_sequence, output_file, start_number)
