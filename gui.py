import tkinter as tk
import fasta_to_seq  # Import your conversion script
from tkinter import filedialog
import file_to_prot
import to_xeasy

# Global variables to store the inputs
output_file = ""  # String to store the Working_space path      /home/biosys/Bureau/python_cy/python_clean/SilB41.seq
start_number = 0 # Integer to store the number of amino acids
sequence = ""  # String to store the FASTA sequence   APEKKVLFWYDPMKPDTKFDKPGKSPFMDMDLVPKYADESG
attribution_prot="" #path to the file outpute of cyana
sequence_code_dict=''#path to the file forma json/dict
file_13=''#path to file NOE 13C
file_15=''#path to file NOE 15N


# Function to execute the conversion and writing when the "Launch" button is clicked
def run_sequence_conversion():
    global output_file, start_number, sequence
    start_number = int(number_entry.get())
    sequence = fasta_entry.get()
    output_file = save_entry.get()
    output = save_entry.get() + '/prot.seq'
    converted_sequence = fasta_to_seq.convert_sequence(sequence)  # Pass the sequence directly
    fasta_to_seq.write_to_file(converted_sequence, output, start_number)  # Pass the output file and start number
    
def fille_to_prot():
    global sequence_code_dict,attribution_prot, output_file
    sequence_code_dict = save_entry.get() + '/proline.txt'
    attribution_prot = save_entry.get() + '/attrib.csv'
    file_to_prot.process_protein_data(attribution_prot,sequence_code_dict,output_file)  # Pass the sequence directly
    
def xeasy_file():
    global choix1, output_file, file_13,file_15
    if choix1 == True:
        version = 2
    else:
        version = 3
    file_13 = output_file + '/13C.csv'
    file_15 = output_file + '/15N.csv'
    to_xeasy.process_files(file_13=file_13, file_15=file_15, version=version, path_save=output_file)


def launch_script():
    run_sequence_conversion()  # First command
    fille_to_prot()            # Second command
    xeasy_file()
    
def select_save_location():
    # Ouvre une boîte de dialogue pour sélectionner un répertoire
    folder_selected = filedialog.askdirectory()  # Demande à l'utilisateur de sélectionner un répertoire
    if folder_selected:  # Si un répertoire a été sélectionné
        save_entry.delete(0, tk.END)  # Supprime le texte actuel dans l'Entry
        save_entry.insert(0, folder_selected)  # Insère le chemin du répertoire dans l'Entry

def file_entery(to_execute, text_):
    # Create a section for save location input
    reasearch_file_label = tk.Label(root, text=text_)  # Label for save location entry
    reasearch_file_label.pack(pady=5)  # Add label with padding
    save_file_entry = tk.Entry(root, width=40)  # Entry field for save path
    save_file_entry.pack(pady=5)  # Add entry field with padding
    save_file_button = tk.Button(root, text="Browse Location", command=to_execute)  # Button to browse directories
    save_file_button.pack(pady=5)  # Add button with padding 
        

# Create the main application window
root = tk.Tk()
root.title("Cycy python")  # Set the window title
root.geometry("400x500")  # Set the window size
choix1 = tk.BooleanVar()

# Create a label to display the status or settings
status_label = tk.Label(root)  # Initial text for the status label
status_label.pack(pady=10)  # Add status label with padding

# Create a section for amino acid count input
number_label = tk.Label(root, text="Amino Acid Start:")  # Label for amino acid count entry
number_label.pack(pady=5)  # Add label with padding
number_entry = tk.Entry(root, width=10) 
number_entry.pack(pady=5)  # Add entry field with padding

# Create a section for amino acid count input
fasta_label = tk.Label(root, text="Sequence :")  # Label for amino acid count entry
fasta_label.pack(pady=5)  # Add label with padding
fasta_entry = tk.Entry(root, width=10) 
fasta_entry.pack(pady=5)  # Add entry field with padding

# Create a section for save location input
working_space_label = tk.Label(root, text="Working_space:")  # Label for save location entry
working_space_label.pack(pady=5)  # Add label with padding
save_entry = tk.Entry(root, width=40)  # Entry field for save path
save_entry.pack(pady=5)  # Add entry field with padding
save_button = tk.Button(root, text="Browse Location", command=select_save_location)  # Button to browse directories
save_button.pack(pady=5)  # Add button with padding

checkbox1 = tk.Checkbutton(root, text="Version 2", variable=choix1)
checkbox2 = tk.Checkbutton(root, text="Version 3")

checkbox1.pack(pady=10)
checkbox2.pack(pady=10)

# Create the "Launch" button to run the sequence conversion
launch_button = tk.Button(root, text="Launch", command=launch_script)  # Button triggers sequence conversion
launch_button.pack(pady=10)  # Add button with padding

# Start the main event loop to run the application
root.mainloop()






