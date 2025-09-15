import sys
import os
import subprocess
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, 
                             QLabel, QLineEdit, QPushButton, QRadioButton, QFileDialog, 
                             QProgressBar, QMessageBox)
from PyQt6.QtCore import Qt, QThread, pyqtSignal

class ScriptThread(QThread):
    """Thread to run the external Python script."""
    finished = pyqtSignal(str)  # Signal to send completion message
    error = pyqtSignal(str)     # Signal to send error message

    def __init__(self, script_path, args):
        super().__init__()
        self.script_path = script_path
        self.args = args

    def run(self):
        try:
            # Run the external script with arguments
            result = subprocess.run(
                [sys.executable, self.script_path] + self.args,
                capture_output=True,
                text=True,
                check=True
            )
            self.finished.emit("Processing complete! Output: " + result.stdout)
        except subprocess.CalledProcessError as e:
            self.error.emit(f"Script error: {e.stderr}")
        except Exception as e:
            self.error.emit(f"Error: {str(e)}")

class CycyApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("CCPNMR to Cyana")
        self.setFixedSize(450, 550)

        # Main widget and layout
        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)
        self.main_layout = QVBoxLayout(self.central_widget)
        self.main_layout.setContentsMargins(20, 20, 20, 20)
        self.main_layout.setSpacing(15)

        # Apply dark theme with purple accents
        self.setStyleSheet("""
            QMainWindow { background-color: #2E2E2E; }
            QLabel { color: #E0E0E0; font: 11pt "Segoe UI"; }
            QLineEdit { 
                background-color: #3C3C3C; 
                color: #E0E0E0; 
                border: 1px solid #7B3FE4; 
                border-radius: 5px; 
                padding: 5px; 
            }
            QPushButton { 
                background-color: #7B3FE4; 
                color: #E0E0E0; 
                font: bold 10pt "Segoe UI"; 
                border-radius: 5px; 
                padding: 8px; 
            }
            QPushButton:hover { background-color: #A78BFA; }
            QPushButton:pressed { background-color: #5F2EE6; }
            QRadioButton { color: #E0E0E0; font: 10pt "Segoe UI"; }
            QRadioButton::indicator:checked { background-color: #7B3FE4; }
            QProgressBar { 
                background-color: #3C3C3C; 
                color: #E0E0E0; 
                border-radius: 5px; 
            }
            QProgressBar::chunk { background-color: #7B3FE4; }
        """)

        # Create UI elements
        self.create_widgets()

    def create_widgets(self):
        # Title label
        title_label = QLabel("CCPNMRV3 to Cyana", alignment=Qt.AlignmentFlag.AlignCenter)
        title_label.setStyleSheet("font: bold 14pt 'Segoe UI'; color: #A78BFA;")
        self.main_layout.addWidget(title_label)

        # Status label and progress bar
        self.status_label = QLabel("Ready", alignment=Qt.AlignmentFlag.AlignCenter)
        self.main_layout.addWidget(self.status_label)
        self.progress = QProgressBar()
        self.progress.setTextVisible(False)
        self.progress.setFixedHeight(10)
        self.main_layout.addWidget(self.progress)
        self.progress.hide()

        # Input section
        self.number_input = QLineEdit()
        self.number_input.setPlaceholderText("1")
        self.main_layout.addWidget(QLabel("Amino Acid Start:"))
        self.main_layout.addWidget(self.number_input)

        self.fasta_input = QLineEdit()
        self.fasta_input.setPlaceholderText("APEKKVLFWYDPMKPDTKFDKPGKSPFMDMDLVPKYADESG")
        self.main_layout.addWidget(QLabel("Enter protein sequence:"))
        self.main_layout.addWidget(self.fasta_input)

        self.save_input = QLineEdit()
        self.save_input.setPlaceholderText("Select working directory")
        self.main_layout.addWidget(QLabel("Working Space:"))
        self.main_layout.addWidget(self.save_input)
        browse_btn = QPushButton("Browse Location")
        browse_btn.clicked.connect(self.select_save_location)
        self.main_layout.addWidget(browse_btn)

        # Version selection
        version_layout = QHBoxLayout()
        self.version_2 = QRadioButton("Version 2")
        self.version_2.setChecked(True)
        self.version_3 = QRadioButton("Version 3")
        version_layout.addWidget(QLabel("XEASY Version:"))
        version_layout.addWidget(self.version_2)
        version_layout.addWidget(self.version_3)
        self.main_layout.addLayout(version_layout)

        # Launch button
        launch_btn = QPushButton("Launch")
        launch_btn.clicked.connect(self.launch_script)
        self.main_layout.addWidget(launch_btn)

        # Add stretch to center content
        self.main_layout.addStretch()

    def select_save_location(self):
        folder = QFileDialog.getExistingDirectory(self, "Select Directory")
        if folder:
            self.save_input.setText(folder)

    def launch_script(self):
        # Validate inputs
        start_number = self.number_input.text().strip()
        fasta_sequence = self.fasta_input.text().strip()
        save_path = self.save_input.text().strip()
        version = "2" if self.version_2.isChecked() else "3"
        file_13C = os.path.join(save_path, "13C.csv")
        file_15N = os.path.join(save_path, "15N.csv")
        cyana_attrib_input = os.path.join(save_path, "attrib.csv")

        # Input validation
        if not all([start_number, fasta_sequence, save_path]):
            QMessageBox.critical(self, "Error", "All fields are required!")
            return

        try:
            start_num = int(start_number)
            if start_num < 1:
                raise ValueError("Start number must be a positive integer")
        except ValueError:
            QMessageBox.critical(self, "Error", "Invalid start number! Must be a positive integer.")
            return

        valid_aas = set("ARNDCQEGHILKMFPSTWYV")
        if not all(aa in valid_aas for aa in fasta_sequence.upper()):
            QMessageBox.critical(self, "Error", "Invalid protein sequence! Use standard one-letter amino acid codes.")
            return

        if not os.path.isdir(save_path):
            QMessageBox.critical(self, "Error", "Invalid save directory!")
            return

        # Define output file path
        output_file = os.path.join(save_path, "prot.seq")
        output_csv = os.path.join(save_path, "attib_cyana.prot")

        # Path to the external Python script
        fasta_to_seq = "fasta_to_seq.py"  # Adjust if the script is in a different directory
        to_xeasy = "to_xeasy.py"
        file_toprot = "file_to_prot.py"
        # Prepare arguments for the script
        args_fasta_to_seq = [start_number, fasta_sequence, output_file]
        args_to_xeasy = [file_13C, file_15N, version, save_path]

        # Start the script in a separate thread
        self.status_label.setText("Processing...")
        self.progress.show()
        self.progress.setRange(0, 0)  # Indeterminate mode
        self.script_thread_fasta_to_seq = ScriptThread(fasta_to_seq, args_fasta_to_seq)
        self.script_thread_fasta_to_seq.finished.connect(self.on_script_finished)
        self.script_thread_fasta_to_seq.error.connect(self.on_script_error)
        self.script_thread_fasta_to_seq.start()

        seq_file = os.path.join(save_path, "prot.seq")
        args_file_to_prot= [cyana_attrib_input, seq_file,save_path, output_csv]


        self.script_thread_xeasy=ScriptThread(to_xeasy,args_to_xeasy)
        self.script_thread_xeasy.finished.connect(self.on_script_finished)
        self.script_thread_xeasy.error.connect(self.on_script_error)
        self.script_thread_xeasy.start()

        self.script_thread_to_prot=ScriptThread(file_toprot,args_file_to_prot)
        self.script_thread_to_prot.finished.connect(self.on_script_finished)
        self.script_thread_to_prot.error.connect(self.on_script_error)
        self.script_thread_to_prot.start()

    def on_script_finished(self, message):
        self.status_label.setText("Processing complete!")
        self.progress.hide()
        #QMessageBox.information(self, "Success", message)

    def on_script_error(self, message):
        self.status_label.setText("Error occurred")
        self.progress.hide()
        QMessageBox.critical(self, "Error", message)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = CycyApp()
    window.show()
    sys.exit(app.exec())