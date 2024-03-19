import os
from PySide6.QtWidgets import QFileDialog
from path_handler import PathHandler
import re
from PySide6.QtGui import QColor

class FileBrowser:
    def __init__(self, button, line_edit):
        self.button = button
        self.line_edit = line_edit
        self.button.clicked.connect(self.open_file_dialog)
        self.pdb_files = []

    def open_file_dialog(self):
        file_dialog = QFileDialog()
        home_dir = os.path.expanduser("~")  # get the path to the home directory
        file_name = file_dialog.getExistingDirectory(None, "Select Directory", home_dir)
        if file_name:
            self.check_how_many_pdb_files(file_name)

    def check_how_many_pdb_files(self, file_name):
        handler = PathHandler(file_name)
        self.pdb_files, number = handler.get_pdb_files()
        self.line_edit.setText(f"{file_name}")
        style_sheet = self.line_edit.styleSheet()
        if number == 0:
            color = "#ff0000"
            print("No pdb files found")
        else:
            color = "#00ff00"
        self.line_edit.setStyleSheet(f"""QLineEdit {{
    background-color:  rgb(30, 30, 30);
    border: 1px solid #ccc;
    border-radius: 10px;
    padding: 8px;
    font-size: 11px;
    color: {color}
}}""")
