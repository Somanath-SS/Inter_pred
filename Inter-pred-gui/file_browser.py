import os
from PySide6.QtWidgets import QFileDialog

class FileBrowser:
    def __init__(self, button):
        self.button = button
        self.button.clicked.connect(self.open_file_dialog)

    def open_file_dialog(self):
        file_dialog = QFileDialog()
        home_dir = os.path.expanduser("~")  # get the path to the home directory
        file_name = file_dialog.getExistingDirectory(None, "Select Directory", home_dir)
        if file_name:
            print(f"Selected directory: {file_name}")