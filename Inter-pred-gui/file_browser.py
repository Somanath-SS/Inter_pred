import os
from PySide6.QtWidgets import QFileDialog
from path_handler import PathHandler
import re
from PySide6.QtGui import QColor
from tabs_functionality import TabHandler
from PySide6.QtWidgets import QWidget

class FileBrowser:
    """
    This class is used to handle the file browser functionality
    
    Args:
        button (QWidget): The button to open the file dialog
        line_edit (QWidget): The line edit to display the path to the file
        tab_handler (TabHandler): The tab handler to handle the tabs
        
    Returns:
        None: This class is used to handle the file browser functionality thus it returns nothing
    
    Functions:
        open_file_dialog: This function is used to open the file dialog
        check_how_many_pdb_files: This function is used to check how many pdb files are in the directory
        add_more_tabs: This function is used to add more tabs
    """
    def __init__(self, button, line_edit, tab_handler: TabHandler):
        self.button = button
        self.line_edit = line_edit
        self.button.clicked.connect(self.open_file_dialog)
        self.pdb_files = []
        self.tab_handler = tab_handler

    def open_file_dialog(self):
        """
        This function is used to open the file dialog

        Args:
            None: This function does not take any arguments

        Returns:
            None: This function does not return anything

        Raises:
            None: This function does not raise any exceptions

        """
        file_dialog = QFileDialog()
        home_dir = os.path.expanduser("~")  # get the path to the home directory
        file_name = file_dialog.getExistingDirectory(None, "Select Directory", home_dir)
        if file_name:
            self.check_how_many_pdb_files(file_name)

    def check_how_many_pdb_files(self, file_name):
        """
        This function is used to check how many pdb files are in the directory
        
        Args:
            file_name (str): The name of the file
            
        Returns:
            None: This function does not return anything

        Raises:
            None: This function does not raise any exceptions

            """
        handler = PathHandler(file_name)
        self.pdb_files, number = handler.get_pdb_files()
        self.line_edit.setText(f"{file_name}")
        if number == 0:
            color = "#ff0000"
            print("No pdb files found")
        else:
            color = "#00ff00"
        self.line_edit.setStyleSheet(
            f"""QLineEdit {{
            background-color:  rgb(30, 30, 30);
            border: 1px solid #ccc;
            border-radius: 10px;
            padding: 8px;
            font-size: 11px;
            color: {color}
            }}""")
        self.add_more_tabs(self.pdb_files, self.tab_handler)
        
    def add_more_tabs(self, all_files: PathHandler, tab_handler: TabHandler):
        """
        This function is used to add more tabs
        
        Args:
            all_files (PathHandler): The path handler object
            tab_handler (TabHandler): The tab handler object
        
        Returns:
            None: This function does not return anything
            
        Raises:
            None: This function does not raise any exceptions
        """
        for file in all_files:
            tab_handler.add_new_tab(file.name)



        
