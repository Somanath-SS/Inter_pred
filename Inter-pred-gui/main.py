# This Python file uses the following encoding: utf-8
import sys

from PySide6.QtWidgets import QApplication, QMainWindow, QWidget, QVBoxLayout, QLabel
from PySide6.QtCore import Qt, QPoint 
from draggable_window import DraggableWindow
from title_window_functions import title_window_style
from menu_functions import style_menu
from tabs_functionality import TabHandler
from ui_form import Ui_Main 
from PySide6.QtWidgets import QFileDialog
from file_browser import FileBrowser
from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtWidgets import QVBoxLayout
from PySide6.QtWebEngineCore import QWebEngineSettings
import os
import json

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Backend import pdb_calc as pdb

class Main(QMainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.ui = Ui_Main()
        self.ui.setupUi(self)

        # set the style sheet main window theme
        title_window_style(self)

        # handle functionality of menu bar 
        menu_bar = self.menuBar() 
        style_menu(menu_bar)

        # File browsare button 

        # handle functionality of tab widgest
        self.tab_handler = TabHandler(self.ui.tabWidget)

        self.file_browser = FileBrowser(
            self.ui.FileBrowser, 
            self.ui.PathToTheFileLine, 
            self.tab_handler)
        
        # add more tabs

        # handle functionality of output window 

        # handle functionality of raw data window

        # handle functionality of analysis button
        self.ui.analizeButton.clicked.connect(self.analyze)
    
    def analyze(self):
        for tab_key in self.file_browser.tabs:
            tab = self.file_browser.tabs[tab_key]
            path = tab['file']
            protein_atoms, ligand_atoms, amino_acids = pdb.parse_pdb(path)
            distances = pdb.calculate_distances(protein_atoms, ligand_atoms)
            sorted_distances = pdb.sort_distances(distances)
            hydrogen_bonds = pdb.find_hydrogen_bonds(sorted_distances)
            visualize = pdb.visualize_ligand_sticks_with_labels(ligand_atoms, hydrogen_bonds, protein_atoms, amino_acids)
            # Get the HTML content of the visualization
            html_content = visualize.write_html()
            

            web_view = QWebEngineView()
            # Set the HTML content of the widget to the py3Dmol view
            web_view.setHtml(html_content)
            settings = web_view.settings()
            settings.setAttribute(
                QWebEngineSettings.WebAttribute.LocalContentCanAccessRemoteUrls,
                True
                )
        
            # # Add the QWebEngineView to the current tab
            tab['tab'].setLayout(QVBoxLayout())
            tab['tab'].layout().addWidget(web_view)
            tab['tab'].layout().setAlignment(Qt.AlignCenter)
            
                
 
if __name__ == "__main__":
    app = QApplication(sys.argv)
    widget = Main()
    widget.show()
    sys.exit(app.exec())
