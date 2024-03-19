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
        
        all_files = self.file_browser.pdb_files
        print(all_files)
        
        # add more tabs

        # handle functionality of output window 

        # handle functionality of raw data window

        

if __name__ == "__main__":
    app = QApplication(sys.argv)
    widget = Main()
    widget.show()
    sys.exit(app.exec())
