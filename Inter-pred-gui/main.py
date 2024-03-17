# This Python file uses the following encoding: utf-8
import sys

from PySide6.QtWidgets import QApplication, QMainWindow, QWidget, QVBoxLayout, QLabel
from PySide6.QtCore import Qt, QPoint 
from draggable_window import DraggableWindow

from ui_form import Ui_Main

class Main(QMainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.ui = Ui_Main()
        self.ui.setupUi(self)
                
        self.setWindowFlags(Qt.FramelessWindowHint)

        # menu bar 
        menu_bar = self.menuBar() 
        
        self.draggable_window = DraggableWindow(menu_bar)
        # set the style sheet menu bar theme
        menu_bar.setStyleSheet('''
            QMenuBar {
                background-color: #333; /* Dark Gray */
                color: white;
            }
            QMenuBar::item {
                background-color: #333; /* Dark Gray */
                color: white;
                padding: 4px 10px;
            }
            QMenuBar::item:selected {
                background-color: #555; /* Dark Gray (Hovered) */
            }
        ''')


if __name__ == "__main__":
    app = QApplication(sys.argv)
    widget = Main()
    widget.show()
    sys.exit(app.exec())
