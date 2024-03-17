from PySide6.QtCore import QPoint, Qt 
from PySide6.QtWidgets import QWidget

class DraggableWindow(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.dragging = False
        self.offset = QPoint()

    def mousePressEvent(self, event):
        if event.buttons() == Qt.LeftButton:
            self.dragging = True
            self.offset = event.globalPos() - self.parent().pos()

    def mouseMoveEvent(self, event):
        if self.dragging:
            self.parent().move(event.globalPos() - self.offset)

    def mouseReleaseEvent(self, event):
        if event.button() == Qt.LeftButton:
            self.dragging = False