from PySide6.QtCore import QPoint, Qt, QEvent 
from PySide6.QtWidgets import QWidget, QMenuBar

class DraggableWindow(QMenuBar):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.moving = False
        self.offset = QPoint()

    def mousePressEvent(self, event):
        if event.buttons() == Qt.LeftButton:
            self.moving = True  # Changed from self.dragging to self.moving
            self.offset = event.globalPos() - self.parent().pos()

    def mouseMoveEvent(self, event):
        if self.moving and event.buttons() == Qt.LeftButton:
            self.parent().move(event.globalPos() - self.offset)

    def mouseReleaseEvent(self, event):
        if event.button() == Qt.LeftButton:
            self.moving = False  # Changed from self.dragging to self.moving