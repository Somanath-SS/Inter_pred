from PySide6.QtWidgets import QWidget, QVBoxLayout, QLabel

class RawAnalysisHandler():
    def __init__(self, scrollview):
        self.scrollview = scrollview
        self.scrollview.setWidget(self.create_scrollable_content())

    def create_scrollable_content(self):
            widget = QWidget()
            layout = QVBoxLayout()
            widget.setLayout(layout)

            # Create a new QLabel widget and store it in the text_widget attribute
            self.text_widget = QLabel("Raw Analysis")

            layout.addWidget(self.text_widget)
            return widget

    def add_text(self, text):
            self.text_widget.setText(text)
        