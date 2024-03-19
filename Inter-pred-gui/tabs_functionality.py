from PySide6.QtWidgets import QWidget, QVBoxLayout, QLabel

class TabHandler:
    def __init__(self, tab_widget):
        super().__init__()
        self.tab_widget = tab_widget
        self.tab_widget.setTabsClosable(True)
        self.tab_widget.tabCloseRequested.connect(self.close_tab)

    def add_new_tab(self, file_name: str):
        self.tab_widget.addTab(QWidget(), file_name)

    def close_tab(self, index):
        self.tab_widget.removeTab(index)

    

    
