from PySide6.QtWidgets import QWidget, QVBoxLayout, QLabel

class TabHandler:
    """
    This class is used to handle the tab functionality
    
    Args:
        tab_widget (QWidget): The tab widget
        
    Functions: 
        add_new_tab: This function is used to add a new tab
        close_tab: This function is used to close a tab
    """
    def __init__(self, tab_widget):
        super().__init__()
        self.tab_widget = tab_widget
        self.tab_widget.setTabsClosable(True)
        self.tab_widget.tabCloseRequested.connect(self.close_tab)

    def add_new_tab(self, file_name: str):
        """
        This function is used to add a new tab
        
        Args:
            file_name (str): The name of the file
            
        Returns:
            None: This function does not return anything
        """
        self.tab_widget.addTab(QWidget(), file_name)

    def close_tab(self, index):
        """
        This function is used to close a tab
        
        Args:
            index (int): The index of the tab
        
        Returns:
            None: This function does not return anything
        """
        self.tab_widget.removeTab(index)

    

    
