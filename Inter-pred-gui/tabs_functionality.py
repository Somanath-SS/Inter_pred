class TabHandler:
    def __init__(self, tab_widget):
        self.tab_widget = tab_widget
        self.tab_widget.setTabsClosable(True)
        self.tab_widget.tabCloseRequested.connect(self.close_tab)

    def close_tab(self, index):
        self.tab_widget.removeTab(index)

    
