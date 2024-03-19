def style_menu(menu_bar):
    # set the style sheet menu bar theme
    menu_bar.setStyleSheet('''
        QMenuBar {
            background-color: #333; 
            color: white;
        }
        QMenuBar::item {
            background-color: #333;
            color: white;
            padding: 4px 10px;
        }
        QMenuBar::item:selected {
            background-color: #555;
        }
    ''')