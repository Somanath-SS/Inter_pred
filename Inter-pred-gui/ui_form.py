# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'form.ui'
##
## Created by: Qt User Interface Compiler version 6.6.2
##
## WARNING! All changes made in this file will be lost when recompiling UI file!
################################################################################

from PySide6.QtCore import (QCoreApplication, QDate, QDateTime, QLocale,
    QMetaObject, QObject, QPoint, QRect,
    QSize, QTime, QUrl, Qt)
from PySide6.QtGui import (QAction, QBrush, QColor, QConicalGradient,
    QCursor, QFont, QFontDatabase, QGradient,
    QIcon, QImage, QKeySequence, QLinearGradient,
    QPainter, QPalette, QPixmap, QRadialGradient,
    QTransform)
from PySide6.QtWidgets import (QApplication, QHBoxLayout, QLineEdit, QMainWindow,
    QMenu, QMenuBar, QPushButton, QScrollArea,
    QSizePolicy, QSplitter, QStatusBar, QTabWidget,
    QVBoxLayout, QWidget)

class Ui_Main(object):
    def setupUi(self, Main):
        if not Main.objectName():
            Main.setObjectName(u"Main")
        Main.resize(931, 510)
        palette = QPalette()
        brush = QBrush(QColor(244, 244, 244, 255))
        brush.setStyle(Qt.SolidPattern)
        palette.setBrush(QPalette.Active, QPalette.WindowText, brush)
        brush1 = QBrush(QColor(22, 95, 219, 255))
        brush1.setStyle(Qt.SolidPattern)
        palette.setBrush(QPalette.Active, QPalette.Button, brush1)
        brush2 = QBrush(QColor(45, 45, 45, 255))
        brush2.setStyle(Qt.SolidPattern)
        palette.setBrush(QPalette.Active, QPalette.Light, brush2)
        brush3 = QBrush(QColor(37, 37, 37, 255))
        brush3.setStyle(Qt.SolidPattern)
        palette.setBrush(QPalette.Active, QPalette.Midlight, brush3)
        brush4 = QBrush(QColor(15, 15, 15, 255))
        brush4.setStyle(Qt.SolidPattern)
        palette.setBrush(QPalette.Active, QPalette.Dark, brush4)
        brush5 = QBrush(QColor(20, 20, 20, 255))
        brush5.setStyle(Qt.SolidPattern)
        palette.setBrush(QPalette.Active, QPalette.Mid, brush5)
        palette.setBrush(QPalette.Active, QPalette.Text, brush)
        palette.setBrush(QPalette.Active, QPalette.BrightText, brush)
        palette.setBrush(QPalette.Active, QPalette.ButtonText, brush)
        brush6 = QBrush(QColor(30, 30, 30, 255))
        brush6.setStyle(Qt.SolidPattern)
        palette.setBrush(QPalette.Active, QPalette.Base, brush6)
        brush7 = QBrush(QColor(40, 40, 40, 255))
        brush7.setStyle(Qt.SolidPattern)
        palette.setBrush(QPalette.Active, QPalette.Window, brush7)
        brush8 = QBrush(QColor(0, 0, 0, 255))
        brush8.setStyle(Qt.SolidPattern)
        palette.setBrush(QPalette.Active, QPalette.Shadow, brush8)
        brush9 = QBrush(QColor(18, 80, 180, 255))
        brush9.setStyle(Qt.SolidPattern)
        palette.setBrush(QPalette.Active, QPalette.Highlight, brush9)
        brush10 = QBrush(QColor(198, 198, 198, 255))
        brush10.setStyle(Qt.SolidPattern)
        palette.setBrush(QPalette.Inactive, QPalette.WindowText, brush10)
        palette.setBrush(QPalette.Inactive, QPalette.Button, brush9)
        palette.setBrush(QPalette.Inactive, QPalette.Light, brush2)
        palette.setBrush(QPalette.Inactive, QPalette.Midlight, brush3)
        palette.setBrush(QPalette.Inactive, QPalette.Dark, brush4)
        palette.setBrush(QPalette.Inactive, QPalette.Mid, brush5)
        palette.setBrush(QPalette.Inactive, QPalette.Text, brush10)
        palette.setBrush(QPalette.Inactive, QPalette.BrightText, brush10)
        palette.setBrush(QPalette.Inactive, QPalette.ButtonText, brush10)
        brush11 = QBrush(QColor(21, 21, 21, 255))
        brush11.setStyle(Qt.SolidPattern)
        palette.setBrush(QPalette.Inactive, QPalette.Base, brush11)
        palette.setBrush(QPalette.Inactive, QPalette.Window, brush6)
        palette.setBrush(QPalette.Inactive, QPalette.Shadow, brush8)
        palette.setBrush(QPalette.Inactive, QPalette.Highlight, brush9)
        brush12 = QBrush(QColor(172, 172, 172, 255))
        brush12.setStyle(Qt.SolidPattern)
        palette.setBrush(QPalette.Disabled, QPalette.WindowText, brush12)
        brush13 = QBrush(QColor(12, 53, 118, 255))
        brush13.setStyle(Qt.SolidPattern)
        palette.setBrush(QPalette.Disabled, QPalette.Button, brush13)
        palette.setBrush(QPalette.Disabled, QPalette.Light, brush2)
        palette.setBrush(QPalette.Disabled, QPalette.Midlight, brush3)
        palette.setBrush(QPalette.Disabled, QPalette.Dark, brush4)
        palette.setBrush(QPalette.Disabled, QPalette.Mid, brush5)
        palette.setBrush(QPalette.Disabled, QPalette.Text, brush12)
        palette.setBrush(QPalette.Disabled, QPalette.BrightText, brush12)
        palette.setBrush(QPalette.Disabled, QPalette.ButtonText, brush12)
        palette.setBrush(QPalette.Disabled, QPalette.Base, brush4)
        palette.setBrush(QPalette.Disabled, QPalette.Window, brush6)
        palette.setBrush(QPalette.Disabled, QPalette.Shadow, brush8)
        palette.setBrush(QPalette.Disabled, QPalette.Highlight, brush13)
        Main.setPalette(palette)
        Main.setTabShape(QTabWidget.Rounded)
        self.actionSave = QAction(Main)
        self.actionSave.setObjectName(u"actionSave")
        self.actionHow_to_use = QAction(Main)
        self.actionHow_to_use.setObjectName(u"actionHow_to_use")
        self.actionCurrent_version = QAction(Main)
        self.actionCurrent_version.setObjectName(u"actionCurrent_version")
        self.actionGive_us_feedback = QAction(Main)
        self.actionGive_us_feedback.setObjectName(u"actionGive_us_feedback")
        self.actionMore_analysis_options = QAction(Main)
        self.actionMore_analysis_options.setObjectName(u"actionMore_analysis_options")
        self.actionCurrent_version_2 = QAction(Main)
        self.actionCurrent_version_2.setObjectName(u"actionCurrent_version_2")
        self.actionAbout_us = QAction(Main)
        self.actionAbout_us.setObjectName(u"actionAbout_us")
        self.actionExcel_file = QAction(Main)
        self.actionExcel_file.setObjectName(u"actionExcel_file")
        self.actionRaw_data = QAction(Main)
        self.actionRaw_data.setObjectName(u"actionRaw_data")
        self.actioncsv = QAction(Main)
        self.actioncsv.setObjectName(u"actioncsv")
        self.actionjson = QAction(Main)
        self.actionjson.setObjectName(u"actionjson")
        self.actionTheme = QAction(Main)
        self.actionTheme.setObjectName(u"actionTheme")
        self.centralwidget = QWidget(Main)
        self.centralwidget.setObjectName(u"centralwidget")
        self.horizontalLayout_3 = QHBoxLayout(self.centralwidget)
        self.horizontalLayout_3.setObjectName(u"horizontalLayout_3")
        self.splitter = QSplitter(self.centralwidget)
        self.splitter.setObjectName(u"splitter")
        self.splitter.setOrientation(Qt.Horizontal)
        self.layoutWidget = QWidget(self.splitter)
        self.layoutWidget.setObjectName(u"layoutWidget")
        self.horizontalLayout_2 = QHBoxLayout(self.layoutWidget)
        self.horizontalLayout_2.setObjectName(u"horizontalLayout_2")
        self.horizontalLayout_2.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_2 = QVBoxLayout()
        self.verticalLayout_2.setObjectName(u"verticalLayout_2")
        self.verticalLayout = QVBoxLayout()
        self.verticalLayout.setObjectName(u"verticalLayout")
        self.tabWidget = QTabWidget(self.layoutWidget)
        self.tabWidget.setObjectName(u"tabWidget")
        self.tabWidget.setEnabled(True)
        font = QFont()
        font.setFamilies([u"Cascadia Mono ExtraLight"])
        self.tabWidget.setFont(font)
        self.tabWidget.setStyleSheet(u"   QTabWidget::tab-bar {\n"
"           alignment: left;\n"
"        }\n"
"        QTabBar::tab {\n"
"            background-color: #rgb(240, 240, 240); \n"
"            color: rgb(240, 240, 240);              \n"
"            padding: 8px 20px;           \n"
"            border-top-left-radius: 4px; \n"
"            border-top-right-radius: 4px;\n"
"            border: 1px solid white;\n"
"        }\n"
"        QTabBar::tab:selected {\n"
"            color: rgb(0, 255, 0);              \n"
"        }")
        self.tabWidget.setElideMode(Qt.ElideLeft)
        self.tabWidget.setDocumentMode(False)
        self.tabWidget.setTabsClosable(False)
        self.tabWidget.setMovable(True)
        self.tabWidget.setTabBarAutoHide(False)
        self.Starting = QWidget()
        self.Starting.setObjectName(u"Starting")
        self.tabWidget.addTab(self.Starting, "")

        self.verticalLayout.addWidget(self.tabWidget)


        self.verticalLayout_2.addLayout(self.verticalLayout)

        self.horizontalLayout = QHBoxLayout()
        self.horizontalLayout.setObjectName(u"horizontalLayout")
        self.PathToTheFileLine = QLineEdit(self.layoutWidget)
        self.PathToTheFileLine.setObjectName(u"PathToTheFileLine")
        font1 = QFont()
        font1.setFamilies([u"Cascadia Code"])
        self.PathToTheFileLine.setFont(font1)
        self.PathToTheFileLine.setAutoFillBackground(False)
        self.PathToTheFileLine.setStyleSheet(u"        QLineEdit {\n"
"            background-color:  rgb(30, 30, 30); /* Light Gray */\n"
"            border: 1px solid #ccc; /* Light Gray */\n"
"            border-radius: 10px;\n"
"            padding: 8px;\n"
"            font-size: 11px;\n"
"            color: rgb(0, 255, 0) /* Dark Gray */\n"
"        }")
        self.PathToTheFileLine.setReadOnly(True)

        self.horizontalLayout.addWidget(self.PathToTheFileLine)

        self.FileBrowser = QPushButton(self.layoutWidget)
        self.FileBrowser.setObjectName(u"FileBrowser")
        font2 = QFont()
        font2.setFamilies([u"Cascadia Code SemiBold"])
        font2.setBold(True)
        self.FileBrowser.setFont(font2)
        self.FileBrowser.setCursor(QCursor(Qt.ArrowCursor))
        self.FileBrowser.setStyleSheet(u"        QPushButton {\n"
"            background-color: #rgb(30, 30, 30);\n"
"            color: white;\n"
"            font-size: 14px;\n"
"            margin: 4px 2px;\n"
"            border-radius: 10px;\n"
"            border: 1px solid white;\n"
"            width: 130px;\n"
"            height: 30px\n"
"\n"
"        }\n"
"        QPushButton:hover {\n"
"            border: 2px #282828;\n"
"        }")

        self.horizontalLayout.addWidget(self.FileBrowser)

        self.analizeButton = QPushButton(self.layoutWidget)
        self.analizeButton.setObjectName(u"analizeButton")
        self.analizeButton.setFont(font2)
        self.analizeButton.setStyleSheet(u"        QPushButton {\n"
"            background-color: #rgb(30, 30, 30);\n"
"            color: white;\n"
"            font-size: 14px;\n"
"            margin: 4px 2px;\n"
"            border-radius: 10px;\n"
"            border: 1px solid white;\n"
"            width: 100px;\n"
"            height: 30px\n"
"\n"
"        }\n"
"        QPushButton:hover {\n"
"            border: 2px #282828;\n"
"        }")

        self.horizontalLayout.addWidget(self.analizeButton)

        self.errorButton = QPushButton(self.layoutWidget)
        self.errorButton.setObjectName(u"errorButton")
        font3 = QFont()
        font3.setFamilies([u"Cascadia Mono SemiBold"])
        font3.setBold(True)
        self.errorButton.setFont(font3)
        self.errorButton.setAutoFillBackground(False)
        self.errorButton.setStyleSheet(u"        QPushButton {\n"
"            background-color: #rgb(30, 30, 30);\n"
"            color: white;\n"
"            font-size: 14px;\n"
"            margin: 4px 2px;\n"
"            border-radius: 10px;\n"
"            border: 1px solid white;\n"
"            width: 100px;\n"
"            height: 30px\n"
"\n"
"        }\n"
"        QPushButton:hover {\n"
"            border: 2px #282828;\n"
"        }")
        self.errorButton.setCheckable(False)

        self.horizontalLayout.addWidget(self.errorButton)


        self.verticalLayout_2.addLayout(self.horizontalLayout)


        self.horizontalLayout_2.addLayout(self.verticalLayout_2)

        self.splitter.addWidget(self.layoutWidget)
        self.layoutWidget1 = QWidget(self.splitter)
        self.layoutWidget1.setObjectName(u"layoutWidget1")
        self.verticalLayout_3 = QVBoxLayout(self.layoutWidget1)
        self.verticalLayout_3.setObjectName(u"verticalLayout_3")
        self.verticalLayout_3.setContentsMargins(0, 0, 0, 0)
        self.CodeView = QScrollArea(self.layoutWidget1)
        self.CodeView.setObjectName(u"CodeView")
        self.CodeView.setFont(font1)
        self.CodeView.setWidgetResizable(True)
        self.scrollAreaWidgetContents = QWidget()
        self.scrollAreaWidgetContents.setObjectName(u"scrollAreaWidgetContents")
        self.scrollAreaWidgetContents.setGeometry(QRect(0, 0, 112, 219))
        self.CodeView.setWidget(self.scrollAreaWidgetContents)

        self.verticalLayout_3.addWidget(self.CodeView)

        self.RawAnalysis = QScrollArea(self.layoutWidget1)
        self.RawAnalysis.setObjectName(u"RawAnalysis")
        self.RawAnalysis.setFont(font1)
        self.RawAnalysis.setAutoFillBackground(False)
        self.RawAnalysis.setWidgetResizable(True)
        self.scrollAreaWidgetContents_2 = QWidget()
        self.scrollAreaWidgetContents_2.setObjectName(u"scrollAreaWidgetContents_2")
        self.scrollAreaWidgetContents_2.setGeometry(QRect(0, 0, 112, 218))
        self.RawAnalysis.setWidget(self.scrollAreaWidgetContents_2)

        self.verticalLayout_3.addWidget(self.RawAnalysis)

        self.splitter.addWidget(self.layoutWidget1)

        self.horizontalLayout_3.addWidget(self.splitter)

        Main.setCentralWidget(self.centralwidget)
        self.menubar = QMenuBar(Main)
        self.menubar.setObjectName(u"menubar")
        self.menubar.setGeometry(QRect(0, 0, 931, 22))
        self.menuFile = QMenu(self.menubar)
        self.menuFile.setObjectName(u"menuFile")
        self.menuSave_as = QMenu(self.menuFile)
        self.menuSave_as.setObjectName(u"menuSave_as")
        self.menuOptions = QMenu(self.menubar)
        self.menuOptions.setObjectName(u"menuOptions")
        self.menuHelp = QMenu(self.menubar)
        self.menuHelp.setObjectName(u"menuHelp")
        self.menuInformation = QMenu(self.menubar)
        self.menuInformation.setObjectName(u"menuInformation")
        self.menuGive_us_feedback = QMenu(self.menubar)
        self.menuGive_us_feedback.setObjectName(u"menuGive_us_feedback")
        Main.setMenuBar(self.menubar)
        self.statusbar = QStatusBar(Main)
        self.statusbar.setObjectName(u"statusbar")
        Main.setStatusBar(self.statusbar)

        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuOptions.menuAction())
        self.menubar.addAction(self.menuHelp.menuAction())
        self.menubar.addAction(self.menuInformation.menuAction())
        self.menubar.addAction(self.menuGive_us_feedback.menuAction())
        self.menuFile.addAction(self.actionSave)
        self.menuFile.addAction(self.menuSave_as.menuAction())
        self.menuFile.addSeparator()
        self.menuFile.addAction(self.actionMore_analysis_options)
        self.menuSave_as.addAction(self.actionExcel_file)
        self.menuSave_as.addAction(self.actionRaw_data)
        self.menuSave_as.addAction(self.actioncsv)
        self.menuSave_as.addAction(self.actionjson)
        self.menuOptions.addAction(self.actionTheme)
        self.menuHelp.addAction(self.actionHow_to_use)
        self.menuInformation.addAction(self.actionCurrent_version_2)
        self.menuInformation.addAction(self.actionAbout_us)

        self.retranslateUi(Main)

        self.tabWidget.setCurrentIndex(1)


        QMetaObject.connectSlotsByName(Main)
    # setupUi

    def retranslateUi(self, Main):
        Main.setWindowTitle(QCoreApplication.translate("Main", u"Main", None))
        self.actionSave.setText(QCoreApplication.translate("Main", u"Save", None))
        self.actionHow_to_use.setText(QCoreApplication.translate("Main", u"How to use", None))
        self.actionCurrent_version.setText(QCoreApplication.translate("Main", u"Current version", None))
        self.actionGive_us_feedback.setText(QCoreApplication.translate("Main", u"Give us feedback!", None))
        self.actionMore_analysis_options.setText(QCoreApplication.translate("Main", u"More analysis options", None))
        self.actionCurrent_version_2.setText(QCoreApplication.translate("Main", u"Current version", None))
        self.actionAbout_us.setText(QCoreApplication.translate("Main", u"About us", None))
        self.actionExcel_file.setText(QCoreApplication.translate("Main", u"Excel file", None))
        self.actionRaw_data.setText(QCoreApplication.translate("Main", u"Raw data", None))
        self.actioncsv.setText(QCoreApplication.translate("Main", u"csv", None))
        self.actionjson.setText(QCoreApplication.translate("Main", u"json", None))
        self.actionTheme.setText(QCoreApplication.translate("Main", u"Theme", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.Starting), QCoreApplication.translate("Main", u"Starting", None))
        self.PathToTheFileLine.setPlaceholderText(QCoreApplication.translate("Main", u"Path to the file", None))
        self.FileBrowser.setText(QCoreApplication.translate("Main", u"File browser", None))
        self.analizeButton.setText(QCoreApplication.translate("Main", u"Analize", None))
        self.errorButton.setText(QCoreApplication.translate("Main", u"Errors", None))
        self.menuFile.setTitle(QCoreApplication.translate("Main", u"File", None))
        self.menuSave_as.setTitle(QCoreApplication.translate("Main", u"Save as...", None))
        self.menuOptions.setTitle(QCoreApplication.translate("Main", u"Options", None))
        self.menuHelp.setTitle(QCoreApplication.translate("Main", u"Help", None))
        self.menuInformation.setTitle(QCoreApplication.translate("Main", u"Information", None))
        self.menuGive_us_feedback.setTitle(QCoreApplication.translate("Main", u"Give us feedback", None))
    # retranslateUi

