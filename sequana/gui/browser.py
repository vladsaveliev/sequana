# coding: utf-8
from PyQt5 import QtCore, QtGui, Qt, QtWidgets
from PyQt5.QtWebKitWidgets import QWebView
from PyQt5 import QtWebKit
from PyQt5.QtWidgets import QProgressBar, QLineEdit


# potential resources for improvements:
# https://github.com/ralsina/devicenzo/blob/master/devicenzo.py

class Browser(Qt.QMainWindow):
    """


    On purpose, there is no caching so that (if re-generated), the 
    new content of an HTML is shown.
    
    """

    def __init__(self, url):
        Qt.QMainWindow.__init__(self)
        self.sb = self.statusBar()

        # Progress bar
        self.pbar = QProgressBar()
        self.pbar.setMaximumWidth(120)

        # Main page
        self.wb = QWebView(
            loadProgress=self.pbar.setValue,
            loadFinished=self.pbar.hide, 
            loadStarted=self.pbar.show, 
            titleChanged=self.setWindowTitle)
        self.setCentralWidget(self.wb)

        # Main menu
        self.tb = self.addToolBar("Main Toolbar")
        for a in (Qt.QWebPage.Back, Qt.QWebPage.Forward, Qt.QWebPage.Reload):
            self.tb.addAction(self.wb.pageAction(a))
        self.tb.addAction

        self.url = QLineEdit(returnPressed =lambda:self.wb.setUrl(
            QtCore.QUrl.fromUserInput(self.url.text())))
        self.tb.addWidget(self.url)

        self.wb.urlChanged.connect(lambda u: self.url.setText(u.toString()))

        self.wb.statusBarMessage.connect(self.sb.showMessage)
        self.wb.page().linkHovered.connect(lambda l: self.sb.showMessage(l, 3000))

        # Search bar
        self.search = QLineEdit(returnPressed = lambda: self.wb.findText(self.search.text()))
        self.search.show()
        self.search.hide() # To make ctrl+F effective, need to show/hide ?

        # THe shortcuts
        self.showSearch = Qt.QShortcut("Ctrl+F", self, 
            activated= lambda: (self.search.show() , self.search.setFocus()))
        self.hideSearch = Qt.QShortcut("Esc", self, 
            activated= lambda: (self.search.hide(), self.wb.setFocus()))
        self.quit = Qt.QShortcut("Ctrl+Q", self, activated= self.close)
        self.zoomIn = Qt.QShortcut("Ctrl++", self, 
            activated= lambda: self.wb.setZoomFactor(self.wb.zoomFactor()+.2))
        self.zoomOut = Qt.QShortcut("Ctrl+-", self, 
            activated= lambda: self.wb.setZoomFactor(self.wb.zoomFactor()-.2))
        self.zoomOne = Qt.QShortcut("Ctrl+=", 
            self, activated = lambda: self.wb.setZoomFactor(1))

        self.wb.settings().setAttribute(QtWebKit.QWebSettings.PluginsEnabled, True)

        self.sb.addPermanentWidget(self.search)
        self.sb.addPermanentWidget(self.pbar)
        self.wb.load(QtCore.QUrl(url))

        # No caching
        self.wb.settings().setObjectCacheCapacities(0,0,0)

