# coding: utf-8
from PyQt5 import QtCore, QtGui, Qt, QtWidgets
try:
    # PyQt 5.6 and above (e.g. conda Mac PyQt5.6 does not work like PyQt5.6
    # linux) 
    from PyQt5.QtWebEngineWidgets import QWebEngineView as QWebView
    from PyQt5 import QtWebEngine as QtWebKit
    from PyQt5.Qt import QWebEnginePage as QWebPage
    from PyQt5.QtWebEngineWidgets import QWebEngineSettings
    from PyQt5.Qt import QTabWidget
except:
    # PyQt 5.6 and below
    from PyQt5.QtWebKitWidgets import QWebView
    from PyQt5 import QtWebKit
    from PyQt5.Qt import QWebPage
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

        # Progress bar
        # ------------------------------------------------------------
        self.pbar = QProgressBar()
        self.pbar.setMaximumWidth(120)


        # Main page QWebView
        # -------------------------------------------------------------
        self.wb = QWebView(
            loadProgress=self.pbar.setValue,
            loadFinished=self.pbar.hide,
            loadStarted=self.pbar.show,
            titleChanged=self.setWindowTitle)
        self.wb.urlChanged.connect(lambda u: self.url.setText(u.toString()))

        self.setCentralWidget(self.wb)

        # Main menu tool bar
        # -------------------------------------------------------------
        self.tb = self.addToolBar("Main Toolbar")
        for a in (QWebPage.Back, QWebPage.Forward, QWebPage.Reload):
            self.tb.addAction(self.wb.pageAction(a))

        self.url = QLineEdit(returnPressed =lambda:self.wb.setUrl(
            QtCore.QUrl.fromUserInput(self.url.text())))
        self.tb.addWidget(self.url)

        # status bar ---------------------------------------------------
        self.sb = self.statusBar()
        try:
            #pyqt5.6
            self.wb.statusBarMessage.connect(self.sb.showMessage)
        except:
            pass
        self.wb.page().linkHovered.connect(lambda l: self.sb.showMessage(l, 3000))

        # Search bar
        # ------------------------------------------------------------
        self.search = QLineEdit(returnPressed = lambda: self.wb.findText(self.search.text()))
        self.search.show()
        self.search.hide() # To make ctrl+F effective, need to show/hide ?

        # The shortcuts
        # ---------------------------------------------------------
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

        # Add alt+left and right keys to navigate backward and forward
        Qt.QShortcut(QtCore.Qt.AltModifier + QtCore.Qt.Key_Left, self,
            activated=lambda: self.wb.back())
        Qt.QShortcut(QtCore.Qt.AltModifier + QtCore.Qt.Key_Right, self,
            activated=lambda: self.wb.forward())

        # Javascript and other settings
        # ------------------------------------------------------------
        try:
            # pyqt 5.6 with WebKit version
            self.wb.settings().setAttribute(QtWebKit.QWebSettings.PluginsEnabled,
                                        True)
        except:
            # Qt with Webengine version
            self.wb.settings().setAttribute(
                QWebEngineSettings.JavascriptEnabled,  True)
            self.wb.settings().setAttribute(
                QWebEngineSettings.Accelerated2dCanvasEnabled, False)

        # Add components on the page
        self.sb.addPermanentWidget(self.search)
        self.sb.addPermanentWidget(self.pbar)

        # Finally, load the URL
        self.wb.load(QtCore.QUrl(url))

        # No caching
        # self.wb.settings().setObjectCacheCapacities(0,0,0)



#if __name__ == "__main__": 
#
#    b = Browser()
