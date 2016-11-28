# coding: utf-8
from PyQt5 import QtCore
from PyQt5.QtWebKitWidgets import QWebView


class MyBrowser(QWebView):
    #closing = QtCore.Signal()

    def __init(self):
        super().__init__()
        self.loadFinished.connec(self._results_available)

    def _results_available(self, ok):
        frame = self.page().mainFrame()
        print(unicode(frame.toHtml()).encode('utf-8'))

    def closeEvent(self, event):
        #print("done")
        pass
        #self.closing.emit()


