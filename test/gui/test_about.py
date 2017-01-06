import sys
from sequana.gui import about
from PyQt5 import QtWidgets as QW
from PyQt5 import Qt

app = QW.QApplication(sys.argv)

def test_about():
    w = about.About()
    e = Qt.QCloseEvent()
    w.event(e)
