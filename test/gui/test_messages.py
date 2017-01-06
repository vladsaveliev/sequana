import sys
from sequana.gui import messages
from PyQt5 import QtWidgets as QW

app = QW.QApplication(sys.argv)

def test_warning():
    w = messages.WarningMessage("test")

def test_info():
    w = messages.InfoMessage("test")

def test_critical():
    w = messages.CriticalMessage("test", details="test")




