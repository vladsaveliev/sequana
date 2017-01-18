import sys
from sequana.gui import about
from PyQt5 import QtWidgets as QW
from PyQt5 import Qt

#app = QW.QApplication(sys.argv)

import pytest
from pytestqt.qt_compat import qt_api


def test_directory_dialog(qtbot):
    #assert qt_api.QApplication.instance() is not None
    widget = about.About()
    widget.show()
    qtbot.addWidget(widget)

