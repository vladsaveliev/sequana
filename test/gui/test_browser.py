from sequana.gui.browser import Browser
import pytest
from pytestqt.qt_compat import qt_api


def test_browser(qtbot):
    widget = Browser(url="sequana.readthedocs.io")
    widget.show()
    qtbot.addWidget(widget)
    #widget.wb.createWindow(None)
