import sys
import os
import unittest

from sequana.gui import file_browser
from sequana import sequana_data
from PyQt5 import QtWidgets as QW
from PyQt5.QtTest import QTest

app = QW.QApplication(sys.argv)


class FileBrowserTest(unittest.TestCase):

    def setUp(self):
        self.form = file_browser.FileBrowser(paired=False, directory=False,
            file_filter=None)

    def tearDown(self):
        self.form.close()

    def test_paths(self):
        self.form.set_empty_path()
        self.form.get_filenames()

    def _test_paired(self):
        self.form.browse_directory()

import pytest
from pytestqt.qt_compat import qt_api
"""
#def test_directory_dialog(qtbot):
    #assert qt_api.QApplication.instance() is not None
    #widget = qt_api.QWidget()
    #qtbot.addWidget(widget)
    #widget.show()
    #assert widget.isVisible()
    def setUp(self):
        self.form = file_browser.DirectoryDialog(None, "", ".", "")

    def tearDown(self):
        self.form.close()

    def test_show(self):
        self.form.show()
        self.form.get_directory_path()
        # w1 = self.form.findChildren(PyQt5.QtWidgets.QPushButton)[0]
        # w1.click()
"""



if __name__ == "__main__":
    unittest.main()
