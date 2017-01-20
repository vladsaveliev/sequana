import sys
import os
import unittest

from sequana.gui.file_browser import FileBrowser, DirectoryDialog

from sequana import sequana_data
from PyQt5 import QtCore

import pytest
#from pytestqt.qt_compat import qt_api

#app = QW.QApplication(sys.argv)


# How to use that ? 
# http://pytest-qt.readthedocs.io/en/latest/tutorial.html
def test_basic_search(qtbot, tmpdir):
    '''
    test to ensure basic find files functionality is working.
    '''
    tmpdir.join('test1_R1.fastq.gz').ensure()
    tmpdir.join('test1_R2.fastq.gz').ensure()

    tmpdir.join('test2_R1.fastq.gz').ensure()
    tmpdir.join('test2_R2.fastq.gz').ensure()


def test_directory_dialog(qtbot, mock):

    widget = FileBrowser(paired=False, directory=False, file_filter=None)
    qtbot.addWidget(widget)
    widget.show()
    assert widget.isVisible()
    widget.setup_color()

    widget.set_empty_path()
    assert widget.get_filenames() == ""
    assert widget.path_is_setup() == False

    # Now, we open the dialog, which pops up. We need to close it...
    widget = FileBrowser(paired=True, directory=False, file_filter=None)
    qtbot.addWidget(widget)
    #qtbot.mouseClick(widget.btn, QtCore.Qt.LeftButton)

    widget._set_paired_filenames(["test1.fastq.gz"])
    widget._set_paired_filenames(["test1.fastq.gz", "test2.fastq.gz"])



def test_directory_dialog(qtbot, tmpdir):
    widget = DirectoryDialog(None, "test", str(tmpdir), "*gz")
    qtbot.addWidget(widget)
