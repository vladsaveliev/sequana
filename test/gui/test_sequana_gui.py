import sys
import unittest

from sequana.gui import sequana_gui

from PyQt5 import QtWidgets as QW
from PyQt5.QtTest import QTest
from PyQt5.QtCore import Qt

app = QW.QApplication(sys.argv)

class SequanaGUITest(unittest.TestCase):

    def setUp(self):
        self.form = sequana_gui.SequanaGUI(ipython=False)

    def tearDown(self):
        self.form.close()

    def test_help(self):
        pass
        # blocking need a way to close it
        # self.form.help() 
        #QTest.mouseClick(okWidget, Qt.LeftButton)

    def test_choose_pipeline(self):
        self.form.on_pipeline_choice("quality_control")
        assert self.form.pipeline_is_chosen == True



if __name__ == "__main__":
    unittest.main()

