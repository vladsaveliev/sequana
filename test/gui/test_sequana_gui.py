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

    def _test_help(self):
        self.form.menuHelp()

    def _test_about(self):
        self.form.menuAbout()
        OK = self.form.preferences_dialog.ui.buttonBox.Ok
        okWidget = self.form.preferences_dialog.ui.buttonBox.button(OK)
        QTest.mouseClick(okWidget, Qt.LeftButton)

    def test_choose_pipeline(self):
        self.form.on_sequana_pipeline_choice("quality_control")
        assert self.form.pipeline_is_chosen == True

    def _test_generic_pipeline(self):
        from sequana import sequana_data
        self.snakefile = sequana_data("test_snakefile")

    def _test_config(self):
        self.form.ui.tabWidget_framework.setCurrentIndex(1)


if __name__ == "__main__":
    unittest.main()

