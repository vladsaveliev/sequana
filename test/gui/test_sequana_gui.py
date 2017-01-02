import sys
import os
import unittest

from sequana.gui import sequana_gui
from sequana import sequana_data
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


class SequanaGUIStandaloneTest(unittest.TestCase):

    def setUp(self):
        from tempfile import TemporaryDirectory
        wkdir = TemporaryDirectory()
        inputdir = os.path.realpath(
                sequana_data("Hm2_GTGAAA_L005_R1_001.fastq.gz")).rsplit(os.sep,1)[0]

        options = {"pipeline": "quality_control", 
                    "--wkdir": wkdir.name, 
                    "input_directory": inputdir}
        from easydev import AttrDict
        options = AttrDict(**options)
        self.form = sequana_gui.SequanaGUI(ipython=False, user_options=options)

    def test(self):
        pass


if __name__ == "__main__":
    unittest.main()
