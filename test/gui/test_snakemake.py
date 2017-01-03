import sys
import os
import unittest

from sequana.gui import snakemake
from sequana import sequana_data
from PyQt5 import QtWidgets as QW
from PyQt5.QtTest import QTest
from PyQt5.QtCore import Qt

app = QW.QApplication(sys.argv)


class SnakemakeTest(unittest.TestCase):

    def setUp(self):
        self.form = snakemake.SnakemakeDialog()

    def tearDown(self):
        self.form.close()

    def test_get_options(self):
        self.form.get_snakemake_local_options()
        self.form.get_snakemake_general_options()
        self.form.get_snakemake_cluster_options()

    def test_accept(self):
        # important to not change the settings in the test
        self.form.accept()

    def test_reject(self):
        # important to not change the settings in the test
        self.form.reject()
    

    def test_generic_pipeline(self):
        self.form.get_settings()




if __name__ == "__main__":
    unittest.main()
