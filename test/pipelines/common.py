import subprocess
import tempfile
import shutil
import os

class Pipeline(object):
    """Base class for all pipeline tests""" 
    def __init__(self, wk=None):
        if wk:
            self.wk = wk
        else:
            self.tempdir = tempfile.TemporaryDirectory()
            self.wk = self.tempdir.name

    def run(self):
        res = subprocess.check_call(["sh", "runme.sh"], cwd=self.wk)
        if res != 0:
            raise Exception

    def clean(self):
        self.tempdir.cleanup()

