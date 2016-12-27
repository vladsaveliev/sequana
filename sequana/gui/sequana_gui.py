#  This file is part of Sequana software
#
#  Copyright (c) 2016 - Sequana Development Team
#
#  File author(s):
#      Thomas Cokelaer <thomas.cokelaer@pasteur.fr>
#      Dimitri Desvillechabrol <dimitri.desvillechabrol@pasteur.fr>,
#          <d.desvillechabrol@gmail.com>
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"""Sequana GUI. Can also be used for any snakemake pipeline


"""
import sys
import os
import re
import time
import subprocess as sp
import shutil

from PyQt5.QtWidgets import QWidget
from ui_mainwindow import Ui_MainWindow
from ui_preferences import Ui_Preferences
from sequana.gui.snakemake import SnakemakeDialog

from PyQt5 import QtCore, QtGui
from PyQt5 import QtWidgets as QW
from PyQt5.Qt import QTemporaryDir, QMainWindow
from PyQt5.QtCore import Qt
from PyQt5.QtSvg import QSvgWidget
from PyQt5.QtWebKitWidgets import QWebView

from sequana import snaketools, sequana_data
from sequana.snaketools import Module
from sequana.gui.browser import MyBrowser
from sequana.gui.ipython import QIPythonWidget
from sequana.gui.about import About
from sequana.gui.file_browser import FileBrowser
from sequana.gui.widgets import *
from sequana.gui.messages import *


import signal
def sigint_handler(*args):
    """Handler for the SIGINT signal."""
    sys.stderr.write('\r')
    if QW.QMessageBox.question(None, '', "Are you sure you want to quit?",
                            QW.QMessageBox.Yes | QW.QMessageBox.No,
                            QW.QMessageBox.No) == QW.QMessageBox.Yes:
        QW.QApplication.quit()


class SequanaGUI(QMainWindow):
    """


    - save_config_file does not seem to work when we press it manually
    - settings in preferences and snakemake options
    - import config file
    - Generic pipelines
    - do not copy again requirements if already there
    - extension of the different widgets
    - cursor in the snakemake output must follow the progress and be at the
      bottom all the time


    Developer Guide
    ------------------

    - The GUI is designed with qt designer as much as possible.
    - All GUI objects are in the **ui** attributes. Additional dialog such as the
      snakemake and preferences dialog have their own modules and stored in attributes
      ending in _dialog

    """


    _not_a_rule = {"requirements", "gatk_bin", "input_directory", "input_samples", "input_pattern"}
    _browser_keyword = {"reference"}

    def __init__(self):
        super(SequanaGUI, self).__init__()

        self._tempdir = QTemporaryDir()
        self.shell = ""
        self.shell_error = ""
        self._colors = {
            'green': QtGui.QColor(0,170,0),
            'red': QtGui.QColor(170,0,0),
            'blue': QtGui.QColor(0,90,154),
        }

        # some global attributes
        self.sequana_config = None
        self.pipeline_is_chosen = False

        self.initUI()

    def initUI(self):

        # Set up the user interface from Designer. This is the general layout
        # without dedicated widgets and connections
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        # Now we add what could not be added in the designer
        # First the preferences and snakemake dialogs
        self.preferences_dialog = PreferencesDialog()
        self.snakemake_dialog = SnakemakeDialog()

        # The IPython dialog, which is very useful for debugging
        self.ipyConsole = QIPythonWidget(
            customBanner="Welcome to the embedded ipython console\n")
        self.ipyConsole.pushVariables({"x": 10})
        self.ipyConsole.printText("The variable 'foo' andion.")
        self.ipyConsole.execute("from sequana import *")
        self.ipyConsole.execute("import sequana")
        self.ipyConsole.execute("")
        self.ipyConsole.pushVariables({"gui": self})
        self.ui.layout_ipython.addWidget(self.ipyConsole)

        # layout for config file parameters
        widget_form = QW.QWidget()
        self.form = QW.QVBoxLayout(widget_form)
        self.form.setSpacing(10)
        self.ui.scrollArea.setWidget(widget_form)
        self.ui.scrollArea.setWidgetResizable(True)
        self.ui.scrollArea.setMinimumHeight(300)

        # layout for the snakemake output
        self.output = QW.QTextEdit()
        self.ui.layout_snakemake.addWidget(self.output)

        self.logger = QW.QTextEdit()
        self.ui.layout_logger.addWidget(self.logger)

        # A file browser for the working directory layout
        self.working_dir = FileBrowser(directory=True)
        self.ui.layout_wk.addWidget(self.working_dir)
        self.working_dir.clicked_connect(self.check_existing_config)
        self.working_dir.clicked_connect(self.switch_run)

        # Connectors to actions related to the menu bar
        self.ui.actionQuit.triggered.connect(self.menuQuit)
        self.ui.actionHelp.triggered.connect(self.menuHelp)
        self.ui.actionAbout.triggered.connect(self.menuAbout)
        self.ui.actionImportConfig.triggered.connect(self.menuImportConfig)
        self.ui.actionSnakemake.triggered.connect(self.snakemake_dialog.exec_)
        self.ui.actionPreferences.triggered.connect(self.preferences_dialog.exec_)

        self.set_sequana_pipeline()

        # The run/save/dag footer buttons
        self.connect_footer_buttons()

        self.process = QtCore.QProcess(self)
        self.process.started.connect(lambda: self.ui.run_btn.setEnabled(False))

        self.process.started.connect(lambda: self.ui.stop_btn.setEnabled(True))
        self.process.started.connect(lambda: self.start_progress)

        self.process.finished.connect(lambda: self.ui.run_btn.setEnabled(True))
        self.process.finished.connect(lambda: self.ui.stop_btn.setEnabled(False))
        self.process.finished.connect(lambda: self.end_progress)

        self.process.readyReadStandardOutput.connect(
            self.snakemake_data_stdout)
        self.process.readyReadStandardError.connect(self.snakemake_data_error)
        self.process.finished.connect(self.end_run)


    #|-----------------------------------------------------|
    #|                       MENU related                  |
    #|-----------------------------------------------------|
    def menuAbout(self):
        from sequana import version
        url = 'sequana.readthedocs.io'
        msg = About()
        msg.setIcon(QW.QMessageBox.Information)
        msg.setText("Sequana version %s " % version)
        msg.setInformativeText("""
            Online documentation on <a href="http://%(url)s">%(url)s</a>
            <br>
            <br>
            Authors: Thomas Cokelaer and Dimitri Desvillechabrol, 2016
            """ % {"url": url})
        msg.setWindowTitle("Sequana")
        msg.setStandardButtons(QW.QMessageBox.Ok)
        retval = msg.exec_()

    def menuHelp(self):
        url = 'sequana.readthedocs.io'
        msg = About()
        msg.setIcon(QW.QMessageBox.Information)
        msg.setText("Sequana GUI help")
        msg.setInformativeText("""<p>

            <ol>
            <li> Select a pipeline</li>
            <li> Select the directory or sample tab</li>
                <ul>
               <li> directory: select all fastq.gz files</li>
               <li> samples: select a single-end or paired-end file(s)</li>
                </ul>
            <li> Select the working directory</li>
            <li> Changet the options if needed</li>
            <li> Save the config file (check the DAG image)</li>
            <li> Run the analysis</li>
            <li> Open the report if successful</li>

            </p>
            """ % {"url": url})
        msg.setWindowTitle("Sequana")
        msg.setStandardButtons(QW.QMessageBox.Ok)
        self._msg_help = msg
        retval = msg.exec_()

    def menuQuit(self):
        quit_msg = WarningMessage("Do you really want to quit ?")
        quit_msg.setStandardButtons(QW.QMessageBox.Yes | QW.QMessageBox.No)
        quit_msg.setDefaultButton(QW.QMessageBox.No)
        quit_answer = quit_msg.exec_()
        if quit_answer == QW.QMessageBox.Yes:
            self.close()

    def menuImportConfig(self):
        if self.pipeline_is_chosen is False:
            msg = CriticalMessage("Please select a pipeline first")
            msg.exec_()
            return
        def _local_read_config():
            config_file = self.tab.get_filenames()
            if config_file:
                self._read_config(config_file)
        self.tab = FileBrowser(paired=False,
                               file_filter="YAML file (*.yaml)")
        self.tab.clicked_connect(_local_read_config)
        self.tab.clicked_connect(self.tab.close)
        self.tab.show()

    # ---------------------------------------------------------------
    #
    # ---------------------------------------------------------------

    def set_sequana_pipeline(self):
        # Fill the pipelines
        snaketools.pipeline_names.sort()
        self.ui.choice_button.addItems(["select pipeline"] +
            snaketools.pipeline_names)
        self.ui.choice_button.currentIndexChanged[str].connect(
            self.on_pipeline_choice)
        self.ui.choice_button.installEventFilter(self)

        # Set the file browser
        fastq_filter = "Fastq file (*.fastq *.fastq.gz *.fq *.fq.gz)"
        paired_tab = FileBrowser(paired=True, file_filter=fastq_filter)
        directory_tab = FileBrowser(directory=True)
        paired_tab.clicked_connect(self.switch_run)
        directory_tab.clicked_connect(self.switch_run)
        # fill the tabs_browser
        self.ui.tabs_browser.addTab(directory_tab, "Directory")
        self.ui.tabs_browser.addTab(paired_tab, "Sample")
        self.ui.tabs_browser.removeTab(0)
        self.ui.tabs_browser.removeTab(0)


    @QtCore.pyqtSlot(str)
    def on_pipeline_choice(self, index):
        """ Change options form when user change the pipeline.
        """
        config_file = snaketools.Module(index)._get_config()
        self.ui.choice_button.removeItem(
            self.ui.choice_button.findText("select pipeline"))

        self._read_config(config_file)
        self.fill_combobox(self.rule_list)
        self.pipeline_is_chosen = True
        self.switch_run()

    def fill_combobox(self, whatever):
        active_list = [w.get_name() for w in self.rule_list if w.get_do_rule()]

        self.ui.until_box.clear()
        self.ui.until_box.addItems([None] + active_list)

        self.ui.starting_box.clear()
        self.ui.starting_box.addItems([None] + active_list)

    def _read_config(self, config_file):
        try:
            self.sequana_config = snaketools.SequanaConfig(config_file)
            self.sequana_config.cleanup() # set all empty strings and %()s to None
        except AssertionError:
            print("Warning: could not parse the config file")
            return
        self.create_base_form()

    def import_config(self):
        if self.pipeline_is_chosen is False:
            msg = CriticalMessage("Please select a pipeline first")
            msg.exec_()
            return
        def _local_read_config():
            config_file = self.tab.get_filenames()
            if config_file:
                self._read_config(config_file)
        self.tab = FileBrowser(paired=False,
                               file_filter="YAML file (*.yaml)")
        self.tab.clicked_connect(_local_read_config)
        self.tab.clicked_connect(self.tab.close)
        self.tab.show()

    def clear_layout(self, layout):
        """ Clean all widget contained in a layout.
        """
        while layout.count():
            child = layout.takeAt(0)
            if child.widget() is not None:
                child.widget().deleteLater()
            elif child.layout() is not None:
                self.clear_layout(child.layout())

        # action import config file
        importAction = QW.QAction("Import", self)
        importAction.setShortcut('Ctrl+I')
        importAction.triggered.connect(self.import_config)

    @QtCore.pyqtSlot(str)
    def on_pipeline_choice(self, index):
        """ Change options form when user change the pipeline.
        """
        config_file = snaketools.Module(index)._get_config()
        self.ui.choice_button.removeItem(
            self.ui.choice_button.findText("select pipeline"))

        self._read_config(config_file)
        self.fill_combobox(self.rule_list)
        self.pipeline_is_chosen = True
        self.switch_run()

    def connect_footer_buttons(self):
        self.ui.run_btn.setEnabled(False)
        self.ui.run_btn.clicked.connect(self.start_sequana)

        self.ui.stop_btn.clicked.connect(self.click_stop)
        self.ui.stop_btn.setEnabled(False)

        self.ui.unlock_btn.setShortcut("Ctrl+U")
        self.ui.unlock_btn.clicked.connect(self.unlock_snakemake)
        self.ui.unlock_btn.setEnabled(True)

        self.ui.report_btn.setEnabled(True)
        self.ui.report_btn.clicked.connect(self.open_report)

        self.ui.save_btn = QW.QPushButton("Save Config")
        self.ui.save_btn.clicked.connect(self.save_config_file)

        self.ui.dag_btn.setEnabled(False)
        self.ui.dag_btn.clicked.connect(self.show_dag)

    def _log(self, text, color="blue"):
        print(text)
        cursor = self.logger.textCursor()
        cursor.movePosition(cursor.End)
        cursor.insertHtml('<p style="color:%s">%s</p><br>' % (color, text))
        cursor.movePosition(cursor.End)

    def snakemake_data_stdout(self):
        """ Read standard output of snakemake process
        """
        cursor = self.output.textCursor()
        cursor.movePosition(cursor.End)
        data = str(self.process.readAllStandardOutput())
        self.shell += data
        self.update_progress_bar(data)

        for this in data.split("\\n"):
            line = this.strip()
            if line and len(line) > 3 and "complete in" not in line: # prevent all b'' strings
                line = line.replace("\\r","")
                line = line.replace("\\t","    ")
                cursor.insertHtml('<p style="color:blue">' + line +'</p><br>')
                cursor.movePosition(cursor.End)

    def snakemake_data_error(self):
        """ Read error output of snakemake process
        """
        cursor = self.output.textCursor()
        cursor.movePosition(cursor.End)
        error = str(self.process.readAllStandardError())
        self.shell_error += error
        self.update_progress_bar(error)

        for this in error.split("\\n"):
            line = this.strip()
            if line and len(line) > 3 and "complete in" not in line: # prevent all b'' strings
                line = line.replace("\\r","")
                line = line.replace("\\t","    ")
                cursor.insertHtml('<p style="color:red">' + line +'</p>')
                cursor.movePosition(cursor.End) 

    def get_until_starting_option(self):
        """ Return list with starting rule and end rule.
        """
        until_rule = self.ui.until_box.currentText()
        starting_rule = self.ui.starting_box.currentText()
        option = []
        if until_rule:
            option += ["--no-hooks", "-U", until_rule]
        if starting_rule:
            option += ["-R", starting_rule]
        return option

    def create_base_form(self):
        """ Create form with all options necessary for a pipeline.
        """
        self._log("Interpreting config file")
        self.clear_layout(self.form)
        rules_list = list(self.sequana_config._yaml_code.keys())
        rules_list.sort()
        self.necessary_dict = {}
        self.rule_list = []
        for count, rule in enumerate(rules_list):
            # Check if this is a dictionnary
            contains = self.sequana_config._yaml_code[rule]
            if isinstance(contains, dict) and (
                    rule not in SequanaGUI._not_a_rule):
                rule_box = Ruleform(rule, contains, count, self._browser_keyword)
                self.form.addWidget(rule_box)
                self.rule_list.append(rule_box)
                rule_box.connect_do(self.fill_combobox)
            else:
                if isinstance(contains, list):
                    self.necessary_dict = dict(self.necessary_dict,
                                           **{rule: contains})
                elif contains is None or contains in ["''", '""']:
                    self.necessary_dict = dict(self.necessary_dict,
                                           **{rule: None})
                else:
                    self.necessary_dict = dict(self.necessary_dict,
                                           **{rule: '{0}'.format(contains)})

    def click_stop(self):
        print("stopped manually")
        pal = self.ui.progressBar.palette()
        pal.setColor(QtGui.QPalette.Highlight, self._colors['red'])
        self.ui.progressBar.setPalette(pal)
        self.process.finished.disconnect() # we do not want to call end_run()
        self.ui.run_btn.setEnabled(True)
        self.ui.stop_btn.setEnabled(False)
        self.process.close()

    def update_progress_bar(self, line):
        """ Parse with a regex to retrieve current step and total step.
        """
        grouprex = self._step_regex.findall(line)
        if grouprex:
            step = int(grouprex[0][0]) / float(grouprex[0][1]) * 100
            self.ui.progressBar.setValue(step)

    def end_run(self):
        pal = self.ui.progressBar.palette()
        if self.ui.progressBar.value() >= 100 :
            pal.setColor(QtGui.QPalette.Highlight, self._colors['green'])
            self.ui.progressBar.setPalette(pal)
            print('Run done. Status: successful')
        else:
            pal.setColor(QtGui.QPalette.Highlight, self._colors['red'])
            self.ui.progressBar.setPalette(pal)
            text = 'Run manually to check the exact error or check the log.'
            if "--unlock" in self.shell_error:
                text += "<br>You may need to unlock the directory. "
                text += "click on Unlock button"
            msg = CriticalMessage(text, self.process.readAllStandardError())
            msg.exec_()
            return

    def start_progress(self):
        self.ui.progressBar.setRange(0,1)

    def end_progress(self):
        self.ui.progressBar.setValue(100)
        QtGui.QMessageBox.information(self, "Done")
        self.ui.run_btn.setEnabled(True)

    def save_config_file(self):
        print("hello")
        try:
            form_dict = dict(self.create_form_dict(self.form),
                                 **self.necessary_dict)
        except AttributeError:
            msg = WarningMessage("You must choose a pipeline before saving.")
            msg.exec_()
            return

        # get samples names or input_directory
        # TODO
        if self.ui.tabs_browser.currentIndex() == 1:
            form_dict["samples"] = (
                self.ui.tabs_browser.currentWidget().get_filenames())
        else:
            form_dict["input_directory"] = (
                self.ui.tabs_browser.currentWidget().get_filenames())

        # Let us update tha attribute with the content of the form
        self.sequana_config._yaml_code.update(form_dict)

        if self.working_dir.path_is_setup():
            print(1)
            yaml_path = self.working_dir.get_filenames() + "/config.yaml"
            self.sequana_config.copy_requirements(target=self.working_dir.get_filenames())
            print(yaml_path)

            if os.path.isfile(yaml_path):
                save_msg = WarningMessage(
                    "The file {0} already exist".format(yaml_path))
                save_msg.setInformativeText(
                    "Do you want to overwrite the file?")
                save_msg.setStandardButtons(
                    QW.QMessageBox.Yes | QW.QMessageBox.Discard |
                    QW.QMessageBox.Cancel)
                save_msg.setDefaultButton(QW.QMessageBox.Yes)
                # Yes == 16384
                # Save == 2048
                if save_msg.exec_() in [16384, 2048]:
                    self.sequana_config.cleanup()
                    self.sequana_config.save(yaml_path)
                    self.ui.dag_btn.setEnabled(True)
            else:
                self.sequana_config.cleanup()
                self.sequana_config.save(yaml_path)
                self.ui.dag_btn.setEnabled(True)
        else:
            msg = WarningMessage("You must set a working directory", self)
            msg.exec_()

    def _get_snakemake_command(self, snakefile):
        dialog = self.snakemake_dialog      # an alias
        snakemake_line = ["-s", snakefile, "--stat", "stats.txt"]

        if self.ui.comboBox_local.currentText() == "local":
            snakemake_line += dialog.get_snakemake_local_options()
        elif self.ui.comboBox_local.currentText() == "cluster":
            snakemake_line += dialog.get_snakemake_cluster_options()

        snakemake_line += dialog.get_snakemake_general_options()

        snakemake_line += self.get_until_starting_option()

        return snakemake_line

    def start_sequana2(self):
        # Is it still working ? Can we switch this option from the GUI ?
        working_dir = self.working_dir.get_filenames()
        self.save_config_file()
        rules = self.get_rules(self.form)
        snakefile = Module(self.ui.choice_button.currentText()).snakefile

        new_snakefile = working_dir + os.sep + os.path.basename(snakefile)
        if os.path.exists(new_snakefile) is False:
            shutil.copy(snakefile, working_dir)
        snakefile = new_snakefile

        snakemake_args = self._get_snakemake_command(snakefile)

        self.cmd = ['snakemake'] + snakemake_args

        #self.
        self.cwd = working_dir
        snakemake_proc = sp.Popen(self.cmd,
            cwd=os.path.basename(working_dir))
        #snakemake_proc.communicate()


    def start_sequana(self):
        pal = self.ui.progressBar.palette()
        pal.setColor(QtGui.QPalette.Highlight, self._colors['blue'])
        self.ui.progressBar.setPalette(pal)
        self.ui.progressBar.setValue(1)

        # Set the regex to catch steps
        self._step_regex = re.compile("([0-9]+) of ([0-9]+) steps")

        # Prepare the command and working directory.
        working_dir = self.working_dir.get_filenames()
        self.save_config_file()
        rules = self.get_rules(self.form)
        snakefile = Module(self.ui.choice_button.currentText()).snakefile

        new_snakefile = working_dir + os.sep + os.path.basename(snakefile)
        if os.path.exists(new_snakefile) is False:
            shutil.copy(snakefile, working_dir)
        snakefile = new_snakefile

        snakemake_args = self._get_snakemake_command(snakefile)

        self._log("Starting process with %s " % snakemake_args)
        self.process.setWorkingDirectory(working_dir)
        self.process.start("snakemake", snakemake_args)

    def switch_run(self):
        if self.working_dir.path_is_setup():
            #TODO
            if self.ui.tabs_browser.currentWidget().path_is_setup():
                if self.working_dir.path_is_setup():
                    if self.pipeline_is_chosen:
                        return self.ui.run_btn.setEnabled(True)
        return self.ui.run_btn.setEnabled(False)


    #@requires("working_dir")
    def unlock_snakemake(self):

        if self.working_dir.get_filenames() == "":
            return

        working_dir = self.working_dir.get_filenames()

        # TODO: this is sequana specific
        snakefile = Module(self.ui.choice_button.currentText()).snakefile
        new_snakefile = working_dir + os.sep + os.path.basename(snakefile)
        if os.path.exists(new_snakefile) is False:
            print("config not found. should not happen")
        else:
            snakefile = new_snakefile
            self.cmd = ['snakemake', "-s", snakefile, "--unlock"]
            self.cwd = working_dir
            print(self.cmd)
            print("Please wait a second")
            print(working_dir)
            snakemake_proc = sp.Popen(self.cmd,
                cwd=working_dir)
#                cwd=os.path.basename(working_dir))
        self.ui.run_btn.setEnabled(True)
        self.ui.stop_btn.setEnabled(False)
        #snakemake_proc.communicate()

    def show_dag(self):
        if self.pipeline_is_chosen:
            snakefile = Module(self.ui.choice_button.currentText()).snakefile
        else:
            msg = CriticalMessage("Please select a pipeline first")
            msg.exec_()
            return

        if self.sequana_config is None:
            msg = CriticalMessage("Please select the input and working directory")
            msg.exec_()
            return

        # The config must have been saved, so we just need to copy it
        working_dir = self.working_dir.get_filenames()

        # copy the pipeline in the working directory
        shutil.copy(snakefile, working_dir)
        snakefile = os.path.basename(snakefile)

        svg_filename = self._tempdir.path() + os.sep + "test.svg"

        snakemake_line = ["snakemake", "-s", snakefile]
        snakemake_line += ["--rulegraph", "--configfile", "config.yaml"]
        snakemake_line += self.get_until_starting_option()

        self.process1 = QtCore.QProcess(self)
        self.process2 = QtCore.QProcess(self)
        self.process1.setWorkingDirectory(working_dir)
        self.process1.setStandardOutputProcess(self.process2)

        self.process1.start("snakemake", snakemake_line[1:])
        self.process2.start("dot", ["-Tsvg", "-o", svg_filename])

        self.process1.waitForFinished(10000)
        self.process2.waitForFinished(10000)

        if os.path.exists(svg_filename):
            diag = SVGDialog(svg_filename)
            diag.exec_()
        else:
            msg = 'Could not create the DAG file.'
            error = str(self.process1.readAllStandardError())
            msg = CriticalMessage(msg, error)
            msg.exec_()
            return

    def open_report(self):
        dialog = self.preferences_dialog.ui # an alias
        if self.pipeline_is_chosen and self.working_dir.get_filenames():
            filename = self.working_dir.get_filenames() + "/multi_summary.html"
            if os.path.exists(filename) is False:
                WarningMessage("""multi_summary.html not found.
                Most probably the analysis did not finish correctly""")
                return
            url = "file://" + filename
            browser = dialog.preferences_options_browser.currentText()

            if browser == "pyqt5":
                self.browser = MyBrowser()
                self.browser.load(QtCore.QUrl(url))
                self.browser.show()
            else:
                from easydev import execute
                try:
                    execute("%s %s" % (browser, url), showcmd=False,
                        verbose=False)
                except:
                    msg = CriticalMessage("Error browser",
                        "Could not open %s with %s. Is %s on your system or path ?" % 
                        (url, browser,browser))
                    msg.exec_()
                    return
        else:
            msg = WarningMessage("no working directory selected yet")
            msg.exec_()

    def create_form_dict(self, layout):
        def _cleaner(value):
            # This is to save the YAML file correctly since the widgets tend to
            # convert None and empty strings as '""' or "''"
            if value in ['None', None, '', '""', "''"]:
                return None
            else:
                return value

        widgets = (layout.itemAt(i).widget() for i in range(layout.count()))
        form_dict = {w.get_name(): _cleaner(w.get_value()) if w.is_option()
                         else self.create_form_dict(w.get_layout())
                         for w in widgets}
        return form_dict

    def check_existing_config(self):
        config_file = self.working_dir.get_filenames() + "/config.yaml"
        if not os.path.isfile(config_file):
            return False

        if not self.pipeline_is_chosen:
            msg = WarningMessage("A config.yaml file already exist in this "
                                 "directory. Please, choose a pipeline to "
                                 "know if the existing config file correspond "
                                 "to your pipeline.", self)
            self.working_dir.set_empty_path()
            msg.exec_()
            return False

        try:
            # Load the existing config
            cfg = snaketools.SequanaConfig(config_file)
            cfg.cleanup() # set all empty strings and %()s to None
            config_dict = cfg.config
        except AssertionError:
            msg = WarningMessage("Could not parse the config file.")
            msg.exec_()
            return False
        except Exception as err:
            msg = WarningMessage(
                "Unexpected error while checking the config file %s."
                % config_file + str(err))
            msg.exec_()
            return False

        if set(self.sequana_config._yaml_code.keys()) == set(config_dict.keys()):
            msg = QW.QMessageBox(
                QW.QMessageBox.Question, "Question",
                "A config file already exist in the working directory.\n"
                "Do you want to import its content ?",
                QW.QMessageBox.Yes | QW.QMessageBox.No,
                self, Qt.Dialog | Qt.CustomizeWindowHint)
            if msg.exec_() == 16384:
                self.sequana_config._yaml_code.update(config_dict)
                self.create_base_form()
                self.fill_combobox(self.rule_list)
        return True


    def get_rules(self, layout):
        widgets = (layout.itemAt(i).widget() for i in range(layout.count()))
        rules = [w.get_name() for w in widgets]
        return rules

    def eventFilter(self, source, event):
        """ Inactivate wheel event of combobox
        """
        if event.type() == QtCore.QEvent.Wheel and source is self.ui.choice_button:
            return True
        return False

    def closeEvent(self, event):
        #Close button (red X)
        self._tempdir.remove()
        try:self.browser.close()
        except:pass

    def close(self):
        # Menu or ctrl+q
        self._tempdir.remove()
        try:self.browser.close()
        except:pass
        print("bye now.")
        super().close()


class PreferencesDialog(QW.QDialog):
    """FIXME May not be required anymore"""
    def __init__(self, parent=None):
        super().__init__(parent=parent)
        self.ui = Ui_Preferences()
        self.ui.setupUi(self)


class SVGDialog(QW.QDialog):
    def __init__(self, filename):
        super().__init__()
        self.main_layout = QW.QVBoxLayout(self)
        self.setWindowTitle("DAG")

        if os.path.exists(filename):
            widget = QSvgWidget(filename)
            self.main_layout.addWidget(widget)


def main():
    signal.signal(signal.SIGINT, sigint_handler)
    app = QW.QApplication(sys.argv)

    filename = sequana_data("drawing.png", "../gui")

    if "--nosplash" in sys.argv:
        app.processEvents()
        sequana = SequanaGUI()
        sequana.show()
        if "--nosplash" in sys.argv is False:
            splash.finish(sequana)
        sys.exit(app.exec_())
    else:
        # Show the splash screen for a few seconds
        splash_pix = QtGui.QPixmap(filename)
        splash = QW.QSplashScreen(splash_pix, Qt.WindowStaysOnTopHint)
        #progressBar = QW.QProgressBar(splash)
        splash.setMask(splash_pix.mask())
        splash.show()
        #progressBar.show()
        """for i in range(0, 100):
            #progressBar.setValue(i)
            t = time.time()
            while time.time() < t + 2./100.:
                app.processEvents()
        """
        app.processEvents()
        sequana = SequanaGUI()
        sequana.show()
        splash.finish(sequana)
        sys.exit(app.exec_())


if __name__ == "__main__":
    main()

