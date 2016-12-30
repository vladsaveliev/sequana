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
"""Sequana GUI. Can also be used for any snakemake pipeline"""
import sys
import os
import re
import time
import subprocess as sp
import shutil

from sequana.gui.ui_mainwindow import Ui_MainWindow
from sequana.gui.browser import MyBrowser
from sequana.gui.ipython import QIPythonWidget
from sequana.gui.about import About
from sequana.gui.file_browser import FileBrowser
from sequana.gui.widgets import Ruleform
from sequana.gui.messages import *
from sequana.gui.preferences import PreferencesDialog
from sequana.gui.snakemake import SnakemakeDialog

from PyQt5 import QtCore, QtGui
from PyQt5 import QtWidgets as QW
from PyQt5.Qt import QTemporaryDir, QMainWindow
from PyQt5.QtCore import Qt
from PyQt5.QtSvg import QSvgWidget

from sequana import snaketools, sequana_data
from sequana.snaketools import Module

import colorlog
import logging

import signal

END = QtGui.QTextCursor.End

def sigint_handler(*args):
    """Handler for the SIGINT signal."""
    sys.stderr.write('\r')
    if QW.QMessageBox.question(None, '', "Are you sure you want to quit?",
                            QW.QMessageBox.Yes | QW.QMessageBox.No,
                            QW.QMessageBox.No) == QW.QMessageBox.Yes:
        QW.QApplication.quit()


class QPlainTextEditLogger(colorlog.StreamHandler):
    def __init__(self, parent):
        super().__init__()
        self.widget = QW.QPlainTextEdit(parent)
        self.widget.setReadOnly(True)

    def emit(self, record):
        msg = self.format(record)
        self.widget.appendPlainText(msg)
        self.widget.moveCursor(END)

class SequanaGUI(QMainWindow):
    """

    If quiet, progress bar cannot work.

    - Generic pipelines
    - do not copy again requirements if already there
    - extension of the different widgets ?
    - cursor in the snakemake output must follow the progress and be at the
      bottom all the time. Same for the logger


    Developer Guide
    ------------------

    - The GUI is designed with qt designer as much as possible.
    - All GUI objects are in the **ui** attributes. Additional dialog such as the
      snakemake and preferences dialog have their own modules and stored in attributes
      ending in _dialog

    """


    _not_a_rule = {"requirements", "gatk_bin", "input_directory", "input_samples", "input_pattern"}
    _browser_keyword = {"reference"}

    def __init__(self, parent=None, ipython=True):
        super(SequanaGUI, self).__init__(parent=parent)

        self._tempdir = QTemporaryDir()
        self.shell = ""
        self.shell_error = ""
        self._colors = {
            'green': QtGui.QColor(0,170,0),
            'red': QtGui.QColor(170,0,0),
            'blue': QtGui.QColor(0,90,154),
        }

        # some global attributes
        self.pipeline_is_chosen = False
        self._undefined_section = "Parameters in no sections/rules"

        self._ipython_tab = ipython
        self.initUI()
        self.read_settings()

        self.setStyleSheet("""QToolTip {
                           background-color: white;
                           color: black;
                           border-style: double;
                           border-width: 3px;
                           border-color: green;
                           margin:3px;
                           }""")

    def initUI(self):

        # Set up the user interface from Designer. This is the general layout
        # without dedicated widgets and connections
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        # 2 more dialogs from designer
        self.preferences_dialog = PreferencesDialog(self)
        self.snakemake_dialog = SnakemakeDialog(self)

        # The IPython dialog, which is very useful for debugging
        if self._ipython_tab is True:
            self.ipyConsole = QIPythonWidget(
                customBanner="Welcome to the embedded ipython console\n")
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
        self.output.setReadOnly(True)
        self.ui.layout_snakemake.addWidget(self.output)

        # Add the new logging box widget to the layout
        self.logTextBox = QPlainTextEditLogger(self)
        self.logTextBox.setFormatter(colorlog.ColoredFormatter(
            '%(log_color)s%(asctime)s - %(levelname)s - %(message)s',
            log_colors={
            'DEBUG': 'cyan',
            'INFO': 'green',
            'WARNING': 'yellow',
            'ERROR': 'red',
            'CRITICAL': 'red,bg_white',
            }
            ))
        colorlog.getLogger().addHandler(self.logTextBox)
        # You can control the logging level
        colorlog.getLogger().setLevel(logging.INFO)
        self.ui.layout_logger.addWidget(self.logTextBox.widget)

        # A file browser for the working directory layout
        self.working_dir = FileBrowser(directory=True)
        self.ui.layout_wk.addWidget(self.working_dir)
        self.working_dir.clicked_connect(self.check_existing_config)
        self.working_dir.clicked_connect(self.switch_run)

        # Connectors to actions related to the menu bar
        self.ui.actionQuit.triggered.connect(self.menuQuit)
        self.ui.actionHelp.triggered.connect(self.menuHelp)
        self.ui.actionAbout.triggered.connect(self.menuAbout)
        #self.ui.actionImportConfig.triggered.connect(self.menuImportConfig)
        self.ui.actionSnakemake.triggered.connect(self.snakemake_dialog.exec_)
        self.ui.actionPreferences.triggered.connect(self.preferences_dialog.exec_)

        self.set_sequana_pipeline()
        self.set_generic_pipeline()

        # The run/save/dag footer buttons
        self.connect_footer_buttons()

        self.process = QtCore.QProcess(self)
        self.process.started.connect(lambda: self.ui.run_btn.setEnabled(False))

        self.process.started.connect(lambda: self.ui.stop_btn.setEnabled(True))
        self.process.started.connect(lambda: self.start_progress)

        self.process.finished.connect(lambda: self.ui.run_btn.setEnabled(True))
        self.process.finished.connect(lambda: self.ui.stop_btn.setEnabled(False))

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
        self.msg = About()
        self.msg.setIcon(QW.QMessageBox.Information)
        self.msg.setText("Sequana version %s " % version)
        self.msg.setInformativeText("""
            Online documentation on <a href="http://%(url)s">%(url)s</a>
            <br>
            <br>
            Authors: Thomas Cokelaer and Dimitri Desvillechabrol, 2016
            """ % {"url": url})
        self.msg.setWindowTitle("Sequana")
        self.msg.setStandardButtons(QW.QMessageBox.Ok)
        retval = self.msg.exec_()

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


            <h2>Format of the config file</h2>
            For Sequana, config file are downloaded automatically so work out of
            the box. If you decide to change them (locally), every rule has only
            one level.

            For the generic file, only one level of heararchy is interpreted for
            now. If more than one level, the raw config file is provided.

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

    # ---------------------------------------------------------------
    # More GUI / reading the snakefile (sequana or generic)
    # ---------------------------------------------------------------
    def set_sequana_pipeline(self):
        # Fill the pipelines
        snaketools.pipeline_names.sort()
        self.ui.choice_button.addItems(snaketools.pipeline_names)
        self.ui.choice_button.currentIndexChanged[str].connect(
            self.on_sequana_pipeline_choice)
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
    def on_sequana_pipeline_choice(self, index):
        """ Change options form when user change the pipeline."""
        self.info("Reading sequana %s pipeline" % index)
        module = snaketools.Module(index)
        config_file = module._get_config()
        self.snakefile = module.snakefile
        self._read_sequana_config(config_file)
        self.fill_until_starting(self.rule_list)
        self.pipeline_is_chosen = True
        self.switch_run()

    def set_generic_pipeline(self):
        # Set the file browser
        self._generic_snakefile = FileBrowser(directory=False)
        self._generic_config = FileBrowser(directory=False,
            file_filter="YAML file (*.json *.yaml)")
        self._generic_snakefile.clicked_connect(self._read_generic_snakefile)
        # no switch run for the config file, which is optional
        self._generic_config.clicked_connect(self._read_generic_config)

        # fill the tabs_browser
        self.ui.tabs_browser_2.addTab(self._generic_snakefile, "Snakefile")
        self.ui.tabs_browser_2.addTab(self._generic_config, "Config file")
        self.ui.tabs_browser_2.removeTab(0)
        self.ui.tabs_browser_2.removeTab(0)

    def _read_generic_snakefile(self):
        filename = self._generic_snakefile.get_filenames()
        if len(filename):
            self.snakefile = filename
        else:
            self.snakefile = None
        self.switch_run()

    # ---------------------------------------------------------------------
    # Fotter connectors
    # ---------------------------------------------------------------------

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

        self.ui.save_btn.clicked.connect(self.save_config_file)

        self.ui.dag_btn.setEnabled(False)
        self.ui.dag_btn.clicked.connect(self.show_dag)

    # -----------------------------------------------------------------
    # functionalities
    # -----------------------------------------------------------------

    def _get_mode(self):
        # figure out if we are dealing with a sequana pipeline
        # or a generic one based solely on focused top tab
        index = self.ui.tabWidget_framework.currentIndex()
        if index == 0:
            return "sequana"
        elif index == 1:
            return "generic"
    mode = property(_get_mode)

    def info(self, text):
        colorlog.info(text)
    def error(self, text):
        colorlog.error(text)
    def debug(self, text):
        colorlog.debug(text)
    def critical(self, text):
        colorlog.critical(text)
    def warning(self, text):
        colorlog.warning(text)

    # ----------------------------------------------------------------------
    # Snakemake related (config, running)
    # ----------------------------------------------------------------------

    def fill_until_starting(self, whatever):
        active_list = [w.get_name() for w in self.rule_list if w.get_do_rule()]

        self.ui.until_box.clear()
        self.ui.until_box.addItems([None] + active_list)

        self.ui.starting_box.clear()
        self.ui.starting_box.addItems([None] + active_list)

    # ----------------------------------------------------------
    #  Config file related
    # ---------------------------------------------------------

    def _set_focus_on_config_tab(self):
        # Set focus on config file
        if self._ipython_tab:
            self.ui.tabs.setCurrentIndex(3)
        else:
            self.ui.tabs.setCurrentIndex(2)

    def _get_configfile(self):
        if self.mode == "sequana" and self._configfile is None:
            msg = CriticalMessage("Please select the input and working directory")
            msg.exec_()
            return
        else:
            return self._configfile
    configfile = property(_get_configfile)

    def _read_sequana_config(self, config_file):
        try:
            configfile = snaketools.SequanaConfig(config_file)
            configfile.cleanup() # set all empty strings and %()s to None
            self._configfile = configfile
        except AssertionError:
            self._configfile = None
            self.clear_layout(self.form)
            print("Warning: could not parse the config file")
        self.create_base_form()
        self._set_focus_on_config_tab()

    def _read_generic_config(self):
        print("selected config file")
        filename = self._generic_config.get_filenames()
        if len(filename):
            try:
                configfile = snaketools.SequanaConfig(filename, mode="generic")
            except AssertionError:
                print("Warning: could not parse the config file")
                return
            self._configfile = configfile
            self.create_base_form()
            self._set_focus_on_config_tab()
        else:
            self._configfile = None
            self.clear_layout(self.form)

    # --------------------------------------------------------------------
    # Others
    # --------------------------------------------------------------------

    def clear_layout(self, layout):
        """ Clean all widgets contained in a layout. """
        while layout.count():
            child = layout.takeAt(0)
            if child.widget() is not None:
                child.widget().deleteLater()
            elif child.layout() is not None:
                self.clear_layout(child.layout())

    # --------------------------------------------------------------------
    # Running snakemake
    # --------------------------------------------------------------------

    def snakemake_data_stdout(self):
        """ Read standard output of snakemake process """
        data = str(self.process.readAllStandardOutput())
        self.shell += data
        self.update_progress_bar(data)

        for this in data.split("\\n"):
            line = this.strip()
            if line and len(line) > 3 and "complete in" not in line: # prevent all b'' strings
                if line.startswith("b'"):
                    line = line[2:]
                line = line.replace("\\r","")
                line = line.replace("\\t", "&nbsp;"*4)
                self.output.append('<p style="color:blue">' + line +'</p>')
            self.output.moveCursor(END)

    def snakemake_data_error(self):
        """ Read error output of snakemake process """
        error = str(self.process.readAllStandardError())
        self.shell_error += error
        self.update_progress_bar(error)
        for this in error.split("\\n"):
            line = this.strip()
            if line and len(line) > 3 and "complete in" not in line: # prevent all b'' strings
                if line.startswith("b'"):
                    line = line[2:]
                line = line.replace("\\r","")
                line = line.replace("\\t","&nbsp;"*4)
                self.output.append('<p style="color:red">' + line +'</p>')
                self.output.moveCursor(END)

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

    def start_sequana(self):
        # set focus on the snakemake output
        self.ui.tabs.setCurrentIndex(0)
        self.shell_error = ""
        self.shell = ""

        # the progress bar
        pal = self.ui.progressBar.palette()
        pal.setColor(QtGui.QPalette.Highlight, self._colors['blue'])
        self.ui.progressBar.setPalette(pal)
        self.ui.progressBar.setValue(1)

        # Set the regex to catch steps
        self._step_regex = re.compile("([0-9]+) of ([0-9]+) steps")

        # Prepare the command and working directory.
        working_dir = self.working_dir.get_filenames()
        snakefile = working_dir + os.sep + "Snakefile"
        assert os.path.exists(snakefile)

        snakemake_args = self._get_snakemake_command(snakefile)

        self.info("Starting process with %s " % snakemake_args)

        pref = self.preferences_dialog.ui
        process = pref.preferences_options_general_process_value.currentText()

        if process == "qt":
            self.process.setWorkingDirectory(working_dir)
            self.process.start("snakemake", snakemake_args)
        else:
            self.cmd = ['snakemake'] + snakemake_args
            self.cwd = working_dir
            snakemake_proc = sp.Popen(self.cmd,
                cwd=working_dir)

    # -------------------------------------------------------------------
    # Create the base form
    # -------------------------------------------------------------------

    def create_base_form(self):
        """ Create form with all options necessary for a pipeline.
        """
        self.info("Interpreting config file")
        self.clear_layout(self.form)
        rules_list = list(self._configfile._yaml_code.keys())
        rules_list.sort()
        self.necessary_dict = {}
        self.rule_list = []

        for count, rule in enumerate(rules_list):
            # Check if this is a dictionnary
            contains = self._configfile._yaml_code[rule]
            if isinstance(contains, dict) and (
                    rule not in SequanaGUI._not_a_rule):
                rule_box = Ruleform(rule, contains, count, self._browser_keyword)
                self.form.addWidget(rule_box)
                self.rule_list.append(rule_box)
                rule_box.connect_do(self.fill_until_starting)
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

        if self.mode == "generic" and len(self.necessary_dict):
            rule_box = Ruleform(self._undefined_section, self.necessary_dict,
                                -1, generic=True)
            self.form.addWidget(rule_box)

    # ----------------------------------------------------------
    # STOP footer button
    # ----------------------------------------------------------

    def click_stop(self):
        """The stop button"""
        print("stopped snakemake manually")
        pal = self.ui.progressBar.palette()
        pal.setColor(QtGui.QPalette.Highlight, self._colors['red'])
        self.ui.progressBar.setPalette(pal)
        self.process.finished.disconnect() # we do not want to call end_run()
        self.ui.run_btn.setEnabled(True)
        self.ui.stop_btn.setEnabled(False)
        self.process.close()

    # --------------------------------------------------------------------
    # Progress bar
    # --------------------------------------------------------------------

    def update_progress_bar(self, line):
        """ Parse with a regex to retrieve current step and total step.
        """
        grouprex = self._step_regex.findall(line)
        if grouprex:
            step = int(grouprex[0][0]) / float(grouprex[0][1]) * 100
            self.ui.progressBar.setValue(step)
        if "Nothing to be done" in line:
            self.ui.progressBar.setValue(100)

    def start_progress(self):
        self.ui.progressBar.setRange(0,1)

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

    def save_config_file(self):
        self.info("Saving config file")

        try:
            form_dict = dict(self.create_form_dict(self.form),
                                 **self.necessary_dict)
        except AttributeError:
            msg = WarningMessage("You must choose a pipeline before saving.")
            msg.exec_()
            return

        if self._undefined_section in form_dict.keys():
            del form_dict[self._undefined_section]

        # get samples names or input_directory
        if self.mode == "sequana":
            if self.ui.tabs_browser.currentIndex() == 1:
                form_dict["samples"] = (
                    self.ui.tabs_browser.currentWidget().get_filenames())
            else:
                form_dict["input_directory"] = (
                    self.ui.tabs_browser.currentWidget().get_filenames())

        # Let us update the attribute with the content of the form
        self.configfile._yaml_code.update(form_dict)

        if self.working_dir.path_is_setup():
            yaml_path = self.working_dir.get_filenames() + "/config.yaml"
            if self.mode == "sequana":
                self.warning("copy requirements (if any)")
                self.configfile.copy_requirements(target=self.working_dir.get_filenames())

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
                    self.configfile.save(yaml_path, cleanup=False)
                    self.ui.dag_btn.setEnabled(True)
            else:
                self.configfile.save(yaml_path, cleanup=False)
                self.ui.dag_btn.setEnabled(True)
        else:
            msg = WarningMessage("You must set a working directory", self)
            msg.exec_()

    def _copy_snakefile(self):
        if self.working_dir.path_is_setup():
            if self.snakefile:
                working_dir = self.working_dir.get_filenames()
                working_dir += os.sep + "Snakefile"
                self.info("Copying snakefile in %s " % working_dir)
                try:
                    shutil.copy(self.snakefile, working_dir)
                except:
                    self.warning("cannot overwrite existing (identical) file")
                    pass

    def switch_run(self):
        # This functions copies the snakefile in the working directory
        # if possible.
        self._copy_snakefile()

        # Run is on if working dir is on and
        # 1. snakefile is present and config file if sequana mode
        # 2. snakefile is present irrespective on config file in generic mode
        if self.working_dir.path_is_setup():
            if self.mode == "sequana":
                # requires the directory or samples to be set
                if self.ui.tabs_browser.currentWidget().path_is_setup():
                    if self.pipeline_is_chosen:
                        return self.ui.run_btn.setEnabled(True)
            else: # generic
                if self.snakefile:
                    return self.ui.run_btn.setEnabled(True)
        return self.ui.run_btn.setEnabled(False)

    # -----------------------------------------------------------------------
    # UNLOCK footer button
    # -----------------------------------------------------------------------

    def unlock_snakemake(self):

        if self.working_dir.get_filenames() == "":
            return

        working_dir = self.working_dir.get_filenames()

        snakefile = self._get_snakefile()

        if os.path.exists(snakefile) is False:
            print("config not found. should not happen")
        else:
            self.cmd = ['snakemake', "-s", snakefile, "--unlock"]
            self.cwd = working_dir
            print(self.cmd)
            print("Please wait a second")
            print(working_dir)
            snakemake_proc = sp.Popen(self.cmd,
                cwd=working_dir)
        self.ui.run_btn.setEnabled(True)
        self.ui.stop_btn.setEnabled(False)

    def _get_snakefile(self):
        if self.mode == "sequana":
            if self.pipeline_is_chosen:
                snakefile = Module(self.ui.choice_button.currentText()).snakefile
            else:
                msg = CriticalMessage("Please select a pipeline first")
                msg.exec_()
                return
        else:
            snakefile = self._generic_config.get_filenames()
        return snakefile

    # -----------------------------------------------------------------------
    # DAG footer button
    # -----------------------------------------------------------------------

    def show_dag(self):
        snakefile = self._get_snakefile()

        # The config must have been saved, so we just need to copy it
        working_dir = self.working_dir.get_filenames()

        # copy the pipeline in the working directory
        shutil.copy(snakefile, working_dir)
        snakefile = os.path.basename(snakefile)

        svg_filename = self._tempdir.path() + os.sep + "test.svg"

        snakemake_line = ["snakemake", "-s", snakefile]
        snakemake_line += ["--rulegraph"]
        if self.configfile:
            # make sure to copy the config file
            snakemake_line += ["--configfile", "config.yaml"]
        snakemake_line += self.get_until_starting_option()

        self.process1 = QtCore.QProcess(self)
        self.process2 = QtCore.QProcess(self)
        self.process1.setWorkingDirectory(working_dir)
        self.process1.setStandardOutputProcess(self.process2)

        self.process1.start("snakemake", snakemake_line[1:])
        self.process2.start("dot", ["-Tsvg", "-o", svg_filename])

        self.process1.waitForFinished(50000)
        self.process2.waitForFinished(50000)

        if os.path.exists(svg_filename):
            diag = SVGDialog(svg_filename)
            diag.exec_()
        else:
            msg = 'Could not create the DAG file.'
            error = str(self.process1.readAllStandardError())
            msg = CriticalMessage(msg, error)
            msg.exec_()
            return

    def open_report(self, filename="multi_summary.html"):
        dialog = self.preferences_dialog.ui # an alias
        if self.pipeline_is_chosen and self.working_dir.get_filenames():
            filename = self.working_dir.get_filenames() + filename
            if os.path.exists(filename) is False:
                WarningMessage("""multi_summary.html not found.
                Most probably the analysis did not finish correctly""")
                return
            url = "file://" + filename
            browser = dialog.preferences_options_general_browser_value.currentText()

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

        if self.mode == "sequana":
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
        else:
            snakefile = self._generic_snakefile.get_filenames()
            config_file = self._generic_config.get_filenames()
            working_dir = self.working_dir.get_filenames()
            if len(config_file) == 0 or len(snakefile) == 0:
                msg = WarningMessage("You must set the pipeline first")
                msg.exec_()
                self.working_dir.set_empty_path()
                return False

        try:
            # Load the existing config
            if self.mode == "sequana":
                cfg = snaketools.SequanaConfig(config_file)
                cfg.cleanup() # set all empty strings and %()s to None
            else:
                cfg = snaketools.SequanaConfig(config_file, mode="generic")
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

        if set(self.configfile._yaml_code.keys()) == set(config_dict.keys()):
            msg = QW.QMessageBox(
                QW.QMessageBox.Question, "Question",
                "A config file already exist in the working directory.\n"
                "Do you want to import its content ?",
                QW.QMessageBox.Yes | QW.QMessageBox.No,
                self, Qt.Dialog | Qt.CustomizeWindowHint)
            if msg.exec_() == 16384:
                self.configfile._yaml_code.update(config_dict)
                self.create_base_form()
                self.fill_until_starting(self.rule_list)
        return True

    def eventFilter(self, source, event):
        """ Inactivate wheel event of combobox
        """
        if event.type() == QtCore.QEvent.Wheel and source is self.ui.choice_button:
            return True
        return False

    # ---------------------------------------------------
    #  settings and close
    # ---------------------------------------------------

    def read_settings(self):
        settings = QtCore.QSettings("sequana_gui", "mainapp")
        if settings.value("tab_position") is None:
            return
        # tab snakemake output/logger/ipython
        index = settings.value("tab_position")
        self.ui.tabs.setCurrentIndex(int(index))

    def write_settings(self):
        settings = QtCore.QSettings("sequana_gui", "mainapp")

        # tab snakemake output/logger/ipython
        index = self.ui.tabs.currentIndex()
        settings.setValue("tab_position", index)

    def _close(self):
        self.write_settings()
        self._tempdir.remove()
        try:self.browser.close()
        except:pass

    def closeEvent(self, event):
        #Close button (red X)
        self._close()

    def close(self):
        # Menu or ctrl+q
        self._close()
        super().close()


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
    else:
        # Show the splash screen for a few seconds
        splash_pix = QtGui.QPixmap(filename)
        splash = QW.QSplashScreen(splash_pix, Qt.WindowStaysOnTopHint)
        splash.setMask(splash_pix.mask())
        splash.show()

        for i in range(0, 100):
            t = time.time()
            while time.time() < t + 0.5/100.:
                app.processEvents()

        app.processEvents()
        sequana = SequanaGUI()
        sequana.show()
        splash.finish(sequana)

    # Make sure the main window is the active one
    sequana.raise_()
    sequana.activateWindow()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()

