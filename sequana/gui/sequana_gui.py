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
import argparse
from optparse import OptionParser

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
from easydev import SmartFormatter

import colorlog
import logging

import signal


# TEST Sequana
# 1 select pipeline
# 2 select working dir
# 3 cancel working dir
# 4 save config
# 5 select working dir
# 6 save config
# 7 run

# Test generic 
# try with and without the config file


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
        self.bgcolor = "black"
        self.widget.setStyleSheet("background-color: %s" % self.bgcolor)
        

    def emit(self, record):
        formatter = """<span style="color:%(color)s;
                        font-weight:%(weight)s">%(msg)s</span>"""
        # "\e[1;31m This is red text \e[0m"
        self.record = record
        msg = self.format(record)
        msg = msg.rstrip("\x1b[0m")
        if msg.startswith('\x1b[31m\x1b[47m'): # critical
            msg = msg.replace("\x1b[31m\x1b[47m", "")
            params = {'msg':msg, 'weight':"bold", "color":"red"}
            self.widget.appendHtml(formatter % params)
        elif msg.startswith('\x1b[32m'): # info
            msg = msg.replace("\x1b[32m", "")
            params = {'msg':msg, 'weight':"normal", "color":"green"}
            self.widget.appendHtml(formatter % params)
        elif msg.startswith('\x1b[33m'): # warning
            msg = msg.replace("\x1b[33m", "")
            params = {'msg':msg, 'weight':"normal", "color":"yellow"}
            self.widget.appendHtml(formatter % params)
        elif msg.startswith('\x1b[31m'): # error
            msg = msg.replace("\x1b[31m", "")
            params = {'msg':msg, 'weight':"normal", "color":"red"}
            self.widget.appendHtml(formatter % params)
        elif msg.startswith('\x1b[36m'): # debug
            msg = msg.replace("\x1b[36m", "")
            params = {'msg':msg, 'weight':"normal", "color":"cyan"}
            self.widget.appendHtml(formatter % params)
        else:
            self.widget.appendHtml(msg)
        self.msg = msg


class SequanaGUI(QMainWindow):
    """

    If quiet, progress bar cannot work.

    - do not copy again requirements if already there
    - extension of the different widgets ?

    Developer Guide
    ------------------

    - The GUI is designed with qt designer as much as possible.
    - All GUI objects are in the **ui** attributes. Additional dialog such as the
      snakemake and preferences dialog have their own modules and stored in attributes
      ending in _dialog

    """
    _not_a_rule = {"requirements", "gatk_bin", "input_directory",
                    "input_samples", "input_pattern", "ignore"}
    _browser_keyword = {"reference"}

    def __init__(self, parent=None, ipython=True, user_options={}):
        super(SequanaGUI, self).__init__(parent=parent)

        self._tempdir = QTemporaryDir()
        self.shell = ""
        self.shell_error = ""
        self._colors = {
            'green': QtGui.QColor(0,170,0),
            'red': QtGui.QColor(170,0,0),
            'orange': QtGui.QColor(170,150,0),
            'blue': QtGui.QColor(0,90,154),
        }

        # some global attributes
        self.pipeline_is_chosen = False
        self._undefined_section = "Parameters in no sections/rules"

        self._ipython_tab = ipython
        self.initUI()
        self.read_settings()
        self._save_tooltips() # save all tooltips

        self.setStyleSheet("""QToolTip {
                           background-color: #aabbcc;
                           color: black;
                           border-style: double;
                           border-width: 3px;
                           border-color: green;
                           border-radius: 5px;
                           margin:3px;
                           opacity: 220;
                           }
                            """)

        # User option
        if "wkdir" in user_options and user_options.wkdir is not None:
            self.info("Setting working directory")
            if os.path.exists(user_options.wkdir) is False:
                easydev.mkdirs(user_options.wkdir)
            # We must use the absolute path
            abspath = os.path.abspath(user_options.wkdir)
            self.working_dir.set_filenames(abspath)

        if "pipeline" in user_options and user_options.pipeline is not None:
            self.info("Setting Sequana pipeline")
            # TODO: check that the pipeline exsits in the combobox
            index = self.ui.choice_button.findText(user_options.pipeline)
            self.ui.choice_button.setCurrentIndex(index)
            self._set_focus_on_pipeline_tab()

        if "input_directory" in user_options and \
                user_options.input_directory is not None:
            directory = user_options.input_directory
            self.info("Setting Sequana input directory")
            if directory and os.path.exists(directory) is False:
                self.warning("%s does not exist" % directory)
            elif directory:
                abspath = os.path.abspath(user_options.input_directory)
                self._sequana_directory_tab.set_filenames(abspath)
        # Can we run the pipeline ?
        self.switch_run()

    def initUI(self):

        # Set up the user interface from Designer. This is the general layout
        # without dedicated widgets and connections
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        # 2 more dialogs from designer
        self.preferences_dialog = PreferencesDialog(self)
        self.snakemake_dialog = SnakemakeDialog(self)


        self.preferences_dialog.ui.buttonBox.accepted.connect(self.set_level)

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
        self.working_dir.clicked_connect(self._copy_snakefile)
        self.working_dir.clicked_connect(self._load_config)

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
        self.process.started.connect(lambda: self.ui.unlock_btn.setEnabled(False))
        self.process.started.connect(lambda: self.start_progress)

        self.process.finished.connect(lambda: self.ui.run_btn.setEnabled(True))
        self.process.finished.connect(lambda: self.ui.stop_btn.setEnabled(False))
        self.process.finished.connect(lambda: self.ui.unlock_btn.setEnabled(True))
        self.process.finished.connect(self.end_run)

        self.process.readyReadStandardOutput.connect(
            self.snakemake_data_stdout)
        self.process.readyReadStandardError.connect(self.snakemake_data_error)

    def _close_last_widget(self):
        try:
            self._last_widget.close()
        except:
            print("nothing to close")
    #|-----------------------------------------------------|
    #|                       MENU related                  |
    #|-----------------------------------------------------|
    def menuAbout(self):
        from sequana import version
        url = 'sequana.readthedocs.io'
        widget = About()
        widget.setIcon(QW.QMessageBox.Information)
        widget.setText("Sequana version %s " % version)
        widget.setInformativeText("""
            Online documentation on <a href="http://%(url)s">%(url)s</a>
            <br>
            <br>
            Authors: Thomas Cokelaer and Dimitri Desvillechabrol, 2016
            """ % {"url": url})
        widget.setWindowTitle("Sequana")
        widget.setStandardButtons(QW.QMessageBox.Ok)
        self._last_widget = widget
        retval = widget.exec_()

    def menuHelp(self):
        url = 'sequana.readthedocs.io'
        msg = About()
        msg.setIcon(QW.QMessageBox.Information)
        msg.setText("<h1>Sequana GUI help</h1>")

        pipelines_text = "<ul>\n"
        url = "http://sequana.readthedocs.io/en/master"
        for pipeline in snaketools.pipeline_names:
            pipelines_text += '    <li><a href="%(url)s/pipeline_%(name)s.html">%(name)s</a></li>\n' %\
              {"url":url,"name":pipeline}
        pipelines_text += "</ul>"

        msg.setInformativeText("""<div>
This GUI can be used to run either Sequana pipelines (see 
<a href="http://sequana.readthedocs.io">Sequana.readthedocs.io</a> for details) or Snakefiles 
(see <a href="http://snakemake.readthedocs.io">snakemake.readthedocs.io</a>for details)

        In both cases, a working directory must be set where the Snakefile 
        and possibly a configuration file will be copied.

        The generic Snakefile must be executable. 

        <h2>Sequana pipelines</h2>
        There are download automatically with their config file from the Sequana
        library. Here is a typical set of actions to run Sequana pipelines:

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
        <li> Open the report if successful</li></ol>

        <h2> Generic pipelines </h2>
        Same as above except that you have to select the config file (if any).


        <h2> Sequana pipeline dedicated help </help>
             %(pipelines)s
        </div>
        """ % {"url": url, "pipelines": pipelines_text})
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

    def set_level(self):
        # Set the level of the logging system
        pref = self.preferences_dialog.ui
        level = pref.preferences_options_general_logging_value.currentText()
        level =  getattr(logging, level)
        colorlog.getLogger().setLevel(level)

    # ---------------------------------------------------------------
    # More GUI / reading the snakefile (sequana or generic)
    # ---------------------------------------------------------------
    def set_sequana_pipeline(self):
        # Fill the pipelines
        snaketools.pipeline_names.sort()
        self.ui.choice_button.addItems(snaketools.pipeline_names)
        self.ui.choice_button.currentIndexChanged[str].connect(
            self.on_sequana_pipeline_choice)
        self.ui.choice_button.activated.connect(self._copy_snakefile)
        self.ui.choice_button.activated.connect(self._load_config)
        self.ui.choice_button.installEventFilter(self)

        # Set the file browser for paired files
        fastq_filter = "Fastq file (*.fastq *.fastq.gz *.fq *.fq.gz)"
        self._sequana_paired_tab = FileBrowser(paired=True, file_filter=fastq_filter)
        self._sequana_paired_tab.clicked_connect(self.switch_run)

        # Set the file browser input_directory tab
        self._sequana_directory_tab = FileBrowser(directory=True)
        self._sequana_directory_tab.clicked_connect(self.switch_run)

        # fill the tabs_browser
        self.ui.tabs_browser.addTab(self._sequana_directory_tab, "Directory")
        self.ui.tabs_browser.addTab(self._sequana_paired_tab, "Sample")
        self.ui.tabs_browser.removeTab(0)
        self.ui.tabs_browser.removeTab(0)

    @QtCore.pyqtSlot(str)
    def on_sequana_pipeline_choice(self, index):
        """ Change options form when user change the pipeline."""
        if self.ui.choice_button.findText(index) == 0:
            self.pipeline_is_chosen = False
            self.clear_form()
            return

        self.info("Reading sequana %s pipeline" % index)
        module = snaketools.Module(index)
        config_file = module._get_config()
        self.snakefile = module.snakefile
        self._read_sequana_config(config_file)
        self.fill_until_starting(self.rule_list)
        self.pipeline_is_chosen = True
        self.info("--------pipeline chosen")

    def set_generic_pipeline(self):
        # Set the file browser
        self._generic_snakefile = FileBrowser(directory=False)
        self._generic_config = FileBrowser(directory=False,
            file_filter="YAML file (*.json *.yaml)")
        self._generic_snakefile.clicked_connect(self._read_generic_snakefile)
        # no switch run for the config file, which is optional
        self._generic_config.clicked_connect(self._read_generic_config)
        self._generic_config.clicked_connect(self._copy_snakefile)

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

    def _set_focus_on_pipeline_tab(self):
        self.ui.tabs_pipeline.setCurrentIndex(0)

    # ---------------------------------------------------------------------
    # Fotter connectors
    # ---------------------------------------------------------------------

    def connect_footer_buttons(self):
        self.ui.run_btn.setEnabled(False)
        self.ui.run_btn.clicked.connect(self.start_sequana)

        self.ui.stop_btn.clicked.connect(self.click_stop)
        self.ui.stop_btn.setEnabled(False)

        self.ui.unlock_btn.clicked.connect(self.ui.run_btn.setEnabled)
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
        index = self.ui.tabs_pipeline.currentIndex()
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
            self.clear_form()
            self.warning("could not parse the config file")
        self.create_base_form()
        self._set_focus_on_config_tab()

    def _read_generic_config(self):
        self.info("selected config file")
        filename = self._generic_config.get_filenames()
        if len(filename):
            try:
                configfile = snaketools.SequanaConfig(filename)
            except AssertionError:
                self.warning("Warning: could not parse the config file")
                return
            self._configfile = configfile
            self.create_base_form()
            self._set_focus_on_config_tab()
        else:
            self._configfile = None
            self.clear_form()

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
                line = line.replace("b'\\r'", "")
                line = line.replace("b'\r'", "")
                line = line.replace("b'\\r '", "")
                line = line.replace("b'\r '", "")
                line = line.replace("b' '", "")
                if len(line.strip()) == 0:
                    continue
                line = line.replace("\\t", "&nbsp;"*4)
                self.output.append('<font style="color:blue">' + line +'</font>')

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
                    line.rstrip("'")
                line = line.replace("\\r","")
                line = line.replace("\\t","&nbsp;"*4)
                grouprex = self._step_regex.findall(line)
                if grouprex:
                    self.output.append('<font style="color:orange">' + line +'</font>')
                elif "Error" in line:
                    self.output.append('<font style="color:red">' + line +'</font>')
                else:
                    self.output.append('<font style="color:green">' + line +'</font>')

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
        snakemake_line = ["-s", snakefile, "--stat", "stats.txt", "-p"]

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

        self.info("Starting process with %s " % " ".join(snakemake_args))

        pref = self.preferences_dialog.ui

        self.process.setWorkingDirectory(working_dir)
        self.process.start("snakemake", snakemake_args)

    # -------------------------------------------------------------------
    # Create the base form
    # -------------------------------------------------------------------

    def create_base_form(self):
        """ Create form with all options necessary for a pipeline.
        """
        self.info("Interpreting config file")
        self.clear_form()
        rules_list = list(self._configfile._yaml_code.keys())
        rules_list.sort()
        self.necessary_dict = {}
        self.rule_list = []

        # A section looks like a large comments :
        #   #========
        #   # valid python docstring to be interepreted by sphinx
        #   # 
        #   section:   # short comment
        #      item: value # comment interne

        # This
        # cfg._yaml_code.ca.items['adapter_removal']
        # is a list of 4 items
        # [None, largecomment, shortcomment, None]

        # For each section, we create a widget (RuleForm). For isntance, first,
        # one is accessible asfollows:
        # gui.form.itemAt(0).widget()

        for count, rule in enumerate(rules_list):
            # Check if this is a dictionnary
            contains = self._configfile._yaml_code[rule]
            comments = self._configfile.get_section_long_comment(rule)
            if isinstance(contains, dict) and (
                    rule not in SequanaGUI._not_a_rule):
                rule_box = Ruleform(rule, contains, count, self._browser_keyword)
                self.comments= comments
                if comments:
                    # Try to interpret it with sphinx
                    from pyquickhelper.helpgen import docstring2html
                    docstring = []
                    for line in comments.split("\n"):
                        if "#############" in line:
                            pass
                        else:
                            if len(line)<2: # an empty line (to keep)
                                docstring.append("")
                            else:
                                docstring.append(line[2:]) # strip the "# " characters

                    docstring = "\n".join(docstring)
                    try:
                        comments = docstring2html(docstring).data
                        comments.replace('class="document"', 'style=""')
                    except:
                        self.warning("Could not interpret docstring of %s" % rule)


                    rule_box.setToolTip(comments)
                else:
                    rule_box.setToolTip("")
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
        self.warning("stopped snakemake process manually. You may need to use unlock")
        pal = self.ui.progressBar.palette()
        pal.setColor(QtGui.QPalette.Highlight, self._colors['orange'])
        self.ui.progressBar.setPalette(pal)

        if self.process.state() != 0:
            self.info("Process running, stopping it... ")
            """try:
                self.process.finished.disconnect() # we do not want to call end_run()
            except Exception as err:
                self.critical(err)
                self.error("process disconnect failed")
            """
            #self.process.close()
            import os
            os.kill(self.process.pid(), signal.SIGINT) 
            self.info("ok ")
        self.ui.run_btn.setEnabled(True)
        self.ui.stop_btn.setEnabled(False)

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
            self.info('Run done. Status: successful')
        else:
            pal.setColor(QtGui.QPalette.Highlight, self._colors['red'])
            self.ui.progressBar.setPalette(pal)
            text = 'Run manually to check the exact error or check the log.'
            if "--unlock" in self.shell_error:
                text += "<br>You may need to unlock the directory. "
                text += "click on Unlock button"
                self.critical(text)
            return

    def save_config_file(self, force=False):
        self.info('DEBUG: SAVE_CONFIG_FILE')
        self.info("Saving config file")

        try:
            form_dict = dict(self.create_form_dict(self.form),
                                 **self.necessary_dict)
        except AttributeError as err:
            logger.error(err)
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
        self.configfile.config.update(form_dict)
        self.configfile._update_yaml()
        # We must update the config and then the yaml to keep the comments. This
        # lose the comments:
        # self.configfile._yaml_code.update(form_dict)

        if self.working_dir.path_is_setup():
            yaml_path = self.working_dir.get_filenames() + "/config.yaml"
            #if self.configfile.mode == "sequana":
            #    self.warning("copy requirements (if any)")
            #    self.configfile.copy_requirements(target=self.working_dir.get_filenames())

            if os.path.isfile(yaml_path):
                if force == False:
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
                self.configfile.save(yaml_path, cleanup=False)
                self.ui.dag_btn.setEnabled(True)
        else:
            msg = WarningMessage("You must set a working directory", self)
            msg.exec_()

    def _copy_snakefile(self):
        # When working dir changes, we try to copy the snakefile
        # if set.
        if self.snakefile is None:
            return # nothing to be done

        if self.working_dir.path_is_setup() is False:
            return

        self.info('DEBUG: COPY SNAKEFILE')

        working_dir = self.working_dir.get_filenames() + os.sep + "Snakefile"
        self.info("Copying snakefile in %s " % working_dir)
        try:
            shutil.copy(self.snakefile, working_dir)
        except:
            self.warning("cannot overwrite existing file. (Probably identical)")
            return

        self.switch_run()


    def switch_run(self):
        # This functions copies the snakefile in the working directory
        # if possible.

        # Run is on if working dir is on AND
        # 1. for sequana: pipeline is set
        # 2. for generic: snakefile is present irrespective of config file in generic mode
        if self.working_dir.path_is_setup():
            self.info('DEBUG: SWITCH RUN')
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

        working_dir = self.working_dir.get_filenames()
        if working_dir == "":
            return
        snakefile = self._get_snakefile()

        # FIXME this does not work as expected
        self.ui.run_btn.setEnabled(False)
        self.update()
        time.sleep(2)

        if os.path.exists(snakefile) is False:
            self.warning("config not found. should not happen")
        else:
            self.cmd = ['snakemake', "-s", snakefile, "--unlock"]
            self.info("Running " + " ".join(self.cmd))
            self.info("Please wait a second. Unlocking working directory")
            snakemake_proc = sp.Popen(self.cmd, cwd=working_dir)
            snakemake_proc.wait()
        self.info("unlocking done")
        self.output.append('<font style="color:brown">Unlocking working directory</font>')

        self.ui.run_btn.setEnabled(True)
        self.ui.stop_btn.setEnabled(False)

    def _get_snakefile(self):
        self.info('DEBUG: Entering _get_snakefile')
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
        self.info("Creating DAG image.")
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

    def open_report(self):

        pref = self.preferences_dialog.ui
        filename = pref.preferences_options_general_htmlpage_value.text()

        if self.working_dir.get_filenames() == "":
            self.error("Working directory not set yet")
            return

        filename = self.working_dir.get_filenames() + os.sep + filename
        if os.path.exists(filename) is False:
            self.error("%s page does not exist. Check the preferences dialog." % filename)
            return
        else:
            self.info("Reading and openning %s" % filename)

        dialog = self.preferences_dialog.ui # an alias
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

    def clear_form(self):
        self.clear_layout(self.form)

    def _read_config_inside_working_dir(self, cfg):
        self.info('DEBUG: Entering _read_config %s' % cfg)
        if self.mode == "sequana":
            # this should always work but who knows
            try:
                cfg = snaketools.SequanaConfig(cfg)
                cfg.cleanup() # set all empty strings and %()s to None
                return cfg.config
            except:
                self.critical("Could not interpret the sequana config file")
        else:
            if len(cfg) == 0:
                return None
            # otherwise, we read it
            try:
                cfg = snaketools.SequanaConfig(cfg)
                return cfg.config
            except AssertionError:
                msg = WarningMessage("Could not parse the config file.")
                msg.exec_()
            except Exception as err:
                self.critical(err)
                msg = WarningMessage(
                    "Unexpected error while checking the config file %s."
                    % cfg + str(err))
                msg.exec_()

    def _load_config(self):
        # Once the working dir is clicked, we may want to update
        # the config file.
        self.info('DEBUG: Entering _load_config')

        if self.mode == "sequana":
            # If working directory is not set yet, nothing to load
            if self.working_dir.get_filenames() == "":
                return
            else:
                config_file = self.working_dir.get_filenames() + "/config.yaml"
                if not os.path.isfile(config_file):
                    self.info("config file %s not found. " % config_file +
                        "Will use the pipeline original config file.")
                    return

            # If pipeline not set, we cannot read the config file
            if not self.pipeline_is_chosen:
                msg = WarningMessage("A config.yaml file already exist in this "
                                 "directory. Please, choose a pipeline to "
                                 "know if the existing config file correspond "
                                 "to your pipeline.", self)
                msg.exec_()
                return
        else: # a generic file
            config_file = self._generic_config.get_filenames()
            if config_file is None:
                self.info("generic case: no config file found")
                self.clear_form()
                return

        # Now, config_file should be a valid file since it was selected by the
        # file browser
        config_dict = self._read_config_inside_working_dir(config_file)

        if config_dict is None:
            # reset the layout
            self.warning("config_dict is None")
            self.clear_form()
            return

        print(config_dict)
        if set(self.configfile._yaml_code.keys()) == set(config_dict.keys()):
            msg = QW.QMessageBox(
                QW.QMessageBox.Question, "Question",
                "A config file already exist in the working directory.\n"
                "Do you want to import its content ?",
                QW.QMessageBox.Yes | QW.QMessageBox.No,
                self, Qt.Dialog | Qt.CustomizeWindowHint)
            # Yes == 16384
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

    #-------------------------------------------------------
    # tooltips
    # ------------------------------------------------------
    def _save_tooltips(self):
        self._tooltips = {}
        self._tooltips["working_directory"] = self.ui.working_directory.toolTip()

    """
    def hide_tooltips(self):
        self.ui.working_directory.setToolTip("")

    def show_tooltips(self):
        self.ui.working_directory.setToolTip(self._tooltips["working_directory"])
    """

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
        # end any process running that may be running
        self.click_stop()

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


class Options(argparse.ArgumentParser):
    def __init__(self, prog="sequana_gui"):
        usage = """dkfj skfk"""
        description = """"""
        super(Options, self).__init__(usage=usage, prog=prog,
            description=description, formatter_class=SmartFormatter)
        group = self.add_argument_group("GENERAL")
        group.add_argument("-w", "--working-directory", dest="wkdir",
            help="Set working directory")
        group = self.add_argument_group("SEQUANA")
        group.add_argument("-p", "--pipeline", dest="pipeline",
            help="A valid sequana pipeline name")
        group.add_argument("-s", "--no-splash", dest="nosplash",
            action="store_true",
            help="No splash screen")
        group.add_argument("-i", "--input-directory", dest="input_directory",
            default=None,
            help="input directory where to find the input data")


def main(args=None):

    if args is None:
        args = sys.argv[:]
    user_options = Options()
    options = user_options.parse_args(args[1:])

    signal.signal(signal.SIGINT, sigint_handler)
    app = QW.QApplication(sys.argv)

    filename = sequana_data("drawing.png", "../gui")

    if options.nosplash:
        app.processEvents()
        sequana = SequanaGUI(user_options=options)
        sequana.show()
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
        sequana = SequanaGUI(user_options=options)
        sequana.show()
        splash.finish(sequana)

    # Make sure the main window is the active one
    sequana.raise_()
    sequana.activateWindow()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()

