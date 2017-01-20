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
import argparse
import signal
from optparse import OptionParser

from sequana.gui.ui_mainwindow import Ui_MainWindow
from sequana.gui.browser import Browser
from sequana.gui.ipython import QIPythonWidget
from sequana.gui.about import About
from sequana.gui.file_browser import FileBrowser
from sequana.gui.widgets import Ruleform, SVGDialog
from sequana.gui.messages import WarningMessage, CriticalMessage
from sequana.gui.preferences import PreferencesDialog
from sequana.gui.snakemake import SnakemakeDialog
from sequana.gui.tools import Tools, QPlainTextEditLogger

from PyQt5 import QtCore, QtGui
from PyQt5 import QtWidgets as QW
from PyQt5.Qt import QTemporaryDir, QMainWindow
from PyQt5.QtCore import Qt

from sequana import snaketools, sequana_data
from sequana.snaketools import Module

import easydev
import colorlog



def sigint_handler(*args):
    """Handler for the SIGINT signal."""
    sys.stderr.write('\r')
    if QW.QMessageBox.question(None, '', "Are you sure you want to quit?",
                            QW.QMessageBox.Yes | QW.QMessageBox.No,
                            QW.QMessageBox.No) == QW.QMessageBox.Yes:
        QW.QApplication.quit()


class BaseFactory(Tools):
    """Tab on top are all based on this abstract class

    It should provide access to a snakefile and its config file as well
    as working directory.

    Currently, the :class:`SequanaFactory` and :class:`GenericFactory` are
    implemented.


    """
    def __init__(self, mode, run_button):
        self.mode = mode
        self._run_button = run_button

        # And finally the working directory
        self._directory_browser = FileBrowser(directory=True)
        self._directory_browser.clicked_connect(self._switch_off_run)

    def _switch_off_run(self):
        self.info("Switching off run button")
        self._run_button.setEnabled(False)

    def copy(self, source, target):
        print("target = " + target)
        print("target = " + os.path.basename(target))
        if os.path.exists(target):
            save_msg = WarningMessage(
                    "The file {0} already exists in the working directory".format(source))
            save_msg.setInformativeText(
                    "Do you want to overwrite it?")
            save_msg.setStandardButtons(
                    QW.QMessageBox.Yes | QW.QMessageBox.Discard |
                    QW.QMessageBox.Cancel)
            save_msg.setDefaultButton(QW.QMessageBox.Yes)
            # Yes == 16384
            # Save == 2048
            retval = save_msg.exec()
            if retval in [16384, 2048]:
                self.warning("Overwritting %s" % target)
                super(BaseFactory, self).copy(source, target)
        else:
                super(BaseFactory, self).copy(source, target)

    def _copy_snakefile(self, force=False):
        if self.snakefile is None:
            self.info("No pipeline selected yet")
            return # nothing to be done

        if self.directory is None:
            self.info("No working directory selected yet")
            return

        target = self.directory + os.sep + os.path.basename(self.snakefile)
        if os.path.basename(self.snakefile) == target:
            self.warning("%s exists already in %s" % (self.snakefile,
                self.directory))
        else:
            self.info("Copying snakefile in %s " % self.directory)
            self.copy(self.snakefile, target)

    def _copy_configfile(self):
        if self.configfile is None:
            self.info("No config selected yet")
            return # nothing to be done

        if self._directory_browser.path_is_setup() is False:
            self.info("No working directory selected yet")
            return

        self.info("Copying config in %s " % self.directory)
        self.copy(self.configfile, self.directory)

    def _get_directory(self):
        filename = self._directory_browser.get_filenames()
        if len(filename):
            return filename
        else:
            return None
    directory = property(_get_directory)

    def __repr__(self):
        return "%s Factory" % self.mode


class SequanaFactory(BaseFactory):
    def __init__(self, run_button, combobox):
        super(SequanaFactory, self).__init__("sequana", run_button)
        self._config = None
        self._choice_button = combobox

        # Some widgets to be used: a file browser for paired files
        fastq_filter = "Fastq file (*.fastq *.fastq.gz *.fq *.fq.gz)"
        self._sequana_paired_tab = FileBrowser(paired=True, file_filter=fastq_filter)

        # Set the file browser input_directory tab
        self._sequana_directory_tab = FileBrowser(directory=True)

        # triggers/connectors
        self._sequana_directory_tab.clicked_connect(self._switch_off_run)
        self._choice_button.activated.connect(self._switch_off_run)
        self._sequana_paired_tab.clicked_connect(self._switch_off_run)

        # TODO: if change tab A-2 from valid input directory to non-set input
        # samples, call switch off


    def _get_pipeline(self):
        index = self._choice_button.currentIndex()
        if index == 0:
            return None
        else:
            return self._choice_button.currentText()
    pipeline = property(_get_pipeline)

    def _get_snakefile(self):
        if self.pipeline:
            module = snaketools.Module(self.pipeline)
            return module.snakefile
    snakefile = property(_get_snakefile)

    def _get_configfile(self):
        if self.pipeline:
            module = snaketools.Module(self.pipeline)
            return module.config
    configfile = property(_get_configfile)

    def _get_config(self):
        if self.configfile:
            try:
                cfg = snaketools.SequanaConfig(self.configfile)
                return cfg
            except AssertionError:
                self.warning("Warning: could not parse the config file")
                return
    config = property(_get_config)

    def __repr__(self):
        in1 = self._sequana_directory_tab.get_filenames()
        in2 = self._sequana_paired_tab.get_filenames()
        txt = super(SequanaFactory, self).__repr__()
        txt += "\npipeline:%s\ninput:\n - %s\n - %s\ndirectory:%s\n"
        return txt % (self.pipeline, in1, in2, self.directory)


class GenericFactory(BaseFactory):

    def __init__(self, run_button):
        super(GenericFactory, self).__init__("generic", run_button)

        # Define the Snakefile browser and config widgets
        self._snakefile_browser = FileBrowser(directory=False)
        self._config_browser = FileBrowser(directory=False,
            file_filter="YAML file (*.json *.yaml)")

        # when a snakefile or config is chosen, switch off run button
        self._config_browser.clicked_connect(self._switch_off_run)
        self._snakefile_browser.clicked_connect(self._switch_off_run)

    def get_case(self):
        filename = self._snakefile_browser.get_filenames()
        if filename == "":
            return "unset"

        # else
        selection = [x for x in open(filename, "r").readlines()
            if x.startswith('configfile:')]
        if len(selection) > 0:
            return "configfile"
        else:
            # Here there could be a config variable, in which case a config is
            # required externally
            return "none"

    def _return_none(self, this):
        if this is None or len(this) == 0:
            return None
        else:
            return this

    def _get_snakefile(self):
        return self._return_none(self._snakefile_browser.get_filenames())
    snakefile = property(_get_snakefile)

    def _get_configfile(self):
        return self._return_none(self._config_browser.get_filenames())
    configfile = property(_get_configfile)

    def _get_config(self):
        filename = self._return_none(self._config_browser.get_filenames())
        if filename:
            try:
                configfile = snaketools.SequanaConfig(filename)
            except AssertionError:
                self.critical("Could not parse the config file %s" % filename)
                return
            except:
                self.critical("Could not parse the config file %s" % filename)
            return configfile
    config = property(_get_config)

    def _configfile_arg(self):
        # if snakefile  uses configfile:
        if self.config and self.get_case() == "configfile":
            # The only reason for the user to provide a configfile
            # is that the variable config is used and so the --configfile
            # argument is probably required
            return True
        else:
            return False

    def is_runnable(self):
        flag1 = self._directory_browser.path_is_setup()
        flag2 = self._snakefile_browser.path_is_setup()
        flag3 = self._config_browser.path_is_setup()

        # flag1 and flag2 are compulsary
        if flag1 and flag2:
            if self.get_case() == "configfile": # a config must be set
                if flag3:
                    return True
            else:
                if flag3:
                    return True
                else:
                    # possibly no config file is required so we also return
                    # True
                    return True
        return False

    def __repr__(self):
        txt = super(GenericFactory, self).__repr__()
        txt += "\nsnakefile:%s\nconfigfile:%s\ndirectory:%s\n"
        return txt % (self.snakefile, self.configfile, self.directory)


class SequanaGUI(QMainWindow, Tools):
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

        colorlog.getLogger().setLevel("INFO")
        colorlog.info("Welcome to Sequana GUI (aka Sequanix)")

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
        self._undefined_section = "Parameters in no sections/rules"
        self._config = None

        self._ipython_tab = ipython
        self.initUI()
        self.read_settings()

        self.setStyleSheet("""QToolTip {
                           background-color: #aabbcc;
                           color: black;
                           border-style: double;
                           border-width: 3px;
                           border-color: green;
                           border-radius: 5px;
                           margin:3px;
                           opacity: 240;
                           }
                            """)

        # User option.
        def isset(options, key):
            if key in options and getattr(options,key):
                return True
            else:
                return False

        if isset(user_options, "wkdir"):
            self.info("Setting working directory using user's argument %s" %
                user_options.wkdir)
            if os.path.exists(user_options.wkdir) is False:
                easydev.mkdirs(user_options.wkdir)
            # We must use the absolute path
            abspath = os.path.abspath(user_options.wkdir)
            self.sequana_factory._directory_browser.set_filenames(abspath)
            self.generic_factory._directory_browser.set_filenames(abspath)

        if isset(user_options, "snakefile"):
            filename = user_options.snakefile
            if os.path.exists(filename) is True:
                self.info("Setting snakefile using user's argument")
                self.generic_factory._snakefile_browser.set_filenames(filename)
            else:
                self.error("%s does not exist" % filename)
            self.ui.tabs_pipeline.setCurrentIndex(1)

        if isset(user_options, "configfile"):
            filename = user_options.configfile
            if os.path.exists(filename) is True:
                self.info("Setting config file using user's argument")
                self.generic_factory._config_browser.set_filenames(filename)
            self.ui.tabs_pipeline.setCurrentIndex(1)

        if isset(user_options, "pipeline"):
            self.info("Setting Sequana pipeline %s ")
            pipelines = self.sequana_factory.valid_pipelines
            if user_options.pipeline in pipelines:
                index = self.ui.choice_button.findText(user_options.pipeline)
                self.ui.choice_button.setCurrentIndex(index)
                # set focus on pipeline tab
                self.ui.tabs_pipeline.setCurrentIndex(0)
            else:
                self.error("unknown pipeline. Use one of %s " % pipelines)

        if isset(user_options, "input_directory"):
            directory = user_options.input_directory
            self.info("Setting Sequana input directory")
            if directory and os.path.exists(directory) is False:
                self.warning("%s does not exist" % directory)
            elif directory:
                abspath = os.path.abspath(user_options.input_directory)
                self.sequana_factory._sequana_directory_tab.set_filenames(abspath)
            self.ui.tabs_pipeline.setCurrentIndex(0)
            self.ui.tabWidget.setCurrentIndex(0)

        if isset(user_options, "input_files"):
            directory = user_options.input_files
            self.info("Setting Sequana input files")
            dirtab = self.sequana_factory._sequana_paired_tab
            dirtab._set_paired_filenames([os.path.abspath(f)
                for f in user_options.input_files])
            self.ui.tabs_pipeline.setCurrentIndex(0)
            self.ui.tabWidget.setCurrentIndex(1)

        # We may have set some pipeline, snakefile, working directory
        self._load_and_merge_config()
        self.create_base_form()
        self.fill_until_starting()
        self.update_footer()

    def initUI(self):
        # The logger is not yet set, so we use the module directly
        colorlog.info("Initialising GUI")

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
                customBanner="Welcome to Sequanix embedded ipython console\n" +
                    "The entire GUI interface is stored in the variable gui\n" +
                    "Note also that you can use this interface as a shell \n" +
                    "command line interface preceding your command with ! character\n")
            #self.ipyConsole.printText("The variable 'foo' andion.")
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
        colorlog.getLogger().setLevel(colorlog.logging.logging.INFO)
        self.ui.layout_logger.addWidget(self.logTextBox.widget)

        # Connectors to actions related to the menu bar
        self.ui.actionQuit.triggered.connect(self.menuQuit)
        self.ui.actionHelp.triggered.connect(self.menuHelp)
        self.ui.actionAbout.triggered.connect(self.menuAbout)
        self.ui.actionSnakemake.triggered.connect(self.snakemake_dialog.exec_)
        self.ui.actionPreferences.triggered.connect(self.preferences_dialog.exec_)

        # connectors related to the pipeline tabs (pipeline/generic)
        self.ui.tabs_pipeline.currentChanged.connect(self.update_footer)

        self.set_sequana_pipeline()
        self.set_generic_pipeline()

        # The run/save/dag footer buttons
        self.connect_footer_buttons()

        self.process = QtCore.QProcess(self)
        self.process.started.connect(lambda: self.ui.run_btn.setEnabled(False))
        self.process.started.connect(lambda: self.ui.stop_btn.setEnabled(True))
        self.process.started.connect(lambda: self.ui.unlock_btn.setEnabled(False))
        self.process.started.connect(lambda: self.start_progress)
        self.process.started.connect(lambda: self.ui.save_btn.setEnabled(False))
        self.process.started.connect(lambda: self.ui.tabs_pipeline.setEnabled(False))

        self.process.finished.connect(lambda: self.ui.run_btn.setEnabled(True))
        self.process.finished.connect(lambda: self.ui.stop_btn.setEnabled(False))
        self.process.finished.connect(lambda: self.ui.unlock_btn.setEnabled(True))
        self.process.finished.connect(lambda: self.ui.save_btn.setEnabled(True))
        self.process.finished.connect(lambda:self.ui.tabs_pipeline.setEnabled(True))
        self.process.finished.connect(self.end_run)

        self.process.readyReadStandardOutput.connect(
            self.snakemake_data_stdout)
        self.process.readyReadStandardError.connect(self.snakemake_data_error)

        # This is for the show dag btn. Created here once for all
        self.process1 = QtCore.QProcess(self)
        self.process2 = QtCore.QProcess(self)


        self.ui.tabWidget.currentChanged.connect(lambda: self.ui.run_btn.setEnabled(False))

    #|-----------------------------------------------------|
    #|                       MENU related                  |
    #|-----------------------------------------------------|
    def menuAbout(self):
        from sequana import version
        url = 'sequana.readthedocs.io'
        widget = About()
        widget.setText("Sequana version %s " % version)
        widget.setInformativeText("""
            Online documentation on <a href="http://%(url)s">%(url)s</a>
            <br>
            <br>
            Authors: Thomas Cokelaer and Dimitri Desvillechabrol, 2016
            """ % {"url": url})
        widget.setWindowTitle("Sequana")
        widget.setStandardButtons(QW.QMessageBox.Ok)
        retval = widget.exec_()
        if retval == QW.QMessageBox.Ok:
            widget.close()

    def menuHelp(self):
        url = 'sequana.readthedocs.io'
        msg = About()
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
        There are downloaded automatically with their config file from the Sequana
        library. Here is a typical set of actions to run Sequana pipelines:

        <ol>
        <li> Select a pipeline</li>
        <li> Select the directory or sample tab</li>
            <ul>
           <li> directory: select all fastq.gz files</li>
           <li> samples: select a single-end or paired-end file(s)</li>
            </ul>
        <li> Select the working directory</li>
        </ol>

        <h2> Generic pipelines </h2>

        <ul>
        <li>The working directory must be set.</li>
        <li>The Snakefile must be set.</li>
        <li>Config file may or may not be required depending on the Snakefile
content:</li>
        <ul><li>If the config file is provided, the --configfile is used. </li>
        <li>If configfile keyword is used, the config file should be provided </li>
        </ul></ul>

        <h2> Sequana pipeline dedicated help </help>
             %(pipelines)s
        </div>
        """ % {"url": url, "pipelines": pipelines_text})
        msg.setWindowTitle("Sequana")
        msg.setStandardButtons(QW.QMessageBox.Ok)
        self._dialog_help = msg
        retval = msg.exec_()
        if retval == QW.QMessageBox.Ok:
            msg.close()

    def menuQuit(self):
        self._quit_msg = WarningMessage("Do you really want to quit ?")
        self._quit_msg.setStandardButtons(QW.QMessageBox.Yes | QW.QMessageBox.No)
        self._quit_msg.setDefaultButton(QW.QMessageBox.No)
        quit_answer = self._quit_msg.exec_()
        if quit_answer == QW.QMessageBox.Yes:
            self.close()

    def set_level(self):
        # Set the level of the logging system
        pref = self.preferences_dialog.ui
        level = pref.preferences_options_general_logging_value.currentText()
        level =  getattr(colorlog.logging.logging, level)
        colorlog.getLogger().setLevel(level)

    # ---------------------------------------------------------------
    # More GUI / reading the snakefile (sequana or generic)
    # ---------------------------------------------------------------
    def set_sequana_pipeline(self):
        #
        # Use choice_button.activated so that if one select the same item again,
        # The pipeline connector
        pipelines = sorted(snaketools.pipeline_names)
        self.ui.choice_button.addItems(pipelines)
        self.ui.choice_button.activated[str].connect(self._update_sequana)
        #self.ui.choice_button.activated.connect(self._load_and_merge_config)
        self.ui.choice_button.installEventFilter(self)

        # populate the factory with the choice button
        self.sequana_factory = SequanaFactory(
            combobox=self.ui.choice_button,
            run_button=self.ui.run_btn)
        self.sequana_factory.valid_pipelines = pipelines

        # a local alias
        saf = self.sequana_factory

        # add widgets
        self.ui.layout_sequana_wkdir.addWidget(saf._directory_browser)
        self.ui.layout_sequana_input_dir.addWidget(saf._sequana_directory_tab)
        self.ui.layout_sequana_input_files.addWidget(saf._sequana_paired_tab)

    @QtCore.pyqtSlot(str)
    def _update_sequana(self, index):
        """ Change options form when user change the pipeline."""
        if self.ui.choice_button.findText(index) == 0:
            self.clear_form()
            self.rule_list = []
            self.fill_until_starting()
            return

        self.info("Reading sequana %s pipeline" % index)
        self.create_base_form()
        self._set_focus_on_config_tab()
        self.fill_until_starting()
        self.update_footer()

    def set_generic_pipeline(self):

        self.generic_factory = GenericFactory(self.ui.run_btn)
        gaf = self.generic_factory

        # The config file connectors
        gaf._config_browser.clicked_connect(self._load_and_merge_config)

        # working directory tab
        gaf._directory_browser.clicked_connect(self._load_and_merge_config)


        # Update the main UI with
        self.ui.layout_generic_snakefile.addWidget(gaf._snakefile_browser)
        self.ui.layout_generic_config.addWidget(gaf._config_browser)
        self.ui.layout_generic_wkdir.addWidget(gaf._directory_browser)

    # ---------------------------------------------------------------------
    # Fotter connectors
    # ---------------------------------------------------------------------

    def connect_footer_buttons(self):
        self.ui.run_btn.setEnabled(False)
        self.ui.run_btn.clicked.connect(self.click_run)

        self.ui.stop_btn.clicked.connect(self.click_stop)
        self.ui.stop_btn.setEnabled(False)

        self.ui.unlock_btn.clicked.connect(self.ui.run_btn.setEnabled)
        self.ui.unlock_btn.clicked.connect(self.unlock_snakemake)
        self.ui.unlock_btn.setEnabled(True)

        self.ui.report_btn.setEnabled(True)
        self.ui.report_btn.clicked.connect(self.open_report)

        self.ui.save_btn.clicked.connect(self.save_configfile)
        self.ui.save_btn.clicked.connect(self.save_snakefile)
        self.ui.save_btn.clicked.connect(self.update_footer)

        self.ui.dag_btn.setEnabled(False)
        self.ui.dag_btn.clicked.connect(self.show_dag)

    # -----------------------------------------------------------------
    # functionalities to switch between sequana and generic pipelines
    # -----------------------------------------------------------------

    def _get_mode(self):
        # figure out if we are dealing with a sequana pipeline or not
        index = self.ui.tabs_pipeline.currentIndex()
        if index == 0:
            return "sequana"
        elif index == 1:
            return "generic"
    mode = property(_get_mode)

    def _get_factory(self):
        return getattr(self, "%s_factory" % self.mode)
    factory = property(_get_factory)

    def _get_config(self):
        return getattr(self, "%s_factory" % self.mode).config
    config = property(_get_config)

    def _get_configfile(self):
        return getattr(self, "%s_factory" % self.mode).configfile
    configfile = property(_get_configfile)

    def _get_snakefile(self):
        return getattr(self, "%s_factory" % self.mode).snakefile
    snakefile = property(_get_snakefile)

    def _get_working_dir(self):
        return getattr(self, "%s_factory" % self.mode).directory
    working_dir = property(_get_working_dir)

    # ----------------------------------------------------------------------
    # Snakemake related (config, running)
    # ----------------------------------------------------------------------

    def fill_until_starting(self):
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

        configfile = os.path.basename(self.factory.configfile)
        snakemake_line += ["--configfile", configfile]

        return snakemake_line

    def click_run(self):
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
        if self.working_dir is None:
            self.warning("Set the working directory first")
            return

        # We copy the sequana and genereic snakefile into a filename called
        # Snakefile
        if self.mode == "sequana":
            snakefile = self.working_dir + os.sep + os.path.basename(self.snakefile)

        elif self.mode == "generic":
            snakefile = self.generic_factory.snakefile

        if os.path.exists(snakefile) is False:
            self.critical("%s does not exist" % snakefile)
            return

        snakemake_args = self._get_snakemake_command(snakefile)

        self.info("Starting process with %s " % " ".join(snakemake_args))
        self.process.setWorkingDirectory(self.working_dir)
        self.process.start("snakemake", snakemake_args)

    # -------------------------------------------------------------------
    # Create the base form
    # -------------------------------------------------------------------

    def create_base_form(self):
        """ Create form with all options necessary for a pipeline.

        ::

            ########################################################
            #   valid python docstring to be interepreted by sphinx
            #
            #   section:
            #      item1: 10
            #      item2: 20


        """
        self.rule_list = []
        if self.config is None:
            self.clear_form()
            return

        self.info("Creating form based on config file")
        self.clear_form()
        rules_list = list(self.config._yaml_code.keys())
        rules_list.sort()
        self.necessary_dict = {}

        # !!!!
        # cfg._yaml_code.ca.items['adapter_removal']
        # is a list of 4 items
        # [None, largecomment, shortcomment, None]

        # For each section, we create a widget (RuleForm). For isntance, first,
        # one is accessible asfollows:
        # gui.form.itemAt(0).widget()

        config = self.config
        for count, rule in enumerate(rules_list):
            # Check if this is a dictionnary
            contains = config._yaml_code[rule]

            comments = config.get_section_long_comment(rule)

            if isinstance(contains, dict) and (
                    rule not in SequanaGUI._not_a_rule):
                rule_box = Ruleform(rule, contains, count, self._browser_keyword)
                rule_box.connect_all_option(
                    lambda: self.ui.run_btn.setEnabled(False))
                self.comments = comments
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
        self._set_focus_on_config_tab()

    # ----------------------------------------------------------
    # STOP footer button
    # ----------------------------------------------------------

    def click_stop(self):
        """The stop button"""
        pal = self.ui.progressBar.palette()
        pal.setColor(QtGui.QPalette.Highlight, self._colors['orange'])
        self.ui.progressBar.setPalette(pal)

        if self.process.state() != 0:
            self.info("Process running, stopping it... ")
            # We must use a ctrl+C interruption like so that snakemake
            # handles the interruption smoothly
            self.info("killing the main snakemake process. This may take a few seconds ")
            try:
                os.kill(self.process.pid(), signal.SIGINT)
                time.sleep(4)
            except:
                pass # already stopped ?
            self.info("Process killed successfully.")
        self.ui.save_btn.setEnabled(True)
        self.ui.run_btn.setEnabled(True)
        self.ui.stop_btn.setEnabled(False)
        self.ui.tabs_pipeline.setEnabled(True)

    # --------------------------------------------------------------------
    # Progress bar
    # --------------------------------------------------------------------

    def update_progress_bar(self, line):
        """ Parse with a regex to retrieve current step and total step.
        """
        grouprex = self._step_regex.findall(line)
        # Use last "x of y" (not the first item at position 0)
        if grouprex:
            step = int(grouprex[-1][0]) / float(grouprex[-1][1]) * 100
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

    def save_snakefile(self):
        self.factory._copy_snakefile()

    def save_configfile(self, force=False):
        self.info('Saving project')
        if self.config is None:
            if self.mode == "generic" and self.generic_factory.is_runnable():
                pass
            else:
                msg = WarningMessage("You must choose a pipeline or config file before saving.")
                msg.exec_()
            return

        try:
            form_dict = dict(self.create_form_dict(self.form),
                                 **self.necessary_dict)
        except AttributeError as err:
            self.error(err)
            msg = WarningMessage("You must choose a pipeline before saving.")
            msg.exec_()
            return

        if self._undefined_section in form_dict.keys() and self.mode != "sequana":
            # if sequana, we ignore input_directory and others but in the
            # generic case, we want the entire config file
            del form_dict[self._undefined_section]

        # get samples names or input_directory
        if self.mode == "sequana":
            self.info("saving config file (sequana)")
            flag1 = self.sequana_factory._sequana_directory_tab.get_filenames()
            flag2 = self.sequana_factory._sequana_paired_tab.get_filenames()
            if len(flag1) == 0 and len(flag2) == 0:
                msg = WarningMessage("You must choose an input first.")
                msg.exec_()
                return

            if self.ui.tabWidget.currentIndex() == 0:
                filename = self.sequana_factory._sequana_directory_tab.get_filenames()
                form_dict["input_directory"] = (filename)
            elif self.ui.tabWidget.currentIndex() == 1:
                filename = self.sequana_factory._sequana_paired_tab.get_filenames()
                form_dict["input_samples"] = (filename)

        if self.mode == "generic":
            self.info("saving config file (generic)")

        # Let us update the attribute with the content of the form
        cfg = self.config
        cfg.config.update(form_dict)
        cfg._update_yaml()
        cfg.cleanup()
        self.cfg = cfg

        # We must update the config and then the yaml to keep the comments. This
        # lose the comments:
        # self.configfile._yaml_code.update(form_dict)

        if self.working_dir:
            # FIXME in generic, config.yaml is not necessary the filename
            if self.mode == "sequana":
                yaml_path = self.working_dir + os.sep + "config.yaml"
                self.warning("copy requirements (if any)")
                cfg.copy_requirements(target=self.working_dir)
            elif self.mode == "generic":
                yaml_path = self.working_dir + os.sep + os.path.basename(self.generic_factory.configfile)

            if os.path.isfile(yaml_path) and force is False:
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
                retval = save_msg.exec()
                if retval in [16384, 2048]:
                    self.info("Saving config file (does not exist)")
                    cfg.save(yaml_path, cleanup=False)
                    self.ui.dag_btn.setEnabled(True)
            else:
                self.critical("Saving config file (exist but force it)")
                cfg.save(yaml_path, cleanup=False)
                self.ui.dag_btn.setEnabled(True)
        else:
            self.critical("Config file not saved (no wkdir)")
            msg = WarningMessage("You must set a working directory", self)
            msg.exec_()

        self.ui.run_btn.setEnabled(True)

    def update_footer(self):
        # This functions copies the snakefile in the working directory
        # if possible.

        # Run is on if working dir is on AND
        # 1. for sequana: pipeline is set
        # 2. for generic: snakefile is present irrespective of config file in generic mode
        if self.mode == "sequana":
            if self.working_dir and self.sequana_factory.pipeline:
                # FIXME how to get the input as well
                # self.ui.tabs_pipeline.currentWidget().path_is_setup():
                #self.ui.run_btn.setEnabled(True)
                self.ui.dag_btn.setEnabled(True)
                self.debug('Switching RUN button on')
            else:
                #self.ui.run_btn.setEnabled(False)
                self.ui.dag_btn.setEnabled(False)
                self.debug('Switching RUN button off (missing working directory or pipeline)')
        else: # generic
            if self.generic_factory.is_runnable():
                #self.ui.run_btn.setEnabled(True)
                self.ui.dag_btn.setEnabled(True)
                self.debug('Switching RUN button on')
            else:
                #self.ui.run_btn.setEnabled(False)
                self.ui.dag_btn.setEnabled(False)
                self.debug('Switching RUN button off (missing working directory or pipeline)')
        #return self.ui.run_btn.setEnabled(False)

    # -----------------------------------------------------------------------
    # UNLOCK footer button
    # -----------------------------------------------------------------------

    def unlock_snakemake(self):

        if self.working_dir is None:
            return

        # FIXME this does not work as expected
        self.ui.run_btn.setEnabled(False)
        self.update()
        time.sleep(2)

        if os.path.exists(self.snakefile) is False:
            self.warning("snakefile not found. should not happen")
            return

        self.cmd = ['snakemake', "-s", self.snakefile, "--unlock"]
        self.info("Running " + " ".join(self.cmd))
        self.info("Please wait a second. Unlocking working directory")
        # focus on tab with snakemake output
        self.ui.tabs.setCurrentIndex(0)

        self.ui.tabs_pipeline.setEnabled(False)
        try:
            snakemake_proc = sp.Popen(self.cmd, cwd=self.working_dir)
            snakemake_proc.wait()
        except:
            self.critical("Issue while unlocking the directory")
        finally:
            self.ui.tabs_pipeline.setEnabled(True)

        self.info("unlocking done")
        self.output.append('<font style="color:brown">Unlocking working directory</font>')

        self.ui.run_btn.setEnabled(True)
        self.ui.stop_btn.setEnabled(False)

    # -----------------------------------------------------------------------
    # DAG footer button
    # -----------------------------------------------------------------------

    def show_dag(self):
        self.info("Creating DAG image.")

        # The config must have been saved, so we just need to copy it

        # copy the pipeline in the working directory
        # FIXME for generic case, not required
        self.copy(self.snakefile, self.working_dir)

        snakefile = os.path.basename(self.snakefile)

        svg_filename = self._tempdir.path() + os.sep + "test.svg"

        snakemake_line = ["snakemake", "-s", snakefile]
        snakemake_line += ["--rulegraph"]
        if self.mode == "generic" and self.configfile:
            # make sure to copy the config file
            snakemake_line += ["--configfile"]
            snakemake_line += [os.path.basename(self.generic_factory.configfile)]

        snakemake_line += self.get_until_starting_option()

        self.info(snakemake_line)
        self.process1.setWorkingDirectory(self.working_dir)
        self.process1.setStandardOutputProcess(self.process2)

        self.process1.start("snakemake", snakemake_line[1:])
        self.process2.start("dot", ["-Tsvg", "-o", svg_filename])

        self.process1.waitForFinished(50000)
        self.process2.waitForFinished(50000)

        if os.path.exists(svg_filename):
            self.diag = SVGDialog(svg_filename)
            self.diag.show()
        else:
            msg = 'Could not create the DAG file.'
            error = str(self.process1.readAllStandardError())
            msg = CriticalMessage(msg, error)
            msg.exec_()
            return

    def open_report(self):

        pref = self.preferences_dialog.ui
        filename = pref.preferences_options_general_htmlpage_value.text()
        if filename == "":
            filename = QW.QFileDialog.getOpenFileNames(self,
               "Select your HTML report", ".",
                "HTML files (*.html)")[0]
            if len(filename) and os.path.exists(filename[0]):
                filename = filename[0]
            else:
                self.warning("No valid HTML selected and none specified in the preferences.")
                return
        else: # we have a filename hardcoded in the preferences
            if self.working_dir is None:
                self.error("Working directory not set yet")
                return

            filename = self.working_dir + os.sep + filename
            if os.path.exists(filename) is False:
                self.error("%s page does not exist. Check the preferences dialog." % filename)
                return
            else:
                self.info("Reading and openning %s" % filename)

        dialog = self.preferences_dialog.ui # an alias
        url = "file://" + filename

        # The browser executable itself
        browser = dialog.preferences_options_general_browser_value.currentText()

        if browser == "pyqt5":
            self.browser = Browser(url)
            self.browser.show()
        else:
            try:
                easydev.execute("%s %s" % (browser, url), showcmd=False,
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

    def _load_and_merge_config(self):
        self.info("Entering Load and Merge method")
        # If working directory is not set yet, nothing to load
        if self.configfile is None:
            self.clear_form()
            return

        # If pipeline not set, we cannot read the config file
        if self.mode == "sequana" and self.sequana_factory.pipeline is None:
            self.critical("Are we here at any time ?")
            msg = WarningMessage("Please, choose a pipeline first", self)
            msg.exec_()
            return

        if self.working_dir is None:
            return

        config_file = self.working_dir + os.sep + os.path.basename(self.configfile)
        if os.path.isfile(config_file):
            self.warning("An existing config.yaml file already exists in the target directory")
        else:
            # no config to merge, we will use the pipeline config file
            return

        # this should always work but who knows
        try:
            cfg = snaketools.SequanaConfig(config_file)
            cfg.cleanup() # set all empty strings and %()s to None
            config_dict = cfg.config
        except Exception as err:
            self.warning(err)
            self.critical("Could not interpret the sequana config file")
            self.warning("config_dict is None")
            self.clear_form()
            return

        if set(self.config._yaml_code.keys()) == set(config_dict.keys()):
            msg = QW.QMessageBox(
                QW.QMessageBox.Question, "Question",
                "A config file already exist in the working directory.\n" +
                "%s.\n" % self.working_dir +
                "Do you want to import its content (if not, it will be replaced)?",
                QW.QMessageBox.Yes | QW.QMessageBox.No,
                self, Qt.Dialog | Qt.CustomizeWindowHint)
            # Yes == 16384
            if msg.exec_() == 16384:
                self.config._yaml_code.update(config_dict)
                self.create_base_form()
                self.fill_until_starting() #self.rule_list)
            else:
                # if we do not want to import the existing file

                self.create_base_form()
                self.fill_until_starting() #self.rule_list)
                self.critical("fixme")
        else:
            self.critical("The config file that already exists is different. Nothing done")
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
        self.info("Reading settings")
        settings = QtCore.QSettings("sequana_gui", "mainapp")
        if settings.value("tab_position") is not None:
            index = settings.value("tab_position")
            self.ui.tabs_pipeline.setCurrentIndex(int(index))

        if settings.value("tab_generic_position") is not None:
            index = settings.value("tab_generic_position")
            self.ui.tabs_generic.setCurrentIndex(int(index))

        if settings.value("tab_sequana_position") is not None:
            index = settings.value("tab_sequana_position")
            self.ui.tabs_sequana.setCurrentIndex(int(index))

        if settings.value("tab_sequana_input_position") is not None:
            index = settings.value("tab_sequana_input_position")
            self.ui.tabWidget.setCurrentIndex(int(index))

    def write_settings(self):
        settings = QtCore.QSettings("sequana_gui", "mainapp")

        # tab snakemake output/logger/ipython
        index = self.ui.tabs_pipeline.currentIndex()
        settings.setValue("tab_position", index)

        index = self.ui.tabs_generic.currentIndex()
        settings.setValue("tab_generic_position", index)

        index = self.ui.tabs_sequana.currentIndex()
        settings.setValue("tab_sequana_position", index)

        index = self.ui.tabWidget.currentIndex()
        settings.setValue("tab_sequana_input_position", index)

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


class Options(argparse.ArgumentParser):
    def __init__(self, prog="sequana_gui"):
        usage = """dkfj skfk"""
        description = """"""
        super(Options, self).__init__(usage=usage, prog=prog,
            description=description, formatter_class=easydev.SmartFormatter)
        group = self.add_argument_group("GENERAL")
        group.add_argument("-w", "--working-directory", dest="wkdir",
            help="Set working directory", default=None)
        group.add_argument("-n", "--no-splash", dest="nosplash",
            action="store_true",
            help="No splash screen")

        group = self.add_argument_group("SEQUANA")
        group.add_argument("-p", "--pipeline", dest="pipeline",
            default=None,
            help="A valid sequana pipeline name")

        group_mut = group.add_mutually_exclusive_group()
        group_mut.add_argument("-i", "--input-directory", dest="input_directory",
            default=None,
            help="input directory where to find the input data")
        group_mut.add_argument("-f", "--input-files", dest="input_files",
            default=None, nargs="*",
            help="input files")

        group = self.add_argument_group("GENERIC PIPELINES")
        group.add_argument("-s", "--snakefile", dest="snakefile",
            default=None,
            help="A valid Snakefile")
        group.add_argument("-c", "--configfile", dest="configfile",
            default=None,
            help="optional config file to be used by the Snakefile")


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

