# coding: utf-8
import os
import sys
import time
import tempfile
import shutil
import subprocess as sp
import multiprocessing

from PyQt5 import QtCore, QtGui
from PyQt5 import QtWidgets as QW
from PyQt5.Qt import QTemporaryDir
from PyQt5.QtCore import Qt
from PyQt5.QtSvg import QSvgWidget
from PyQt5.QtWebKitWidgets import QWebView

from sequana import snaketools, sequana_data
from sequana.snaketools import Module
from sequana.gui.browser import MyBrowser
from sequana.gui.ipython import QIPythonWidget


import signal
def sigint_handler(*args):
    """Handler for the SIGINT signal."""
    sys.stderr.write('\r')
    if QW.QMessageBox.question(None, '', "Are you sure you want to quit?",
                            QW.QMessageBox.Yes | QW.QMessageBox.No,
                            QW.QMessageBox.No) == QW.QMessageBox.Yes:
        QW.QApplication.quit()


class SequanaGUI(QW.QWidget):
    """ Sequana GUI !
    """

    _not_a_rule = {"requirements", "gatk_bin", "input_directory", "input_samples", "input_pattern",}
    _browser_keyword = {"reference"}

    def __init__(self, ipython=True):
        super().__init__()
        # some variables
        self._setup_ipython = ipython
        self._tempdir = QTemporaryDir()
        self.shell = ""
        self.shell_error = ""

        # The UI layout
        self.initUI()

        self._colors = {
            'green': QtGui.QColor(0,102,0),
            'red': QtGui.QColor(255,0,0),
            'blue': QtGui.QColor(0,0,255),
        }

    def initUI(self):
        # snakemake cluster/local option dialog windows
        self.snakemake_dialog = SnakemakeOptionDialog(self)

        # create menu bar
        self.create_menu_bar()

        # box to choose the pipeline
        self.sequana_config = None
        self.pipeline_is_chosen = False
        choice_layout = QW.QHBoxLayout()
        self.create_choice_button()
        choice_layout.addWidget(self.choice_button)

        # create browser tab to choice directory or files
        self.create_tabs_browser()
        self.tabs_browser.currentChanged.connect(self.switch_run)

        # select the working directory
        groupbox_layout = QW.QHBoxLayout()
        self.working_dir = FileBrowser(directory=True)
        groupbox_layout.addWidget(self.working_dir)
        groupbox = QW.QGroupBox("Working directory")
        groupbox.setContentsMargins(0, 5, 0, 0)
        groupbox.setLayout(groupbox_layout)

        # "until" and "starting" combobox
        control_widget = QW.QGroupBox("Pipeline control")
        control_widget.setContentsMargins(0, 3, 0, -3)
        control_layout = QW.QVBoxLayout(control_widget)
        control_layout.setSpacing(0)
        self.until_box = ComboBoxOption("Until")
        self.starting_box = ComboBoxOption("Starting")
        control_layout.addWidget(self.starting_box)
        control_layout.addWidget(self.until_box)

        # connect fuction on working_dir browser
        self.working_dir.clicked_connect(self.check_existing_config)
        self.working_dir.clicked_connect(self.switch_run)

        # formular which contains all options the pipeline chosen
        widget_formular = QW.QWidget()
        self.formular = QW.QVBoxLayout(widget_formular)
        self.formular.setSpacing(10)
        scroll_area = QW.QScrollArea()
        scroll_area.setWidget(widget_formular)
        scroll_area.setWidgetResizable(True)
        scroll_area.setMinimumHeight(300)


        # main layout
        vlayout = QW.QVBoxLayout(self)
        vlayout.insertSpacing(0, 10)

        # ipython widget
        if self._setup_ipython:
            self.ipyConsole = QIPythonWidget(
                customBanner="Welcome to the embedded ipython console\n")
            self.ipyConsole.pushVariables({"x": 10})
            self.ipyConsole.printText("The variable 'foo' andion.")
            self.ipyConsole.execute("from sequana import *")
            self.ipyConsole.pushVariables({"gui": self})

        # add widgets in layout
        vlayout.addLayout(choice_layout)
        vlayout.addWidget(self.tabs_browser)
        vlayout.addWidget(groupbox)
        vlayout.addWidget(control_widget)
        vlayout.addWidget(scroll_area)
        vlayout.addWidget(self.create_footer_button())
        if self._setup_ipython:
            vlayout.addWidget(self.ipyConsole)

        # Run snakemake/sequana
        self.process = QtCore.QProcess(self)
        self.process.readyRead.connect(self.snakemake_data)

        self.process.started.connect(lambda: self.run_btn.setEnabled(False))
        self.process.started.connect(lambda: self.stop_btn.setEnabled(True))
        self.process.started.connect(lambda: self.start_progress)

        self.process.finished.connect(lambda: self.run_btn.setEnabled(True))
        self.process.finished.connect(lambda: self.stop_btn.setEnabled(False))
        self.process.finished.connect(lambda: self.end_progress)

        self.output = QW.QTextEdit()
        self.output.setWindowTitle("logger")

        self.progressBar = QW.QProgressBar(self)
        self.progressBar.setToolTip("""<p>Progress of the pipeline. color codes:
            <ul>
                <li style="color:red">Red: an error occured</li>
                <li style="color:green">Green: completed with success</li>
                <li style="color:blue">Blue: in progress</li>
            </ul>
            </p>""")

        vlayout.addWidget(self.progressBar)
        vlayout.addWidget(self.output)

        self._vlayout = vlayout
        # main window options
        self.setGeometry(100, 100, 500, 600)
        self.setWindowTitle("Sequana")

        self.show()

        if self._setup_ipython:
            if self.ipyConsole.isHidden() is False:
                self.ipyConsole.hide()

    def quit(self):
        quit_msg = WarningMessage("Do you really want to quit ?")
        quit_msg.setStandardButtons(QW.QMessageBox.Yes | QW.QMessageBox.No)
        quit_msg.setDefaultButton(QW.QMessageBox.No)
        quit_answer = quit_msg.exec_()
        if quit_answer == QW.QMessageBox.Yes:
            self.close()

    def help(self):
        url = 'gdsctools.readthedocs.io'
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
        # msg.setDetailedText("The details are as follows:")
        msg.setStandardButtons(QW.QMessageBox.Ok)
        # msg.buttonClicked.connect(self.msgbtn)
        self._msg_help = msg
        retval = msg.exec_()
        #msg.show()

    def about(self):
        from sequana import version
        url = 'gdsctools.readthedocs.io'
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
        # msg.setDetailedText("The details are as follows:")
        msg.setStandardButtons(QW.QMessageBox.Ok)
        # msg.buttonClicked.connect(self.msgbtn)
        retval = msg.exec_()

    def _read_config(self, config_file):
        try:
            self.sequana_config = snaketools.SequanaConfig(config_file)
            self.sequana_config.cleanup() # set all empty strings and %()s to None
        except AssertionError:
            print("Warning: could not parse the config file")
            return
        self.create_base_formular()

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

    def create_menu_bar(self):
        """ Create menu bar.
        """
        # action quit
        quitAction = QW.QAction("Quit", self)
        quitAction.setShortcut('Ctrl+Q')
        quitAction.setStatusTip('Quit Sequana')
        quitAction.triggered.connect(self.quit)

        # action About
        aboutAction = QW.QAction("About", self)
        aboutAction.setShortcut('Ctrl+A')
        aboutAction.triggered.connect(self.about)

        # action Help
        helpAction = QW.QAction("Help", self)
        helpAction.setShortcut('Ctrl+H')
        helpAction.triggered.connect(self.help)

        # action Options
        optionAction = QW.QAction("Snakemake Option", self)
        optionAction.setShortcut('Ctrl+O')
        optionAction.triggered.connect(self.snakemake_dialog.exec_)

        # ipython switch
        if self._setup_ipython:
            optionIPython = QW.QAction("Show/Hide IPython dialog", self,
                checkable=True)
            optionIPython.setShortcut('Ctrl+D')
            optionIPython.triggered.connect(self.switch_ipython)

        # action import config file
        importAction = QW.QAction("Import", self)
        importAction.setShortcut('Ctrl+I')
        importAction.triggered.connect(self.import_config)

        # set menu bar
        menu_bar = QW.QMenuBar(self)
        menu_bar.setMinimumWidth(500)
        options_menu = menu_bar.addMenu("&File")
        options_menu.addAction(importAction)
        options_menu.addAction(quitAction)

        options_menu = menu_bar.addMenu("&Option")
        options_menu.addAction(optionAction)
        if self._setup_ipython:
            options_menu.addAction(optionIPython)

        options_menu = menu_bar.addMenu("&Help")
        options_menu.addAction(helpAction)

        options_menu = menu_bar.addMenu("&About")
        options_menu.addAction(aboutAction)

    @QtCore.pyqtSlot(bool)
    def switch_ipython(self):
        if self.ipyConsole.isHidden():
            self.ipyConsole.show()
        else:
            self.ipyConsole.hide()

    def create_choice_button(self):
        """ Create button to select the wished pipeline.
        """
        self.choice_button = QW.QComboBox()
        snaketools.pipeline_names.sort()
        self.choice_button.addItems(["select pipeline"] +
                                    snaketools.pipeline_names)
        self.choice_button.currentIndexChanged[str].connect(
            self.on_pipeline_choice)
        self.choice_button.installEventFilter(self)

    @QtCore.pyqtSlot(str)
    def on_pipeline_choice(self, index):
        """ Change options formular when user change the pipeline.
        """
        config_file = snaketools.Module(index)._get_config()
        self.choice_button.removeItem(
            self.choice_button.findText("select pipeline"))

        self._read_config(config_file)
        self.fill_combobox(self.rule_list)
        self.pipeline_is_chosen = True
        self.switch_run()

    def fill_combobox(self, rules_list):
        """ Fill combobox with available rules.
        """
        active_list = [w.get_name() for w in self.rule_list if w.get_do_rule()]
        self.until_box.add_items(active_list)
        self.starting_box.add_items(active_list)

    def create_base_formular(self):
        """ Create formular with all options necessary for a pipeline.
        """
        self.clear_layout(self.formular)
        rules_list = list(self.sequana_config._yaml_code.keys())
        rules_list.sort()
        self.necessary_dict = {}
        self.rule_list = []
        for count, rule in enumerate(rules_list):
            # Check if this is a dictionnary
            contains = self.sequana_config._yaml_code[rule]
            if isinstance(contains, dict) and (
                    rule not in SequanaGUI._not_a_rule):
                rule_box = RuleFormular(rule, contains, count)
                self.formular.addWidget(rule_box)
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

    def create_import_config(self):
        pass

    def create_tabs_browser(self):
        """ Generate file browser widget.
        """
        # create browser file
        fastq_filter = "Fastq file (*.fastq *.fastq.gz *.fq *.fq.gz)"
        paired_tab = FileBrowser(paired=True, file_filter=fastq_filter)
        directory_tab = FileBrowser(directory=True, file_filter=fastq_filter)
        paired_tab.clicked_connect(self.switch_run)
        directory_tab.clicked_connect(self.switch_run)
        # create tab box
        self.tabs_browser = QW.QTabWidget()
        self.tabs_browser.addTab(directory_tab, "Directory")
        self.tabs_browser.addTab(paired_tab, "Sample")

    def create_footer_button(self):
        """ Create Run/Save/Quit buttons
        """
        self.run_btn = QW.QPushButton("Run")
        self.run_btn.setEnabled(False)
        self.run_btn.clicked.connect(self.start_sequana)
        self.run_btn.setShortcut("Ctrl+R")
        self.run_btn.setToolTip("<p>Run the pipeline (shortcut: Ctrl+)</p>")

        self.stop_btn = QW.QPushButton("Stop")
        self.stop_btn.clicked.connect(self.click_stop)
        self.stop_btn.setEnabled(False)
        self.run_btn.setToolTip("<p>Stop the running pipeline</p>")

        self.unlock_btn = QW.QPushButton("Unlock")
        self.unlock_btn.setShortcut("Ctrl+U")
        self.unlock_btn.clicked.connect(self.unlock_snakemake)
        self.unlock_btn.setEnabled(True)
        self.run_btn.setToolTip("<p>Unlock the directory where the pipeline is run</p>")

        self.report_btn = QW.QPushButton("Open Report")
        self.report_btn.setEnabled(True)
        self.report_btn.clicked.connect(self.open_report)

        self.save_btn = QW.QPushButton("Save Config")
        self.save_btn.clicked.connect(self.save_config_file)

        self.dag_btn = QW.QPushButton("Show DAG")
        self.dag_btn.setEnabled(False)
        self.dag_btn.setToolTip("""<p>Pressing this button, a DAG is created
                                 and shown. This is a good way to check your
                                 config file </p>""" )
        self.dag_btn.clicked.connect(self.show_dag)

        footer_widget = QW.QWidget()
        footer_layout = QW.QHBoxLayout(footer_widget)
        footer_layout.addWidget(self.run_btn)
        footer_layout.addWidget(self.stop_btn)
        footer_layout.addWidget(self.unlock_btn)
        footer_layout.addWidget(self.report_btn)
        footer_layout.addWidget(self.save_btn)
        footer_layout.addWidget(self.dag_btn)
        return footer_widget

    def show_dag(self):
        if self.pipeline_is_chosen:
            snakefile = Module(self.choice_button.currentText()).snakefile
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
        if self.pipeline_is_chosen and self.working_dir.get_filenames():
            filename = self.working_dir.get_filenames() + "/multi_summary.html"
            if os.path.exists(filename) is False:
                WarningMessage("""multi_summary.html not found.
                Most probably the analysis did not finish correctly""")
                return
            url = "file://" + filename

            self.browser = MyBrowser()
            self.browser.load(QtCore.QUrl(url))
            self.browser.show()
        else:
            msg = WarningMessage("no working directory selected yet")
            msg.exec_()

    def snakemake_data(self):
        cursor = self.output.textCursor()
        cursor.movePosition(cursor.End)
        data = str(self.process.readAllStandardOutput())
        error = str(self.process.readAllStandardError())
        self.shell += data
        self.shell_error += error

        for this in data.split("\\n"):
            line = this.strip()
            if line and len(line) > 3 and "complete in" not in line: # prevent all b'' strings
                line = line.replace("\\r","")
                line = line.replace("\\t","    ")
                cursor.insertHtml('<p style="color:blue">' + line +'</p><br>')
                cursor.movePosition(cursor.End)

        for this in error.split("\\n"):
            line = this.strip()
            if line and len(line) > 3 and "complete in" not in line: # prevent all b'' strings
                line = line.replace("\\r","")
                line = line.replace("\\t","    ")
                cursor.insertHtml('<p style="color:red">' + line +'</p><br>')
                cursor.movePosition(cursor.End)

        step = [x for x in self.shell_error.split("\\n") if "steps" in x]
        if len(step):
            step = step[-1]
            start, end = step.split(" of ")
            try:
                start = int(start.strip() )
                end = int(end.strip().split(" ")[0].strip())
                step = int(start) / float(end) * 100
                if step<1:
                    step=1
                self.progressBar.setValue(step)
            except:
                pass

    def get_until_starting_option(self):
        """ Return list with starting rule and end rule.
        """
        until_rule = self.until_box.get_value()
        starting_rule = self.starting_box.get_value()
        option = []
        if until_rule:
            option += ["--no-hooks", "-U", until_rule]
        if starting_rule:
            option += ["-R", starting_rule]
        return option

    def _get_snakemake_command(self, snakefile):
        snakemake_line = ["-s", snakefile, "--stat", "stats.txt"]
        snakemake_line += self.snakemake_dialog.get_snakemake_options()
        options = self.get_until_starting_option()
        snakemake_line += options
        return snakemake_line

    def start_sequana(self):
        pal = self.progressBar.palette()
        pal.setColor(QtGui.QPalette.Highlight, self._colors['blue'])
        self.progressBar.setPalette(pal)
        self.progressBar.setValue(1)

        # Prepare the command and working directory.
        working_dir = self.working_dir.get_filenames()
        self.save_config_file()
        rules = self.get_rules(self.formular)
        snakefile = Module(self.choice_button.currentText()).snakefile

        new_snakefile = working_dir + os.sep + os.path.basename(snakefile)
        if os.path.exists(new_snakefile) is False:
            shutil.copy(snakefile, working_dir)
        snakefile = new_snakefile

        snakemake_args = self._get_snakemake_command(snakefile)

        self.process.setWorkingDirectory(working_dir)
        self.process.start("snakemake", snakemake_args)
        self.process.readyRead.connect(self.snakemake_data)
        self.process.finished.connect(self.end_run)
        #self.process.waitForFinished()

    def unlock_snakemake(self):
        working_dir = self.working_dir.get_filenames()
        snakefile = Module(self.choice_button.currentText()).snakefile
        new_snakefile = working_dir + os.sep + os.path.basename(snakefile)
        if os.path.exists(new_snakefile) is False:
            print("config not found. should not happen")
        else:
            snakefile = new_snakefile
            self.cmd = ['snakemake', "-s", snakefile, "--unlock"]
            self.cwd = working_dir
            print(self.cmd)
            print("Please wait a second")
            snakemake_proc = sp.Popen(self.cmd,
                cwd=os.path.basename(working_dir))
        self.run_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)
        #snakemake_proc.communicate()

    def click_stop(self):
        print("stopped manually")
        pal = self.progressBar.palette()
        pal.setColor(QtGui.QPalette.Highlight, self._colors['red'])
        self.progressBar.setPalette(pal)
        self.process.finished.disconnect() # we do not want to call end_run()
        self.run_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)
        self.process.close()

    def end_run(self):
        pal = self.progressBar.palette()
        if self.progressBar.value() >= 100 :
            pal.setColor(QtGui.QPalette.Highlight, self._colors['green'])
            self.progressBar.setPalette(pal)
            print('Run done. Status: successful')
        else:
            pal.setColor(QtGui.QPalette.Highlight, self._colors['red'])
            self.progressBar.setPalette(pal)
            text = 'Run manually to check the exact error or check the log'
            msg = CriticalMessage(text, self.process.readAllStandardError())
            msg.exec_()
            return

    def start_sequana2(self):
        working_dir = self.working_dir.get_filenames()
        self.save_config_file()
        rules = self.get_rules(self.formular)
        snakefile = Module(self.choice_button.currentText()).snakefile

        new_snakefile = working_dir + os.sep + os.path.basename(snakefile)
        if os.path.exists(new_snakefile) is False:
            shutil.copy(snakefile, working_dir)
        snakefile = new_snakefile

        snakemake_args = self._get_snakemake_command(snakefile)

        self.cmd = ['snakemake'] + snakemake_args
        self.cwd = working_dir
        snakemake_proc = sp.Popen(self.cmd,
            cwd=os.path.basename(working_dir))
        #snakemake_proc.communicate()

    def start_progress(self):
        self.progressBar.setRange(0,1)

    def end_progress(self):
        self.progressBar.setValue(100)
        QtGui.QMessageBox.information(self, "Done")

    def switch_run(self):
        if self.working_dir.path_is_setup():
            if self.tabs_browser.currentWidget().path_is_setup():
                if self.working_dir.path_is_setup():
                    if self.pipeline_is_chosen:
                        return self.run_btn.setEnabled(True)
        return self.run_btn.setEnabled(False)

    def save_config_file(self):
        try:
            formular_dict = dict(self.create_formular_dict(self.formular),
                                 **self.necessary_dict)
        except AttributeError:
            msg = WarningMessage("You must choose a pipeline before saving.")
            msg.exec_()
            return

        # get samples names or input_directory
        if self.tabs_browser.currentIndex() == 1:
            formular_dict["samples"] = (
                self.tabs_browser.currentWidget().get_filenames())
        else:
            formular_dict["input_directory"] = (
                self.tabs_browser.currentWidget().get_filenames())

        # Let us update tha attribute with the content of the form
        self.sequana_config._yaml_code.update(formular_dict)

        if self.working_dir.path_is_setup():
            yaml_path = self.working_dir.get_filenames() + "/config.yaml"
            self.sequana_config.copy_requirements(target=self.working_dir.get_filenames())

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
                    self.dag_btn.setEnabled(True)
            else:
                self.sequana_config.cleanup()
                self.sequana_config.save(yaml_path)
                self.dag_btn.setEnabled(True)
        else:
            msg = WarningMessage("You must set a working directory", self)
            msg.exec_()

    def create_formular_dict(self, layout):
        def _cleaner(value):
            # This is to save the YAML file correctly since the widgets tend to
            # convert None and empty strings as '""' or "''"
            if value in ['None', None, '', '""', "''"]:
                return None
            else:
                return value

        widgets = (layout.itemAt(i).widget() for i in range(layout.count()))
        formular_dict = {w.get_name(): _cleaner(w.get_value()) if w.is_option()
                         else self.create_formular_dict(w.get_layout())
                         for w in widgets}
        return formular_dict

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
                "Do you want to merge this file ?",
                QW.QMessageBox.Yes | QW.QMessageBox.No,
                self, Qt.Dialog | Qt.CustomizeWindowHint)
            if msg.exec_() == 16384:
                self.sequana_config._yaml_code.update(config_dict)
                self.create_base_formular()
                self.fill_combobox(self.rule_list)
        return True

    def get_rules(self, layout):
        widgets = (layout.itemAt(i).widget() for i in range(layout.count()))
        rules = [w.get_name() for w in widgets]
        return rules

    def eventFilter(self, source, event):
        """ Inactivate wheel event of combobox
        """
        if event.type() == QtCore.QEvent.Wheel and source is self.choice_button:
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


class About(QW.QMessageBox):
    """A resizable QMessageBox for the About dialog"""
    def __init__(self, *args, **kwargs):
        super(About, self).__init__(*args, **kwargs)
        self.setSizeGripEnabled(True)

    def event(self, e):
        result = super(About, self).event(e)

        self.setMinimumHeight(0)
        self.setMaximumHeight(16777215)
        self.setMinimumWidth(500)
        self.setMaximumWidth(16777215)
        self.setSizePolicy(QW.QSizePolicy.Expanding, QW.QSizePolicy.Expanding)

        textEdit = self.findChild(QW.QTextEdit)
        if textEdit is not None:
            textEdit.setMinimumHeight(0)
            textEdit.setMaximumHeight(16777215)
            textEdit.setMinimumWidth(0)
            textEdit.setMaximumWidth(16777215)
            textEdit.setSizePolicy(QW.QSizePolicy.Expanding,
                                   QW.QSizePolicy.Expanding)

        return result

#    # We only need to extend resizeEvent, not every event.
#    def resizeEvent(self, event):
#        result = super(About, self).resizeEvent(event)
#        details_box = self.findChild(QW.QTextEdit)
#        if details_box is not None:
#            details_box.setFixedSize(details_box.sizeHint())
#        return result



class FileBrowser(QW.QWidget):
    """ Class to create a file browser in PyQT5.
    """
    def __init__(self, paired=False, directory=False, file_filter=None):
        super().__init__()
        # Set filter for file dialog
        self.filter = "Any file (*)"
        if file_filter is not None:
            self.filter = file_filter + ";;" + self.filter
        self.empty_msg = "No file selected"
        self.btn = QW.QPushButton("Browse")
        self.btn.setFixedSize(100,20)

        # Add default color
        self.btn.setStyleSheet("QPushButton {background-color: #AA0000; "
                               "color: #EEEEEE}")

        if directory:
            self.empty_msg = "No directory selected"
            self.btn.clicked.connect(self.browse_directory)
        elif paired:
            self.btn.clicked.connect(self.browse_paired_file)
        else:
            self.btn.clicked.connect(self.browse_file)
        self.btn_filename = QW.QLabel(self.empty_msg)
        self.set_empty_path()
        widget_layout = QW.QHBoxLayout(self)
        widget_layout.setContentsMargins(0, 3, 0, 3)
        widget_layout.addWidget(self.btn)
        widget_layout.addWidget(self.btn_filename)

    def _setup_true(self):
        self.setup = True

    def setup_color(self):
        if self.path_is_setup():
            self.btn.setStyleSheet("QPushButton {background-color: #00AA00; "
                                   "color: #EEEEEE}")
        else:
            self.btn.setStyleSheet("QPushButton {background-color: #AA0000; "
                                   "color: #EEEEEE}")

    def browse_paired_file(self):
        file_path = QW.QFileDialog.getOpenFileNames(self, "Select a sample", ".",
                                                 self.filter)[0]
        if not file_path:
            self.set_empty_path()
        elif len(file_path) > 2:
            msg = WarningMessage("You must pick only one sample", self)
            self.set_empty_path()
            msg.exec_()
        else:
            self.paths = {"file{0}".format(i+1): file_path[i]
                          for i in range(0, len(file_path))}
            self.btn_filename.setText("\n".join([key + ": " + value
                                      for key, value in self.paths.items()]))
            self._setup_true()
            self.setup_color()

    def browse_directory(self):
        dialog = DirectoryDialog(self, "Select a directory", ".", self.filter)
        directory_path = dialog.get_directory_path()
        if directory_path:
            self.set_filenames(directory_path)
        else:
            self.set_empty_path()

    def browse_file(self):
        try:
            file_path = QW.QFileDialog.getOpenFileNames(self, "Single File", ".",
                                                     self.filter)[0][0]
            self.set_filenames(file_path)
        except IndexError:
            self.set_empty_path()

    def get_filenames(self):
        return self.paths

    def set_filenames(self, filename):
        self.paths = filename
        if len(filename) > 20:
            self.btn_filename.setText("...." + filename[-20:])
        else:
            self.btn_filename.setText(filename)
        self._setup_true()
        self.setup_color()

    def set_empty_path(self):
        self.btn_filename.setText(self.empty_msg)
        self.paths = ""
        self.setup = False
        self.setup_color()

    def set_enable(self, switch_bool):
        if switch_bool:
            self.setup_color()
        else:
            self.btn.setStyleSheet("QPushButton {background-color: #AAAAAA; "
                                   "color: #222222}")
        self.btn.setEnabled(switch_bool)

    def path_is_setup(self):
        return self.setup

    def clicked_connect(self, function):
        """ Connect additionnal function on browser button. It is used to
        activate run button in Sequana GUI.
        """
        self.btn.clicked.connect(function)


class RuleFormular(QW.QGroupBox):
    do_option = "do"
    def __init__(self, rule_name, rule_dict, count=0):
        super().__init__(rule_name)

        # to handle recursive case
        self.do_widget = None

        self.rule_name = rule_name
        self.rule_dict = rule_dict
        self.layout = QW.QVBoxLayout(self)
        self.layout.setSpacing(2)
        self.setAutoFillBackground(True)

        for option, value in self.rule_dict.items():
            if option.endswith("_directory"):
                option_widget = FileBrowserOption(option, value,
                                                  directory=True)
            elif option.endswith("_file"):
                option_widget = FileBrowserOption(option, value,
                                                  directory=False)
            elif option in SequanaGUI._browser_keyword:
                option_widget = FileBrowserOption(option, value,
                                                  directory=False)
            elif isinstance(value, bool) or option=="do":
                # for the do option, we need to check its value
                if value in ["yes", "YES", "True", "TRUE"]:
                    value = True
                elif value in ["no", "NO", "False", "FALSE"]:
                    value = False
                else:
                    print("Incorrect value found in config file for 'do'")
                option_widget = BooleanOption(option, value)
                if option == RuleFormular.do_option:
                    self.do_widget = option_widget
                    option_widget.connect(self._widget_lock)
            else:
                try:
                    option_widget = NumberOption(option, value)
                except TypeError:
                    try:
                        option_widget = TextOption(option, value)
                    except TypeError:
                        option_widget = RuleFormular(option, value)
            self.layout.addWidget(option_widget)
        try:
            self._widget_lock(self.do_widget.get_value())
        except AttributeError:
            pass

    def get_name(self):
        return self.rule_name

    def get_layout(self):
        return self.layout

    def get_do_rule(self):
        """ If there are no "do widget", rules must be done. Else, it gets value
        of check box.
        """
        if self.do_widget is None:
            return True
        else:
            return self.do_widget.get_value()

    def is_option(self):
        return False

    def connect_do(self, task):
        if self.do_widget:
            self.do_widget.connect(task)

    def _widget_lock(self, switch_bool):
        widget_list = (self.layout.itemAt(i).widget() for i in
                       range(self.layout.count()))
        for w in widget_list:
            if w is self.do_widget:
                continue
            w.set_enable(switch_bool)

    def set_enable(self, switch_bool):
        self._widget_lock(switch_bool)


class GeneralOption(QW.QWidget):
    """ Parent class for Options. It defines design of options
    """
    def __init__(self, option):
        super().__init__()
        self.option = option
        self.layout = QW.QHBoxLayout(self)
        self.layout.setContentsMargins(0, 3, 0, 3)
        self.layout.addWidget(QW.QLabel(option))

    def get_name(self):
        return self.option

    def is_option(self):
        return True

    def get_value(self):
        pass

    def get_tuple(self):
        return (self.get_name(), self.get_value())

    def set_enable(self):
        pass


class BooleanOption(GeneralOption):
    """ Wrapp QCheckBox class
    """
    def __init__(self, option, value):
        super().__init__(option)

        self.check_box = QW.QCheckBox()
        self.check_box.setChecked(value)

        self.answer = QW.QLabel()
        self.switch_answer()

        self.check_box.clicked.connect(self.switch_answer)

        self.layout.addWidget(self.check_box)
        self.layout.addWidget(self.answer)

    def get_value(self):
        return self.check_box.isChecked()

    def set_enable(self, switch_bool):
        self.check_box.setEnabled(switch_bool)

    def switch_answer(self):
        value = self.get_value()
        if value:
            self.answer.setText("<b> yes <\b>")
        else:
            self.answer.setText("<b> no <\b>")

    def connect(self, task):
        self.check_box.clicked.connect(task)


class TextOption(GeneralOption):
    def __init__(self, option, value):
        super().__init__(option)
        self.text = QW.QLineEdit(value)

        self.layout.addWidget(self.text)

    def get_value(self):
        if not self.text.text():
            return "''"
        return self.text.text()

    def set_value(self, text):
        self.text.setText(text)

    def set_enable(self, switch_bool):
        self.text.setEnabled(switch_bool)


class NumberOption(GeneralOption):
    def __init__(self, option, value):
        super().__init__(option)

        if isinstance(value, float):
            self.number = QW.QDoubleSpinBox()
        else:
            self.number = QW.QSpinBox()
        self.number.setRange(-1000000, 1000000)
        self.number.setValue(value)
        self.number.installEventFilter(self)

        self.layout.addWidget(self.number)

    def get_value(self):
        return self.number.value()

    def set_value(self, value):
        self.number.setValue(value)

    def set_range(self, min_value, max_value):
        self.number.setRange(min_value, max_value)

    def set_enable(self, switch_bool):
        self.number.setEnabled(switch_bool)

    def eventFilter(self, source, event):
        if event.type() == QtCore.QEvent.Wheel and source is self.number:
            return True
        return False


class FileBrowserOption(GeneralOption):
    def __init__(self, option, value=None, directory=False):
        super().__init__(option)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.browser = FileBrowser(directory=directory)
        self.layout.addWidget(self.browser)
        if value:
            self.browser.set_filenames(value)

    def get_value(self):
        if not self.browser.get_filenames():
            return "''"
        return self.browser.get_filenames()

    def set_enable(self, switch_bool):
        self.browser.set_enable(switch_bool)


class ComboBoxOption(GeneralOption):
    def __init__(self, option):
        super().__init__(option)
        self.combobox = QW.QComboBox()
        self.layout.addWidget(self.combobox)

    def add_items(self, items_list):
        """ Fill the combobox with a list of string
        """
        self.combobox.clear()
        self.combobox.addItems([None] + items_list)

    def get_value(self):
        return self.combobox.currentText()


class SVGDialog(QW.QDialog):
    def __init__(self, filename):
        super().__init__()
        self.main_layout = QW.QVBoxLayout(self)
        self.setWindowTitle("DAG")

        if os.path.exists(filename):
            widget = QSvgWidget(filename)
            self.main_layout.addWidget(widget)


class SnakemakeOptionDialog(QW.QDialog):
    """ Widget to set up options of snakemake and launch pipeline. It provides
    a progress bar to know how your jobs work.
    """
    def __init__(self, parent=None):
        super().__init__(parent=parent)
        self.main_layout = QW.QVBoxLayout(self)
        self.setWindowTitle("Snakemake options")

        # Settings option
        settings = QtCore.QSettings("Sequana", "Snakemake_option")
        self._cores = settings.value("cores", 2, type=int)
        self._jobs = settings.value("jobs", 2, type=int)
        self._cluster = settings.value("cluster", "", type=str)
        self._tab_pos = settings.value("tab_pos", 0, type=int)

        self.set_launch_options()

        footer = self.footer_button()
        self.main_layout.addWidget(self.tabs)
        self.main_layout.addWidget(footer)

    def ok_event(self):
        settings = QtCore.QSettings("Sequana", "Snakemake_option")
        settings.setValue("cores", self.cores_option.get_value())
        settings.setValue("jobs", self.jobs_option.get_value())
        settings.setValue("cluster", self.cluster_options.get_value())
        settings.setValue("tab_pos", self.tabs.currentIndex())
        self.close()

    def cancel_event(self):
        self.cores_option.set_value(self._cores)
        self.jobs_option.set_value(self._jobs)
        self.cluster_options.set_value(self._cluster)
        self.tabs.setCurrentIndex(self._tab_pos)
        self.close()

    def footer_button(self):
        ok_btn = QW.QPushButton("Ok")
        ok_btn.clicked.connect(self.ok_event)
        cancel_btn = QW.QPushButton("Cancel")
        cancel_btn.clicked.connect(self.cancel_event)

        footer_widget = QW.QWidget()
        footer_layout = QW.QHBoxLayout(footer_widget)
        footer_layout.addWidget(ok_btn)
        footer_layout.addWidget(cancel_btn)

        return footer_widget

    def set_launch_options(self):
        self.cores_option = NumberOption("--cores", self._cores)
        self.jobs_option = NumberOption("--jobs", self._jobs)
        self.cluster_options = TextOption("--cluster", self._cluster)

        launch_cluster_widget = QW.QWidget()
        launch_cluster_layout = QW.QVBoxLayout(launch_cluster_widget)
        launch_cluster_layout.addWidget(self.cluster_options)
        launch_cluster_layout.addWidget(self.jobs_option)

        general_options_widget = QW.QWidget()
        general_options_layout = QW.QVBoxLayout(general_options_widget)

        self.jobs_option.set_range(1, 10000)
        self.cores_option.set_range(1, multiprocessing.cpu_count())

        # create tab
        self.tabs = QW.QTabWidget()
        self.tabs.addTab(self.cores_option, "Local")
        self.tabs.addTab(launch_cluster_widget, "Cluster")
        self.tabs.addTab(general_options_widget, "General")
        self.tabs.setCurrentIndex(self._tab_pos)

    def get_snakemake_options(self):
        # cluster/local option
        current_tab = self.tabs.currentWidget()
        try:
            option_list = [current_tab.get_name(),
                           str(current_tab.get_value())]
        except AttributeError:
            current_layout = current_tab.layout()
            widgets = (current_layout.itemAt(i).widget() for i in
                       range(current_layout.count()))
            option_list = [str(x) for w in widgets for x in w.get_tuple()]

        # drop options with no arguments !
        new_options = []
        for i, this in enumerate(option_list):
            if this is not None and this not in ["", '', "''", '""']:
                new_options.append(this)
            else:
                # remove the options that has no value associated
                # for instance --cluster '' should be removed
                _ = new_options.pop()
        return new_options


class CriticalMessage(QW.QMessageBox):
    def __init__(self, msg, details="", parent=None):
        super().__init__(parent=parent)
        self.setWindowTitle("Error message")
        self.setIcon(QW.QMessageBox.Critical)

        # Force a minimum width ! Cannot use setFixedWidth. This is a trick
        # found on http://www.qtcentre.org/threads/22298-QMessageBox-Controlling-the-width
        layout = self.layout()
        spacer = QW.QSpacerItem(600,0)
        layout.addItem(spacer, layout.rowCount(), 0,1,layout.columnCount())

        msg = '<b style="color:red">' + msg + "</b><br><br>"
        try: details = str(details).replace("\\n", "<br>")
        except: pass
        self.setText(msg + details)


class WarningMessage(QW.QMessageBox):
    def __init__(self, msg, parent=None):
        super().__init__(parent=parent)
        self.setWindowTitle("Warning message")
        self.setIcon(QW.QMessageBox.Warning)
        self.setText(msg)


class InfoMessage(QW.QMessageBox):
    def __init__(self, msg, parent=None):
        super().__init__(parent=parent)
        self.setWindowTitle("Info")
        self.setIcon(QW.QMessageBox.Information)
        self.setText(msg)


class DirectoryDialog(QW.QFileDialog):
    def __init__(self, parent, title, directory, file_filter):
        super().__init__(parent)
        self.setAcceptMode(QW.QFileDialog.AcceptOpen)
        self.setFileMode(QW.QFileDialog.Directory)
        self.setViewMode(QW.QFileDialog.Detail)
        self.setWindowTitle(title)
        self.setDirectory(directory)
        self.setNameFilter(file_filter)

    def get_directory_path(self):
        if self.exec_():
            return self.selectedFiles()[0]
        return None


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
        for i in range(0, 100):
            #progressBar.setValue(i)
            t = time.time()
            while time.time() < t + 2./100.:
                app.processEvents()

        app.processEvents()
        sequana = SequanaGUI()
        sequana.show()
        splash.finish(sequana)
        sys.exit(app.exec_())


if __name__ == "__main__":
    main()
