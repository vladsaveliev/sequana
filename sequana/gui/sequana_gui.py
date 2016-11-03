# coding: utf-8

import os
import sys
import shutil
import subprocess as sp
import multiprocessing

from PyQt5.QtCore import Qt, QCoreApplication, pyqtSlot, QEvent
from PyQt5.QtWidgets import (QApplication, QPushButton, QComboBox, QWidget,
                             QLineEdit, QTabWidget, QLabel, QFrame, QGroupBox,
                             QCheckBox, QSpinBox, QDoubleSpinBox, QScrollArea,
                             QMenuBar, QAction)
from PyQt5.QtWidgets import QFileDialog, QDialog
from PyQt5.QtWidgets import QFormLayout, QHBoxLayout, QVBoxLayout

from sequana import snaketools
from sequana.snaketools import Module


class SequanaGUI(QWidget):
    """ Sequana GUI !
    """

    _not_a_rule = {"requirements", "gatk_bin", "input_directory",
                   "pattern", "samples"}

    def __init__(self):
        super().__init__()

        # snakemake dialog windows
        self.snakemake_dialog = SnakemakeOptionDialog()

        # set menu bar
        self.menu_bar = QMenuBar(self)
        options_menu = self.menu_bar.addMenu("File")


        # box to choose the pipeline
        self.choice_flag = False
        choice_layout = QHBoxLayout()
        self.create_choice_button()
        choice_layout.addWidget(self.choice_button)

        # create browser tab to choice directory or files
        self.create_tabs_browser()
        self.tabs_browser.currentChanged.connect(self.switch_run)

        # select the working directory
        groupbox_layout = QHBoxLayout()
        self.working_dir = FileBrowser(directory=True)
        self.working_dir.clicked_connect(self.switch_run)
        groupbox_layout.addWidget(self.working_dir)
        groupbox = QGroupBox("Working directory")
        groupbox.setLayout(groupbox_layout)

        # formular which contains all options the pipeline chosen
        widget_formular = QWidget()
        self.formular = QVBoxLayout(widget_formular)
        self.formular.setSpacing(0)
        scroll_area = QScrollArea()
        scroll_area.setWidget(widget_formular)
        scroll_area.setWidgetResizable(True)
        scroll_area.setMinimumHeight(250)

        # main layout
        vlayout = QVBoxLayout(self)
        vlayout.insertSpacing(0, 10)

        # add widgets in layout
        vlayout.addLayout(choice_layout)
        vlayout.addWidget(self.tabs_browser)
        vlayout.addWidget(groupbox)
        vlayout.addWidget(scroll_area)
        vlayout.addWidget(self.create_footer_button())

        # main window options
        self.setGeometry(200, 200, 500, 500)
        self.setWindowTitle("Sequana")
        self.show()

    def clearLayout(self, layout):
        """ Clean all widget contained in a layout.
        """
        while layout.count():
            child = layout.takeAt(0)
            if child.widget() is not None:
                child.widget().deleteLater()
            elif child.layout() is not None:
                self.clearLayout(child.layout())

    def create_choice_button(self):
        """ Create button to select the wished pipeline.
        """
        self.choice_button = QComboBox()
        snaketools.pipeline_names.sort()
        self.choice_button.addItems(["select pipeline"] +
                                    snaketools.pipeline_names)
        self.choice_button.currentIndexChanged[str].connect(
            self.on_pipeline_choice)
        self.choice_button.installEventFilter(self)

    @pyqtSlot(str)
    def on_pipeline_choice(self, index):
        """ Change options formular when user change the pipeline.
        """
        config_file = snaketools.Module(index)._get_config()
        self.choice_button.removeItem(
            self.choice_button.findText("select pipeline"))
        self.create_base_formular(config_file)
        self.choice_flag = True
        self.switch_run()

    def create_base_formular(self, config_file):
        """ Create formular with all options necessary for a pipeline.
        """
        self.clearLayout(self.formular)
        config_dict = snaketools.SequanaConfig(config_file).config
        rules_list = list(config_dict.keys())
        rules_list.sort()
        self.necessary_dict = {}
        for rule in rules_list:
            # Check if dictionnary or not
            contains = config_dict[rule]
            if isinstance(contains, dict) and (
                    rule not in SequanaGUI._not_a_rule):
                rule_box = RuleFormular(rule, contains)
                self.formular.addWidget(rule_box)
            else:
                self.necessary_dict = dict(self.necessary_dict,
                                           **{rule: contains})

    def create_tabs_browser(self):
        """ Generate file browser widget.
        """
        # create browser file
        paired_tab = FileBrowser(paired=True)
        directory_tab = FileBrowser(directory=True)
        paired_tab.clicked_connect(self.switch_run)
        directory_tab.clicked_connect(self.switch_run)
        # create tab box
        self.tabs_browser = QTabWidget()
        self.tabs_browser.addTab(paired_tab, "Paired end")
        self.tabs_browser.addTab(directory_tab, "Directory")

    def create_footer_button(self):
        """ Create Run/Save/Quit buttons
        """
        self.run_btn = QPushButton("Run")
        self.run_btn.setEnabled(False)
        self.run_btn.clicked.connect(self.run_sequana)

        save_btn = QPushButton("Save")
        save_btn.clicked.connect(self.save_config_file)

        quit_btn = QPushButton("Quit")
        quit_btn.clicked.connect(self.close)

        footer_widget = QWidget()
        footer_layout = QHBoxLayout(footer_widget)
        footer_layout.addWidget(self.run_btn)
        footer_layout.addWidget(save_btn)
        footer_layout.addWidget(quit_btn)
        return footer_widget

################

    def run_sequana(self):
        working_dir = self.working_dir.get_filenames()
        self.save_config_file()
        rules = self.get_rules(self.formular)
        module = Module(self.choice_button.currentText())
        self.snakemake_dialog.fill_options(rules)
        self.snakemake_dialog.exec_()

###############

    def switch_run(self):
        if self.working_dir.path_is_setup():
            if self.tabs_browser.currentWidget().path_is_setup():
                if self.working_dir.path_is_setup():
                    if self.choice_flag:
                        return self.run_btn.setEnabled(True)
        return self.run_btn.setEnabled(False)

    def save_config_file(self):
        formular_dict = dict(self.create_formular_dict(self.formular),
                             **self.necessary_dict)
        # get samples names or input_directory
        if self.tabs_browser.currentIndex() == 0:
            formular_dict["samples"] = (
                self.tabs_browser.currentWidget().get_filenames())
        else:
            formular_dict["input_directory"] = (
                self.tabs_browser.currentWidget().get_filenames())
        yaml = snaketools.SequanaConfig(data=formular_dict,
                                        test_requirements=False)
        if self.working_dir.path_is_setup():
            yaml.save(self.working_dir.get_filenames() + "/config.yaml")
        else:
            yaml.save()

    def create_formular_dict(self, layout):
        widgets = (layout.itemAt(i).widget() for i in range(layout.count()))
        formular_dict = {w.get_name(): w.get_value() if w.is_option()
                         else self.create_formular_dict(w.get_layout())
                         for w in widgets}
        return formular_dict

    def get_rules(self, layout):
        widgets = (layout.itemAt(i).widget() for i in range(layout.count()))
        rules = [w.get_name() for w in widgets]
        return rules

    def eventFilter(self, source, event):
        """ Inactivate wheel event of combobox
        """
        if event.type() == QEvent.Wheel and source is self.choice_button:
            return True
        return False


class FileBrowser(QWidget):
    """ Class to create a file browser in PyQT5.
    """
    def __init__(self, paired=False, directory=False):
        super().__init__()

        self.paths = ""
        self.setup = False
        self.btn = QPushButton("Browse")
        self.btn_filename = QLabel("No file selected")
        if directory:
            self.btn_filename.setText("No directory selected")
            self.btn.clicked.connect(self.browse_directory)
        elif paired:
            self.btn.clicked.connect(self.browse_paired_file)
        else:
            self.btn.clicked.connect(self.browse_file)

        widget_layout = QHBoxLayout(self)
        widget_layout.addWidget(self.btn)
        widget_layout.addWidget(self.btn_filename)

    def browse_paired_file(self):
        file_path = QFileDialog.getOpenFileNames(self,
                                                 "Paired end File", ".")[0]
        self.setup = False
        if not file_path:
            self.btn_filename.setText("No file selected")
            self.paths = {"file1": "", "file2": ""}
        elif len(file_path) > 2:
            self.btn_filename.setText("No file selected")
        else:
            self.paths = {"file{0}".format(i+1): file_path[i]
                          for i in range(0, len(file_path))}
            self.btn_filename.setText("\n".join([key + ": " + value
                                      for key, value in self.paths.items()]))
            self.setup = True
        
    def browse_directory(self):
        directory_path = QFileDialog.getExistingDirectory(self,
                                                          "Directory", ".")
        try:
            self.btn_filename.setText(directory_path)
            self.paths = directory_path
            self.setup = True
            if not directory_path:
                self.setup = False
        except IndexError:
            self.btn_filename.setText("No directory selected")
            self.paths = ""
            self.setup = False

    def browse_file(self):
        try:
            file_path = QFileDialog.getOpenFileNames(self,
                                                     "Single File", ".")[0][0]
            self.btn_filename.setText(file_path)
            self.paths = file_path
            self.setup = True
        except IndexError:
            self.btn_filename.setText("No file selected")
            self.paths = ""
            self.setup = False

    def get_filenames(self):
        return self.paths

    def path_is_setup(self):
        return self.setup

    def clicked_connect(self, function):
        """ Connect additionnal function on browser button. It is used to
        activate run button in Sequana GUI.
        """
        self.btn.clicked.connect(function)


class RuleFormular(QGroupBox):
    def __init__(self, rule_name, rule_dict):
        super().__init__(rule_name)

        self.rule_name = rule_name
        self.rule_dict = rule_dict
        self.layout = QVBoxLayout(self)

        for option, value in self.rule_dict.items():
            if option == "reference":
                option = FileBrowserOption(option)
            elif isinstance(value, bool):
                option = BooleanOption(option, value)
            else:
                try:
                    option = NumberOption(option, value)
                except TypeError:
                    try:
                        option = TextOption(option, value)
                    except TypeError:
                        option = RuleFormular(option, value)
            self.layout.addWidget(option)

    def get_name(self):
        return self.rule_name

    def get_layout(self):
        return self.layout

    def is_option(self):
        return False


class GeneralOption(QWidget):
    """ Parent class for Options. It define design of options
    """
    def __init__(self, option):
        super().__init__()
        self.option = option

        self.layout = QHBoxLayout(self)
        self.layout.addWidget(QLabel(option))

    def get_name(self):
        return self.option

    def is_option(self):
        return True

    def get_value(self):
        pass


class BooleanOption(GeneralOption):
    """ Wrapp QCheckBox class
    """
    def __init__(self, option, value):
        super().__init__(option)

        self.check_box = QCheckBox()
        self.check_box.setChecked(value)

        self.layout.addWidget(self.check_box)

    def get_value(self):
        return self.check_box.isChecked()


class TextOption(GeneralOption):
    def __init__(self, option, value):
        super().__init__(option)
        self.text = QLineEdit(value)

        self.layout.addWidget(self.text)

    def get_value(self):
        if not self.text.text():
            return " "
        return self.text.text()


class NumberOption(GeneralOption):
    def __init__(self, option, value):
        super().__init__(option)

        if isinstance(value, float):
            self.number = QDoubleSpinBox()
        else:
            self.number = QSpinBox()
        self.number.setRange(-1000000, 1000000)
        self.number.setValue(value)
        self.number.installEventFilter(self)

        self.layout.addWidget(self.number)

    def get_value(self):
        return self.number.value()

    def set_range(self, min_value, max_value):
        self.number.setRange(min_value, max_value)

    def eventFilter(self, source, event):
        if event.type() == QEvent.Wheel and source is self.number:
            return True
        return False


class FileBrowserOption(GeneralOption):
    def __init__(self, option):
        super().__init__(option)
        self.browser = FileBrowser()
        self.layout.addWidget(self.browser)

    def get_value(self):
        if not self.browser.get_filenames():
            return " "
        return self.browser.get_filenames()


class ComboBoxOption(GeneralOption):
    def __init__(self, option):
        super().__init__(option)
        self.combobox = QComboBox()
        self.layout.addWidget(self.combobox)

    def add_items(self, items_list):
        self.combobox.clear()
        self.combobox.addItems([None] + items_list)

    def get_value(self):
        return self.combobox.currentText()


class SnakemakeOptionDialog(QDialog):
    """ Widget to set up options of snakemake and launch pipeline. It provides
    a progress bar to know how your jobs work.
    """
    def __init__(self):
        super().__init__()
        self.main_layout = QVBoxLayout(self)
        self.setWindowTitle("Snakemake options")
        self.cores_option = NumberOption("--cores", 2)
        self.forcerun_option = ComboBoxOption("--forcerun")
        self.until_option = ComboBoxOption("--until")

        self.cores_option.set_range(1, multiprocessing.cpu_count())

        self.main_layout.addWidget(self.cores_option)
        self.main_layout.addWidget(self.forcerun_option)
        self.main_layout.addWidget(self.until_option)

    def fill_options(self, rules):
        self.forcerun_option.add_items(rules)
        self.until_option.add_items(rules)

    def cluster_options(self):
        self.cluster_check_box = BooleanOption("On a cluster ?", False)



if __name__ == "__main__":
    app = QApplication(sys.argv)
    sequana = SequanaGUI()
    sys.exit(app.exec_())
