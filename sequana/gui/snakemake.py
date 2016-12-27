# coding: utf-8
#
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
"""Snakemake Dialog for the main GUI application"""
from sequana.gui.ui_snakemake import Ui_Snakemake

from PyQt5 import QtWidgets as QW
from PyQt5 import QtCore


class SnakemakeDialog(QW.QDialog):
    """ Widget to set up options of snakemake and launch pipeline. It provides
    a progress bar to know how your jobs work.
    """
    def __init__(self, parent=None):
        super().__init__(parent=parent)
        self.ui = Ui_Snakemake()
        self.ui.setupUi(self)

        settings = QtCore.QSettings("Sequana", "snakemake_options")
        #self._cores = settings.value("cores", 2, type=int)
        #self._jobs = settings.value("jobs", 2, type=int)
        #self._cluster = settings.value("cluster", "", type=str)
        self._tab_pos = settings.value("tab_pos", 0, type=int)

        print(self._tab_pos)
        self.ui.tabs.setCurrentIndex(self._tab_pos)
        #footer = self.footer_button()

    def ok_event(self):
        settings = QtCore.QSettings("Sequana", "snakemake_options")
        items = self.get_items()
        for k,v in self.get_items().items():
            print(k,v)
            settings.setValue(k, v)
        settings.setValue("tab_pos", self.ui.tabs.currentIndex())
        self.close()

    def cancel_event(self):
        #self.cores_option.set_value(self._cores)
        #self.jobs_option.set_value(self._jobs)
        #self.cluster_options.set_value(self._cluster)
        self.ui.tabs.setCurrentIndex(self._tab_pos)
        self.close()

    def _footer_button(self):
        ok_btn = QW.QPushButton("Ok")
        ok_btn.clicked.connect(self.ok_event)
        cancel_btn = QW.QPushButton("Cancel")
        cancel_btn.clicked.connect(self.cancel_event)

        footer_widget = QW.QWidget()
        footer_layout = QW.QHBoxLayout(footer_widget)
        footer_layout.addWidget(ok_btn)
        footer_layout.addWidget(cancel_btn)

        return footer_widget

    def _set_launch_options(self):
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

    def get_items(self):
        # get all items to save in settings
        items = {}
        names = [this for this in dir(self.ui) if this.startswith("snakemake_options")]
        names = [this for this in widgets if this.endswith('_value')]
        values = [getattr(gui.snakemake_dialog.ui, x).text() for x in values]

        for name, val in zip(names, values):
            items[name] = val
        return items

    def set_items(self):
        # set all items 
        pass

    def get_widgets(self, prefix):
        # identify names of the widget objects
        widgets = [this for this in dir(self.ui) if this.startswith(prefix)]
        widgets = [this for this in widgets if this.endswith('_value')]

        # Get the objects themselves
        names = [this for this in widgets]
        widgets = [getattr(self.ui, this) for this in widgets]

        names = [this.replace(prefix + "_", "") for this in names]
        names = [this.rstrip("_value") for this in names]

        options = []
        for name, widget in zip(names, widgets):
            options.append(SOptions(name, widget))
        return options

    def _get_options(self, prefix):
        options = []
        for widget in self.get_widgets(prefix):
            option = widget.get_option()
            if option:
                options.extend(option)
        return options

    def get_snakemake_local_options(self):
        return self._get_options("snakemake_options_local")

    def get_snakemake_cluster_options(self):
        return self._get_options("snakemake_options_cluster")

    def get_snakemake_general_options(self):
        return self._get_options("snakemake_options_general")


class SOptions(object):

    def __init__(self, name, widget):
        self.name = name
        self.widget = widget

    def _get_value(self):
        """Return a string"""
        if isinstance(self.widget, QW.QSpinBox):
            value = self.widget.text()
        elif isinstance(self.widget, QW.QLabel):
            value = self.widget.text()
        else:
            try:
                value = self.widget.text()
            except:
                print("unknown widget" + str(type(self.widget)))
                value = ""
        return value

    def get_option(self):
        """Return option and its value as a list

        An option may be without any value, but we still return a list

        e.G. ["--verbose"] or ["--cores", "2"]
        """
        if isinstance(self.widget, QW.QCheckBox):
            if self.widget.isChecked() is True:
                return ["--" + self.name]
        else:
            value = self._get_value()
            if value is None or value in ["", '', "''", '""']:
                return []
            else:
                return ["--" + self.name, value]

