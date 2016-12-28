# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'preferences.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Preferences(object):
    def setupUi(self, Preferences):
        Preferences.setObjectName("Preferences")
        Preferences.resize(400, 300)
        self.buttonBox = QtWidgets.QDialogButtonBox(Preferences)
        self.buttonBox.setGeometry(QtCore.QRect(30, 240, 341, 32))
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.verticalLayoutWidget = QtWidgets.QWidget(Preferences)
        self.verticalLayoutWidget.setGeometry(QtCore.QRect(19, 19, 294, 95))
        self.verticalLayoutWidget.setObjectName("verticalLayoutWidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.verticalLayoutWidget)
        self.verticalLayout.setSizeConstraint(QtWidgets.QLayout.SetDefaultConstraint)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.tabs = QtWidgets.QTabWidget(self.verticalLayoutWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.tabs.sizePolicy().hasHeightForWidth())
        self.tabs.setSizePolicy(sizePolicy)
        self.tabs.setObjectName("tabs")
        self.tab_general = QtWidgets.QWidget()
        self.tab_general.setObjectName("tab_general")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.tab_general)
        self.horizontalLayout.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.formLayout = QtWidgets.QFormLayout()
        self.formLayout.setObjectName("formLayout")
        self.label = QtWidgets.QLabel(self.tab_general)
        self.label.setObjectName("label")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.label)
        self.preferences_options_general_browser_value = QtWidgets.QComboBox(self.tab_general)
        self.preferences_options_general_browser_value.setObjectName("preferences_options_general_browser_value")
        self.preferences_options_general_browser_value.addItem("")
        self.preferences_options_general_browser_value.addItem("")
        self.preferences_options_general_browser_value.addItem("")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.preferences_options_general_browser_value)
        self.horizontalLayout.addLayout(self.formLayout)
        self.tabs.addTab(self.tab_general, "")
        self.verticalLayout.addWidget(self.tabs, 0, QtCore.Qt.AlignHCenter)

        self.retranslateUi(Preferences)
        self.tabs.setCurrentIndex(0)
        self.buttonBox.accepted.connect(Preferences.accept)
        self.buttonBox.rejected.connect(Preferences.reject)
        QtCore.QMetaObject.connectSlotsByName(Preferences)
        Preferences.setTabOrder(self.tabs, self.preferences_options_general_browser_value)

    def retranslateUi(self, Preferences):
        _translate = QtCore.QCoreApplication.translate
        Preferences.setWindowTitle(_translate("Preferences", "Dialog"))
        self.label.setText(_translate("Preferences", "Select the browser to be used"))
        self.preferences_options_general_browser_value.setItemText(0, _translate("Preferences", "pyqt5"))
        self.preferences_options_general_browser_value.setItemText(1, _translate("Preferences", "firefox"))
        self.preferences_options_general_browser_value.setItemText(2, _translate("Preferences", "safari"))
        self.tabs.setTabText(self.tabs.indexOf(self.tab_general), _translate("Preferences", "General"))

