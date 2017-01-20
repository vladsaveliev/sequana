from PyQt5 import QtWidgets as QW
from sequana.gui.messages import WarningMessage


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
        self.btn.setFixedSize(100, 20)

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
        widget_layout.setContentsMargins(3, 3, 3, 3)
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
        file_path = QW.QFileDialog.getOpenFileNames(self, 
            "Select a sample", ".", self.filter)[0]

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
            self.set_green()

    def browse_directory(self):
        self.dialog = DirectoryDialog(self, "Select a directory", ".", 
                                 self.filter)
        directory_path = self.dialog.get_directory_path()
        if directory_path:
            self.set_filenames(directory_path)
        else:
            self.set_empty_path()

    def browse_file(self):
        try:
            file_path = QW.QFileDialog.getOpenFileNames(self, 
                "Single File", ".",  self.filter)[0][0]
            self.set_filenames(file_path)
        except IndexError:
            self.set_empty_path()

    def get_filenames(self):
        return self.paths

    def set_filenames(self, filename):
        self.paths = filename
        if len(filename) > 30:
            self.btn_filename.setText("...." + filename[-30:])
            self.btn_filename.setToolTip(filename)
        else:
            self.btn_filename.setText(filename)
        self.set_green()

    def set_empty_path(self):
        self.btn_filename.setText(self.empty_msg)
        self.paths = ""
        self.set_red()

    def set_enable(self, switch_bool):
        if switch_bool:
            self.setup_color()
        else:
            self.btn.setStyleSheet("QPushButton {background-color: #AAAAAA; "
                                   "color: #222222}")
        self.btn.setEnabled(switch_bool)

    def path_is_setup(self):
        return self.setup

    def set_red(self):
        self.setup = False
        self.setup_color()

    def set_green(self):
        self.setup = True
        self.setup_color()

    def clicked_connect(self, function):
        """ Connect additionnal function on browser button. It is used to
        activate run button in Sequana GUI.
        """
        self.btn.clicked.connect(function)

    def changed_connect(self, function):
        """ Trigger function when user change the location.
        """
        self.btn_filename.textChanged.connect(function)


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
        if self.isVisible() or self.exec_():
            return self.selectedFiles()[0]
        return None
