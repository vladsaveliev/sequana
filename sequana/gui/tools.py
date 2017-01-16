# -*- coding: utf-8 -*-
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
import colorlog
import shutil

from PyQt5 import QtWidgets as QW


__all__ = ['Logger', 'Tools', 'QPlainTextEditLogger']

class Logger(object):
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


class Tools(Logger):
    def copy(self, source, target):
        try:
            shutil.copy(source, target)
        except Exception as err:
            self.error(err)
            self.warning("Cannot overwrite existing file. (Probably identical)")


class QPlainTextEditLogger(colorlog.StreamHandler):
    def __init__(self, parent):
        super().__init__()
        self.widget = QW.QPlainTextEdit(parent)
        self.widget.setReadOnly(True)
        self.bgcolor = "#aabbcc"
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

