"""



"""
import os
import sys
import json

from os.path import isdir
from easydev import get_package_location as gpl

import pandas as pd
import pylab


class SnakeMakeProfile(object):
    def __init__(self, filename):
        self.filename = filename
        
    def parse(self):
        data = json.loads(self.filename)


class SnakeMakeStats(object):
    def __init__(self, filename):
        self.filename = filename
        
    def parse_data(self):
        with open(self.filename, 'r') as fin:
            data = json.load(fin)
        return data

    def plot(self, fontsize=16):
        df = pd.DataFrame(self.parse_data()['rules'])
        ts = df.ix['mean-runtime']
        ts.plot.barh(fontsize=fontsize)
        pylab.grid(True)
        pylab.xlabel("Seconds (s)", fontsize=fontsize)
        try:pylab.tight_layout()
        except:pass
        

class RuleBase(object):
    def __init__(self):
        self.basedir = gpl("sequana") + os.sep + "pipelines"


class Rules(RuleBase):
    def __init__(self):
        super(Rules, self).__init__()
        self.names = [this for this in os.listdir(self.basedir)
            if isdir(self.basedir + os.sep + this)]

        for this in ["__pycache__"]:
            try:self.names.remove(this)
            except:pass



    def isvalid(self, name):
        if name in self.names:
            return True
        else:
            return False


class Rule(RuleBase):
    """

    Rules provides a simple way to retrieve the path of a Snakefile
    for a given rule. Snakefiles are stored in sequana/pipelines.
    For instance, in the following example, we wish to known the path
    of the Snakefile to ne found in sequana/pipelines/dag

    ::

        from sequana import Rules
        filename = Rules('dag').filename

    this returns the full path of the Snakefile.

    """
    def __init__(self, name):
        """

        :param str snakefile: name of a registered rule

        """
        super(Rule, self).__init__()
        self._rules = Rules()

        if self._rules.isvalid(name) is False:
            msg = "The rule %s is not part of the sequana workflows"
            raise ValueError(msg  % name)

        self.name = name

        self.location = self.basedir + os.sep + self.name

        if os.path.exists(self.location + os.sep + "Snakefile"):
            self.location = self.location + os.sep + "Snakefile"
        elif os.path.exists(self.location + os.sep + "Snakefile." + self.name):
            self.location = self.location + os.sep + "Snakefile." + self.name
        else:
            print("Snakefile for %s not found" % self.name)


        self.description = None
        try:
            with open(self.location + os.sep + "README.rst", "r") as fh:
                self.description = fh.read()
        except:
            self.description = "no description" 

    def __str__(self):
        txt = "Rule **" + self.name + "**:\n" + self.description
        return txt


#: define a dictionary to be used in Snakefile to include sequana's snakefiles
rules = {}
for name in Rules().names:
    rules[name] = Rule(name).location


class ValidateConfig(object):
    """

    Converts json or yaml into a dictionary with keys accessible as attributes
    With config.json file content as::

        {'e':1}

    type::

        >>> vc = ValidateConfig(config)
        >>> config = vc()
        >>> config.e == 1
        True

    """
    def __init__(self, filename): 
        """Could be a json or a yaml"""
        self.filename = filename

        if isinstance(filename, str):
            if filename.endswith('json'):
                self.config = json.load(open(self.filename, 'r'))
            else:
                raise NotImplementedError
        else:
            self.config = filename 

    def __call__(self):
        from easydev import AttrDict
        config = AttrDict(**self.config)
        return config


def message(mes):
    from easydev.console import purple
    print("// -- " + purple(mes))



class DOTParser(object):
    def __init__(self, filename):
        self.filename = filename

    def add_urls(self):
        with open(self.filename, "r") as fh:
            data = fh.read()

        with open(self.filename.replace(".dot", ".ann.dot"), "w") as fout:
            for line in data.split("\n"):
                if "[label =" not in line:
                    fout.write(line + "\n")
                else:
                    separator = "color ="
                    lhs, rhs = line.split(separator)
                    name = lhs.split("label =")[1]
                    name = name.replace(",","")
                    name = name.replace('"',"")
                    name = name.strip()
                    line = lhs + ' URL="%s.html" target="_blank", ' % name 
                    line += separator + rhs
                    fout.write(line + "\n")

    #  label="cutadapt.html", URL="cutadapt.html", target="_blank",





