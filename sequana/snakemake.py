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
    """

    Run the Snakemake with this option::

        -- stats stats.txt

    Then:

    .. plot::
        :include-source:

        from sequana.snakemake import SnakeMakeStats
        from sequana import sequana_data
        filename = sequana_data("test_snakemake_stats.txt")
        s = SnakeMakeStats(filename)
        s.plot()

    """
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
        self.basedir = gpl("sequana") + os.sep.join(["sequana", "rules"])


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
    """Utility to parse the dot returned by Snakemake and add URLs automatically


    """
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


class Modules(object):
    """Class to get information about a module/pipeline

    ::

        from sequana.snakemake import Modules
        m = Modules()
        m.onweb('dag')
        m.info('dag')


    .. todo::

    """
    def __init__(self):
        self.rules = {}
        for name in Rules().names:
            self.rules[name] = Rule(name).location
        self.registered = rules.keys()

    def onweb(self, name):
        """Open web page with the README file corresponding to the module

        """
        assert name in self.names
        from easydev import onweb
        url = "https://github.com/sequana/sequana/blob/master/rules/" 
        url += name + "/README.rst"
        onweb(url)

    def info(self, name):
        """print the README of the module"""
        assert name in self.names
        filename = self.rules[name]
        lhs, rhs = filename.rsplit("/", 1)
        filename = lhs + "/README.rst"
        with open(filename, "r") as fin:
            print(fin.read())

    def _get_names(self):
        return self.registered
    names = property(_get_names)


modules = Modules()



import glob


class GetInOutFiles(object):
    """
    Given list of input files, replaces the input directory with output one:


    Given this file : **input/test.csv**::

        >>> inout = GetInOutFiles("input/*csv", "output")
        >>> ins, outs, wkdir = inout.getinfo()
        >>> ins
        ["input/test.csv"]
        >>> outs
        ["output/test.csv"]
        >>> wdir
        "output"

    """
    def __init__(self, wildcard, wkdir):
        self.wildcard = wildcard
        self.wkdir = wkdir

        # store the input_filenames
        self.input_filenames = glob.glob(wildcard)

        # create the output filenames
        self.output_filenames = []
        for input_filename in self.input_filenames:
            path, filename = os.path.split(input_filename)
            output_filename = wkdir + os.sep + filename
            self.output_filenames.append(output_filename)
            

def get_inout(inputs, wkdir):
    s = GetInOutFiles(inputs, wkdir)
    return s.input_filenames, s.output_filenames, s.wkdir 


def get_filenames(wildcard):
    import os
    import glob
    filenames = glob.glob(wildcard)
    filenames = [os.path.split(filename)[1] for filename in filenames]
    return filenames


def get_cleanup_rules(filename):
    import snakemake
    s = snakemake.Workflow(filename)
    s.include(filename)
    names = [rule.name for rule in list(s.rules) if rule.name.endswith('_cleanup')]
    return names



