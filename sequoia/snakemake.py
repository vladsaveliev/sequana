"""



"""
import os
from os.path import isdir
from easydev import get_package_location as gpl





class RuleBase(object):
    def __init__(self):
        self.basedir = gpl("sequoia") + os.sep + "pipelines"


class Rules(RuleBase):
    def __init__(self):
        super(Rules, self).__init__()
        self.names = [this for this in os.listdir(self.basedir)
            if isdir(self.basedir + os.sep + this)]

    def isvalid(self, name):
        if name in self.names:
            return True
        else:
            return False


class Rule(RuleBase):
    """

    Rules provides a simple way to retrieve the path of a Snakefile
    for a given rule. Snakefiles are stored in sequoia/pipelines.
    For instance, in the following example, we wish to known the path
    of the Snakefile to ne found in sequoia/pipelines/dag

    ::

        from sequoia import Rules
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
            msg = "The rule %s is not part of the sequoia workflows"
            raise ValueError(msg  % name)

        self.name = name

        self.location = self.basedir + os.sep + self.name

        if os.path.exists(self.location + os.sep + "Snakefile"):
            self.location = self.location + os.sep + "Snakefile"
        elif os.path.exists(self.location + os.sep + "Snakefile." + self.name):
            self.location = self.location + os.sep + "Snakefile." + self.name
        else:
            print("not found")


        self.description = None
        try:
            with open(self.location + os.sep + "README.rst", "r") as fh:
                self.description = fh.read()
        except:
            self.description = "no description" 

    def __str__(self):
        txt = "Rule **" + self.name + "**:\n" + self.description
        return txt



rules = {}
for name in Rules().names:
    rules[name] = Rule(name).location 
