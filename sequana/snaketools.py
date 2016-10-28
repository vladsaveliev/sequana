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
"""Set of tools to manipulate Snakefile and config files

Here is an overview (see details here below)

.. autosummary::
    :nosignatures:

    sequana.snaketools.DOTParser
    sequana.snaketools.Module
    sequana.snaketools.ModuleFinder
    sequana.snaketools.SnakeMakeProfile
    sequana.snaketools.SnakeMakeStats
    sequana.snaketools.SequanaConfig
    sequana.snaketools.get_cleanup_rules
    sequana.snaketools.message
    sequana.snaketools.modules


"""
import os
import json
import glob

from easydev import get_package_location as gpl
from easydev import load_configfile, AttrDict
from easydev.tools import touch
import pandas as pd
import pylab


# __all__ = ["SequanaConfig"]

try:
    # This is for python2.7
    import snakemake
except:
    print("Snakemake must be installed. Available for Python3 only")

    class MockSnakeMake(object):
        def __init__(self):
            pass

        def Workflow(self, filename):
            raise ImportError
    snakemake = MockSnakeMake()



class _ExpandedSnakeFile(object):
    """Read a Snakefile and its dependencies (include) and create single file

    **Motivation**

    a Snakefile may look like::

        from sequana import snaketools as sm
        include: sm.modules['bwa_phix']
        include: sm.modules['fastqc']

    This is nice and compact but we do not see anymore what the Snakefile does.
    This class will recreate the Snakefile without this compactness so that one
    can see the entire structure. The expansion is performed by :meth:`expand`.

    .. warning:: experimental. docstrings should be removed. lines starting
        with **include** should also be removed.

    """
    def __init__(self, snakefile):
        # read the file and include it (Snakemake API)
        # This import is the package not the sequana.snaketools file
        self.snakefile = snakemake.Workflow(snakefile)
        self.snakefile.include(snakefile)

    def expand(self, output_filename="Snakefile.expanded"):
        """The expansion of Snakefile

        :param str filename: name of the output file

        """
        # open a file to save all included rules and workflow in one place
        with open(output_filename, "w") as workflow_expanded:
            # scan all included workflow
            for include in self.snakefile.included:
                with open(include, "r") as wk_toinclude:
                    data = wk_toinclude.read()
                    workflow_expanded.write(data)


class SnakeMakeProfile(object):
    """Interpret the snakemake profile file

    Run the Snakemake with this option::

        --profile profile.json

    """
    def __init__(self, filename):
        self.filename = filename

    def parse(self):
        data = json.loads(self.filename)


class SnakeMakeStats(object):
    """Interpret the snakemake stats file

    Run the Snakemake with this option::

        -- stats stats.txt

    Then:

    .. plot::
        :include-source:

        from sequana.snaketools import SnakeMakeStats
        from sequana import sequana_data
        filename = sequana_data("test_snakemake_stats.txt", "testing")
        s = SnakeMakeStats(filename)
        s.plot()

    """
    def __init__(self, filename, N=1):
        """.. rubric:: Cosntructor"""
        self.filename = filename
        self.N = N

    def _parse_data(self):
        with open(self.filename, 'r') as fin:
            data = json.load(fin)
        return data

    def plot(self, fontsize=16):
        """Create the barplot from the stats file"""
        pylab.clf()
        df = pd.DataFrame(self._parse_data()['rules'])
        ts = df.ix['mean-runtime']
        ts['total'] = self._parse_data()['total_runtime'] / float(self.N)
        ts.sort_values(inplace=True)

        ts.plot.barh(fontsize=fontsize)
        pylab.grid(True)
        pylab.xlabel("Seconds (s)", fontsize=fontsize)
        try:
            pylab.tight_layout()
        except:
            pass

    def plot_and_save(self, filename="snakemake_stats.png",
                      outputdir="report"):
        self.plot()
        pylab.savefig(outputdir + os.sep + filename)


def plot_stats(inputdir=".", outputdir=".",
               filename="snakemake_stats.png", N=1):
    print("Workflow finished. Creating stats image")
    try:
        SnakeMakeStats("%s/stats.txt" % inputdir, N=N).plot_and_save(
            outputdir=outputdir, filename=filename)
    except Exception as err:
        print(err)
        print("INFO: Could not process %s/stats.txt file" % inputdir)


class ModuleFinder(object):
    """Data structure to hold the :class:`Module` names"""
    def __init__(self):
        """.. rubric:: constructor

        :param list extra_paths:

        .. doctest::

            >>> from sequana import ModuleFinder
            >>> modnames = ModuleFinder()
            >>> modnames.isvalid('dag')
            True
            >>> modnames.isvalid('dummy')
            False

        """

        # names for each directory
        self._paths = {}
        self._type = {}

        # scan the official paths
        self._add_names("rules")
        self._add_names("pipelines")

    def _add_names(self, path):
        sepjoin = os.sep.join
        fullpath = sepjoin([gpl("sequana"), "sequana", path])

        fullpaths = self._iglob(fullpath)
        for this in fullpaths:
            whatever, module_name, filename = this.rsplit(os.sep, 2)
            if module_name in self._paths.keys():
                raise ValueError("Found duplicated name %s. "
                                 "Overwrites previous rule " % module_name)
            self._paths[module_name] = whatever + os.sep + module_name
            self._type[module_name] = path[:-1]

    def _iglob(self, path, extension="rules"):
        try:
            from glob import iglob
            matches = list(iglob("%s/**/*.%s" % (path, extension),
                                 recursive=True))
        except:
            # iglob use recursivity with ** only in py3.5 (snakemake)
            import fnmatch
            import os
            matches = []
            for root, dirnames, filenames in os.walk(path):
                for filename in fnmatch.filter(filenames, '*.' + extension):
                    matches.append(os.path.join(root, filename))
        return matches

    def _get_names(self):
        return sorted(list(self._paths.keys()))
    names = property(_get_names, doc="list of existing module names")

    def isvalid(self, name):
        """Check that a name is an existing and valid module"""
        if name not in self.names:
            return False
        return True

    def is_pipeline(self, name):
        return self._type[name] == "pipeline"


class Module(object):
    """Data structure that holds metadata about a **Module**

    A **Module** in sequana's parlance is a directory that contains
    the following files:

        - A **snakemake** file named after the directory with the extension
          **.rules**
        - A **README.rst** file in restructured text format
        - An optional config file in YAML format named config.yaml.
          Although json format is possible, we use YAML throughout
          **sequana** for consistency. If not found, the config.yaml is
          taken from the parent directory. Rules do not have any but pipelines
          do. So if a pipeline does not provide a config.yaml, the one found
          in ./sequana/sequana/pipelines will be used.
        - a **requirements.txt**

    The name of the module is the name of the directory where the files are
    stored. The **Modules** are stored in sequana/rules and sequana/pipelines
    directories. Those main directories may have sub-directories main names
    should not but duplicated

    The Snakefile should be named **<module_name>.rules**.

    Before explaining the different type of **Modules**, let us remind
    that a **rule** in **Snakemake** terminology looks like::

        rule <name>:
            :input: file1
            :output: file2
            :shell: "cp file1 file2"

    However, in **sequana**, the module may contain several rules.

    Here is the terminolgy. Within a module directory, :

        - if the **Snakefile** includes other **Snakefile** then
          we consider that that Module is a **Pipeline** module.
        - if the **Snakefile** does not include other **Snakefile** and
          all internal snakemake rules start with the same name, we consider
          that that module is a **Rule** module.

    This data structure ease the retrieval of metadata and Snakefiles stored
    within Sequana.

    """
    def __init__(self, name):
        self._mf = ModuleFinder()
        self._mf.isvalid(name)
        if name not in self._mf.names:
            raise ValueError("""Sequana error: unknown rule or pipeline. Check
the source code at:

    https://github.com/sequana/sequana/tree/develop/sequana/pipelines and
    https://github.com/sequana/sequana/tree/develop/sequana/rules

or open a Python shell and type::

    import sequana
    sequana.modules.keys()""")
        else:
            self._path = self._mf._paths[name]
        self._name = name

        # could look into ./rules or ./pipelines
        self._snakefile = None
        self._description = None
        self._requirements = None

    def is_pipeline(self):
        return self._mf.is_pipeline(self._name)

    def _get_file(self, name):
        filename = self._path + os.sep + name
        if os.path.exists(filename):
            return filename

    def __str__(self):
        txt = "Rule **" + self.name + "**:\n" + self.description
        return txt

    def _get_path(self):
        return self._path
    path = property(_get_path, doc="full path to the module directory")

    def _get_config(self):
        # The default config file for that module
        filename = self._get_file("config.yaml")
        if filename is None:
            # or the sequana default config file
            filename = self._get_file("../config.yaml")
        return filename
    config = property(_get_config,
                      doc="full path to the config file of the module")

    def _get_readme(self):
        return self._get_file("README.rst")
    readme = property(_get_readme,
                      doc="full path to the README file of the module")

    def _get_overview(self):
        result = "no information. For developers: please fix the pipeline "
        result += "README.rst file by adding an :Overview: field"
        for this in self.description.split("\n"):
            if this.startswith(':Overview:'):
                try:
                    result = this.split(":Overview:")[1].strip()
                except:
                    result += "Bad format in :Overview: field"
        return result
    overview = property(_get_overview)

    def _get_snakefile(self):
        if self._snakefile is not None:
            return self._snakefile

        if self._get_file("Snakefile"):
            self._snakefile = self._get_file("Snakefile")
        elif self._get_file("Snakefile." + self.name):
            self._snakefile = self._get_file("Snakefile." + self.name)
        elif self._get_file(self.name + '.rules'):
            self._snakefile = self._get_file(self.name + ".rules")
        else:
            print("//Snakefile for %s not found" % self.name)
        return self._snakefile
    snakefile = property(_get_snakefile,
                         doc="full path to the Snakefile file of the module")

    def _get_name(self):
        return self._name
    name = property(_get_name, doc="name of the module")

    def _get_requirements(self):
        if self._requirements is not None:
            return self._requirements
        if self._get_file("requirements.txt"):
            self._requirements = self._get_file("requirements.txt")
            return self._requirements
    requirements = property(_get_requirements, doc="list of requirements")

    def is_executable(self, verbose=False):
        if self.requirements is None:
            return True

        executable = True
        import easydev

        # reads the file and interpret it to figure out the
        # executables/packages and pipelines required
        pipelines = []
        with open(self.requirements, "r") as fh:
            data = fh.read()
            datalist = [this.strip() for this in data.split("\n")
                        if this.strip() not in [""]]
            reqlist = []
            for this in datalist:
                if this.startswith('-'):
                    req = this.split("-")[1].split()[0].strip()
                    if req.startswith("["):
                        req = req.replace("[", "")
                        req = req.replace("]", "")
                        pipelines.append(req)
                    else:
                        reqlist.append(req)

        # Check the pipelines independently
        for pipeline in pipelines:
            Module(pipeline).check()

        for req in reqlist:
            # It is either a Python package or an executable
            try:
                easydev.shellcmd("which %s" % req)
                if verbose:
                    print("Found %s executable" % req)
            except:
                # is this a Python code ?
                if len(easydev.get_dependencies(req)) == 0:
                    if verbose:
                        print("%s not found !!" % req)
                    executable = False
                else:
                    if verbose:
                        print("%s python package" % req)
        return executable

    def check(self, mode="warning"):
        if self.is_executable(verbose=False):
            #print("All requirements fulfilled for %s pipeline/rule." % self._name)
            return
        else:
            self.is_executable(verbose=True)
            txt = "Some executable or Python packages are not available:\n"
            txt += "Some functionalities may not work "
            if mode == "warning":
                print(txt)
            elif mode == "error":
                txt += "Use \n conda install missing_package_name"
                raise ValueError(txt)

    def _get_description(self):
        try:
            with open(self.readme) as fh:
                self._description = fh.read()
        except:
                self._description = "no description"
        return self._description
    description = property(_get_description,
                           doc=("Content of the README file associated with "
                                "the module.\n\n::\n"
                                "   from sequana import Module\n"
                                "   m = Module('dag')\n"
                                "   print(m.description)\n"))

    def onweb(self):
        # TOD: automatic switch
        from easydev import onweb
        if "rules" in self._path:
            suffix = self.snakefile.split("rules/")[1]
            suffix = suffix.rsplit("/", 1)[0]
            onweb("http://github.com/sequana/sequana/tree/"
                  "master/sequana/rules/%s" % suffix)
        else:
            suffix = self.snakefile.split("pipelines/")[1]
            suffix = suffix.rsplit("/", 1)[0]
            onweb("http://github.com/sequana/sequana/tree/"
                  "master/sequana/pipelines/%s" % self.name)


def _get_modules_snakefiles():
    modules = {}
    for name in ModuleFinder().names:
        if Module(name).snakefile:
            modules[name] = Module(name).snakefile
    return modules

#: dictionary with module names as keys and fullpath to the Snakefile as values
modules = _get_modules_snakefiles()

#: list of pipeline names found in the list of modules
pipeline_names = [m for m in modules if Module(m).is_pipeline()]


class SequanaConfig(object):
    """Reads YAML (or json) config file and ease access to its contents

    This can also be used to check the validity of the config file

    ::

        >>> sc = SequanaConfig(config)
        >>> sc.config.pattern == "*.fastq.gz"
        True

    Input files should be stored into::

        samples:
            - file1: FILE1
            - file2: FILE2

    The second file may be optional.


    Empty strings in a config are interpreted as None but SequanaConfig will
    replace  None with empty strings, which is probably what was expected from
    the user. Similarly, in snakemake when settings the config file, one
    can override a value with a False but this is interepted as "False"
    This will transform back the "True" into True

    """
    def __init__(self, data=None, test_requirements=True,
                 converts_none_to_str=True, mode="NGS"):
        """Could be a json or a yaml

        :param str filename: filename to a config file in json or yaml format.

        mode NGS means the following fields are expected: project, sample,
        file1, file2
        """
        if data is None:
            config = None
            self.config = AttrDict()
            mode = "others"
        elif isinstance(data, str):
            config = load_configfile(data)
            self.config = AttrDict(**config)
        else:
            self.config = AttrDict(**data)

        if converts_none_to_str and mode == "NGS":
            self._set_none_to_empty_string(self.config)

        self._converts_boolean(self.config)

        if test_requirements and mode == "NGS":
            requirements = ["input_directory", "samples", "samples:file1", "samples:file2"]
            # converts to dictionary ?
            for this in requirements:
                this = this.split(":")[0]
                assert this in self.config.keys(),\
                    "Your config must contain %s" % this



    def save(self, filename="config.yaml"):
        """Export config into YAML file

        Save the config file into a config file in YAML format.
        The initial input file is in JSON format
        """
        # Using yaml.dump, the cases with arguments starting with - prevents
        # yaml.dump to understand that it is a field and so quotes are missing.
        # So we do it by hand assuming two levels at most

        # We do not put quotes except when we find a space, a % (e.g for
        # templating with "%(param)s" case, a \\ (e.g. bwa_ref case)
        txt = ""

        def special(string):
            if isinstance(string, str) is False:
                return False
            if " " in string:
                return True
            if "%" in string or "\\" in string:
                return True
            return False

        for k1 in sorted(self.config.keys()):
            v1 = self.config[k1]
            if isinstance(v1, dict):
                txt += "%s:\n" % k1
                for k2, v2 in v1.items():
                    if isinstance(v2, dict):
                        txt += "    %s:\n" % k2
                        for k3, v3 in v2.items():
                            if special(v3):
                                txt += "        %s: '%s'\n" % (k3, v3)
                            else:
                                txt += "        %s: %s\n" % (k3, v3)
                    elif isinstance(v2, list):
                        txt += "    %s: '%s'\n" % (k2, self.config[k1][k2])
                    else:
                        if special(v2):
                            txt += '    %s: "%s"\n' % (k2, v2)
                        else:
                            txt += '    %s: %s\n' % (k2, v2)
            elif isinstance(v1, list):
                txt += "%s:\n" % k1
                for item in self.config[k1]:
                    if special(item):
                        txt += "    - '%s'\n" % item
                    else:
                        txt += "    - %s\n" % item
            else:
                if special(v1):
                    txt += '%s: %s\n' % (k1, v1)
                else:
                    txt += '%s: "%s"\n' % (k1, v1)
            txt += "\n"
            with open(filename, "w") as fout:
                fout.write(txt)

    def _converts_boolean(self, subdic):
        for key, value in subdic.items():
            if isinstance(value, dict):
                subdic[key] = self._converts_boolean(value)
            else:
                if value in ['True', 'true', 'TRUE']:
                    subdic[key] = True
                if value in ['False', 'false', 'FALSE']:
                    subdic[key] = False
        return subdic

    def _set_none_to_empty_string(self, subdic):
        # recursively set parameter (None) to ""
        for key, value in subdic.items():
            if isinstance(value, dict):
                subdic[key] = self._set_none_to_empty_string(value)
            else:
                if value is None:
                    subdic[key] = ""
        return subdic

    @staticmethod
    def from_dict(dic):
        return SequanaConfig(dic)

    def check(self, requirements_dict):
        """a dcitionary in the form

        ::

            rule_name:["param1", ":param2", ":kwargs:N"]

        in a config::

            param1: 1
            param2:
                - subparam1
                - subparam2:
                    - subparam3: 10


        """
        dd = requirements_dict

        # check that each field request is present
        for k, vlist in dd.items():
            k = k.replace('__sequana__', '')
            if k not in self.config.keys():
                raise KeyError("Expected section %s not found" % k)
            else:
                for item in vlist:
                    # with : , this means a sub field field
                    if item.startswith(":") and item.count(":") == 1:
                        assert item[1:] in self.config[k].keys()
                    elif item.count(":") > 1:
                        raise NotImplementedError(
                            "2 hierarchy checks in config not implemented ")
                    else:
                        # without : , this means a normal field so item is
                        # actually a key here
                        assert item in self.config.keys()

    def get(self, field, default=None, cfg=None):
        if cfg is None:
            cfg = self.config
        if ":" in field:
            level1, level2 = field.split(":", 1)
            if level1 not in cfg.keys():
                raise ValueError("first level key (%s) not found " % level1)
            print("Found :")
            return self.get(level2, default=default, cfg=cfg[level1])
        else:
            if isinstance(cfg, dict) and field in cfg.keys():
                return cfg[field]
            else:
                return default


class PipelineManager(object):
    def __init__(self, name, config, pattern="*.fastq.gz"):
        cfg = SequanaConfig(config)
        cfg.config.pipeline_name = name

        # Default mode is the glob/pattern. 
        glob_dir = cfg.config.input_directory + "/" + \
                   cfg.get('pattern', pattern)
        self.ff = FastQFactory(glob_dir)

        R1 = [1 for this in self.ff.filenames if "_R1_" in this]
        R2 = [1 for this in self.ff.filenames if "_R2_" in this]
        if len(R2) == 0:
            self.paired = False
        else:
            if R1 == R2:
                self.paired = True
            else:
                raise ValueError("Mix of paired and single-end data sets not implemented yet")

        # Note, however, that another mode is the samples.file1/file2 . If
        # provided, filenames is not empty and superseeds
        # the previous results (with the glob). Here only 2 files are provided
        # at most
        if not self.ff:
            filenames = self._get_filenames(cfg.config)
            if len(filenames):
                self.ff = FastQFactory(filenames)
                if len(filenames) == 2:
                    self.paired = True
                else:
                    self.paired = False

        ff = self.ff  # an alias
        self.samples = {tag: [ff.get_file1(tag), ff.get_file2(tag)]
                        if ff.get_file2(tag) else [ff.get_file1(tag)]
                        for tag in ff.tags}

        if len(ff.tags) == 0:
            raise ValueError("Could not find fastq.gz files with valid format "
                             "(NAME_R1_<SUFFIX>.fastq.gz where <SUFFIX> is "
                             "optional")
        elif len(ff.tags) == 1:
            self.mode = "nowc"  # no wildcard
            self.sample = ff.tags[0]
            self.basename = self.sample + "/%s/" + self.sample
        else:
            self.mode = "wc"
            self.sample = "{sample}"
            self.basename = "{sample}/%s/{sample}"

        # finally, keep track of the config file
        self.config = cfg.config

    def getname(self, rulename, suffix=None):
        if suffix is None:
            suffix = ""
        return self.basename % rulename + suffix

    def getreportdir(self, acronym):
        return self.sample + "/report_" + acronym + "_" + self.sample + "/"

    def getwkdir(self, rulename):
        return self.sample + "/" + rulename

    def getrawdata(self):
        if self.mode == "nowc":
            return self.ff.realpaths
        else:
            return lambda wildcards: self.samples[wildcards.sample]

    def _get_filenames(self, cfg):
        filenames = []
        file1 = cfg.samples.file1
        file2 = cfg.samples.file2
        if file1:
            if os.path.exists(file1):
                filenames.append(file1)
            else:
                raise FileNotFoundError("%s not found" % file1)
        if file2:
            if os.path.exists(file2):
                filenames.append(file2)
            else:
                raise FileNotFoundError("%s not found" % file2)
        return filenames






def sequana_check_config(config, globs):
    s = SequanaConfig.from_dict(config)
    dic = dict([(k, v) for k, v in globs.items()
                if k.startswith("__sequana__")])
    s.check(dic)


def message(mes):
    """Dedicated print function to include in Snakefiles

    In a Snakefile, the stand print function may interfer with other process
    An example is the creation of the dag file. Not sure this is a bug but
    meanwhile, one must use this function to print information.

    This adds the // -- characters in front of the prin statements."""
    from easydev.console import purple
    print("// -- " + purple(mes))


class DOTParser(object):
    """Utility to parse the dot returned by Snakemake and add URLs automatically

    ::

        from sequana import sequana_data
        from sequana.snaketools import DOTParser

        filename = sequana_data("test_dag.dot", "testing")
        dot = DOTParser(filename, {"fastqc": "fastqc.html"})

        # creates test_dag.ann.dot locally
        dot.add_urls()

    """
    def __init__(self, filename):
        self.filename = filename

    def add_urls(self, output_filename=None, mapper={}):
        with open(self.filename, "r") as fh:
            data = fh.read()

        if output_filename is None:
            import os
            output_filename = os.path.basename(self.filename)

        with open(output_filename.replace(".dot", ".ann.dot"), "w") as fout:
            indices_to_drop = []
            for line in data.split("\n"):
                if "[label =" not in line:
                    if " -> " in line:
                        label = line.split(" -> ")[0].strip()
                        if label not in indices_to_drop:
                            fout.write(line + "\n")
                    else:
                        line = line.replace("dashed", "")
                        fout.write(line + "\n")
                else:
                    separator = "color ="
                    lhs, rhs = line.split(separator)
                    name = lhs.split("label =")[1]
                    name = name.replace(",", "")
                    name = name.replace('"', "")
                    name = name.strip()
                    if name in ['dag', 'conda', "rulegraph"]:
                        index = lhs.split("[")[0]
                        indices_to_drop.append(index.strip())
                    elif name in mapper.keys():
                        url = mapper[name]
                        newline = lhs + (' URL="%s"'
                                         ' target="_parent", ') % url
                        newline += separator + rhs
                        newline = newline.replace("dashed", "")
                        newline = newline.replace('];', 'color="blue"];')
                        fout.write(newline + "\n")
                    else:
                        newline = line.replace('];', 'color="orange"];')
                        fout.write(newline + "\n")


def __get_tagname(filename):
    """Given a fullpath name, remove extension and prefix and return the name

    .. deprecated::

    """
    raise ValueError("deprecated")
    import os
    # This should always work
    name = os.path.split(filename)[1].split('.', 1)[0]
    return name


class FileFactory(object):
    """

    ::

        ff = FileNameFactory("../*")
        ff.dataset

    special case if ends in .fastq.gz


    definition to use::

        >>> fullpath = /home/user/test/A.fastq.gz
        >>> dirname(fullpath)
        '/home/user/test'
        >>> basename(fullpath)
        'A.fastq.gz'
        >>> realpath(fullpath) # is .., expanded to /home/user/test
        >>> all_extensions
        "fastq.gz"
        >>> extensions
        ".gz"

    """
    def __init__(self, pattern):
        self.pattern = pattern
        if isinstance(pattern, str):
            self._glob = glob.glob(pattern)
        elif isinstance(pattern, list):
            self._glob = pattern[:]

    def _get_realpath(self):
        return [os.path.realpath(filename) for filename in self._glob]
    realpaths = property(_get_realpath)

    def _get_basenames(self):
        return [os.path.split(filename)[1] for filename in self._glob]
    basenames = property(_get_basenames)

    def rstrip(self, ext):
        return [filename.rstrip(ext) for filename in self.basenames]

    def _get_filenames(self):
        return [this.split(".")[0] for this in self.basenames]
    filenames = property(_get_filenames)

    def _pathnames(self):
        pathnames = [os.path.split(filename)[0] for filename in self._glob]
        return pathnames
    pathnames = property(_pathnames)

    def _pathname(self):
        pathname = set(self.pathnames)
        if len(pathname) == 1:
            return list(pathname)[0] + os.sep
        else:
            raise ValueError("found more than one pathname")
    pathname = property(_pathname)

    def _get_extensions(self):
        filenames = [os.path.splitext(filename)[1] for filename in self._glob]
        return filenames
    extensions = property(_get_extensions)

    def _get_all_extensions(self):
        filenames = [this.split('.', 1)[1] if "." in this else ""
                     for this in self.basenames]
        return filenames
    all_extensions = property(_get_all_extensions)


class FastQFactory(FileFactory):
    """


    """
    def __init__(self, pattern, extension="fastq.gz", strict=True,
                 rstrip_underscore=True):
        """

        :param strict: if true, the pattern _R1_ or _R2_ must be found
        """
        super(FastQFactory, self).__init__(pattern)

        if len(self.filenames) == 0:
            raise ValueError("No files found with the requested pattern (%s)"
                             % pattern)

        # Check the extension of each file (fastq.gz by default)
        for this in self.all_extensions:
            assert this.endswith(extension), \
                "Expecting file with %s extension. Found %s" % (extension,
                                                                this)

        # identify a possible tag
        self.tags = []
        for this in self.filenames:
            if '_R1_' in this:
                self.tags.append(this.split('_R1_', 1)[0])
            elif '_R2_' in this:
                self.tags.append(this.split('_R2_', 1)[0])
            elif strict is True:
                # Files must have _R1_ and _R2_
                raise ValueError('FastQ filenames must contain _R1_ or _R2_')
        self.tags = list(set(self.tags))

        self.short_tags = [x.split("_")[0] for x in self.tags]


    def _get_file(self, tag, rtag):
        assert rtag in ["_R1_", "_R2_"]

        if tag is None:
            if len(self.tags) == 1:
                tag = self.tags[0]
            elif len(self.tags) > 1:
                raise ValueError("Ambiguous tag. You must provide one "
                                 "(sequana.FastQFactory)")
        else:
            assert tag in self.tags, 'invalid tag'

        # Note that the rtag + "_" is to prevent the two following pairs to
        # be counted as one sample:
        # project-name_1
        # project-name-bis_1
        candidates = [realpath for filename, realpath in
                      zip(self.filenames, self.realpaths)
                      if rtag in filename and filename.startswith(tag+"_")]

        if len(candidates) == 0 and rtag == "_R2_":
            # assuming there is no R2
            return None
        elif len(candidates) == 1:
            return candidates[0]
        else:
            print('Found too many candidates: %s ' % candidates)
            raise ValueError("Files must have the tag _R1_ or _R2_ (%s)" % tag)

    def get_file1(self, tag=None):
        return self._get_file(tag, "_R1_")

    def get_file2(self, tag=None):
        return self._get_file(tag, "_R2_")



def init(filename, namespace):
    """Defines the global variable __snakefile__ inside snakefiles

    If not already defined, __snakefile__ is created to hold the name of the
    pipeline. We also define initialise these variables :

        * expected_output as an empty list
        * toclean as an empty  list

    """
    # Create global name for later
    if "__snakefile__" in namespace.keys():
        pass
    else:
        namespace['__snakefile__'] = filename
        namespace['expected_output'] = []
        namespace['toclean'] = []

    # check requirements
    Module(filename.replace(".rules", "")).check('error')


def create_recursive_cleanup(filename=".sequana_cleanup.py"):
    with open(filename, "w") as fh:
        fh.write("""
import subprocess
import glob
import os
for this in glob.glob("*"):
    if os.path.isdir(this) and this not in ["fastq_sampling", "report"]:
        print(" --- Cleaning up %s directory" % this)
        subprocess.Popen(["python", ".sequana_cleanup.py"], cwd=this)

from easydev import shellcmd
shellcmd("rm  README runme.sh config.yaml" )
shellcmd("rm  stats.txt *fa" )
shellcmd("rm  dag.svg" )
shellcmd("rm  *.rules" )

# We can further clean additional files

""")



def create_cleanup(targetdir):
    """A script to include in directory created by the different pipelines"""
    with open(targetdir + os.sep + ".sequana_cleanup.py", "w") as fout:
        fout.write("""
import glob
import os
import shutil
from easydev import shellcmd
import time

directories = glob.glob("*")

for this in directories:
    if os.path.isdir(this) and this not in ['logs'] and 'report' not in this:
        print('Deleting %s' % this)
        time.sleep(0.1)
        shellcmd("rm -rf %s" % this)
shellcmd("rm -f  README config.yaml snakejob.*")
shellcmd("rm -f .sequana_cleanup.py")
shellcmd("rm -rf .snakemake")
""")
