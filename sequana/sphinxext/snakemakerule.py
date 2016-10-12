# -*- coding: utf-8 -*-
r"""
    sphinx.ext.inheritance_diagram
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Defines a docutils directive for inserting simple docstring extracting from 
    a snakemake rule (from sequana project).

    .. snakemakerule: test

"""

from __future__ import absolute_import, print_function
from six import class_types
import re
from docutils.parsers.rst import  directives
from docutils import nodes


def get_rule_doc(name):
    """Decode and return lines of the docstring(s) for the object."""
    from sequana import Module
    import snakemake
    rule = Module(name)
    wf = snakemake.Workflow(rule.path + "/%s.rules" %  name)
    wf.include(rule.path + "/%s.rules" % name)
    docstring = list(wf.rules)[0].docstring
    return docstring


from docutils.nodes import Body, Element


class snakemake_base(Body, Element):
    def dont_traverse(self, *args, **kwargs):
        return []


class snakemake_rule(snakemake_base):
    pass


def run(content, node_class, state, content_offset):
    text = get_rule_doc(content[0])
    node = node_class("")  # shall we add something here ?
    node.rule_docstring = get_rule_doc(content[0])
    state.nested_parse(content, content_offset, node)
    return [node]


def snakemake_rule_directive(name, arguments, options, content, lineno,
                         content_offset, block_text, state, state_machine):
    return run(content, snakemake_rule, state, content_offset)


def setup(app):
    #app.add_autodocumenter(RuleDocumenter)
    app.add_directive('snakemakerule', snakemake_rule_directive, True, (0,0,0))

    # Add visit/depart methods to HTML-Translator:
    def visit_perform(self, node):
        # Ideally, we should use sphinx but this is a simple temporary solution
        from docutils import core
        from docutils.writers.html4css1 import Writer
        w = Writer()
        res = core.publish_parts(node.rule_docstring, writer=w)['html_body']
        self.body.append('<div class="snakemake">' + res + '</div>' )
        node.children = []

    def depart_perform(self, node):
        node.children = []

    def visit_ignore(self, node):
        node.children = []

    def depart_ignore(self, node):
        node.children = []

    app.add_node(snakemake_rule,
                 html=(visit_perform, depart_perform),
                 latex=(visit_ignore, depart_ignore))







