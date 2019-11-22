"""Configuration for SPHINX to generate documentation."""

import sys
import inspect

# -*- coding: utf-8 -*-
#
# Configuration file for the Sphinx documentation builder.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os

sys.path.insert(
    0,
    os.path.abspath(
        os.path.join(
            os.path.dirname(__file__),
            '../../../PythonPackage'
        )
    )
)


# -- Project information -----------------------------------------------------

project = '(P)lanetary (O)rbital (E)volution due to (T)ides'
#pylint: disable=redefined-builtin
copyright = '2019, Kaloyan Penev'
#pylint: enable=redefined-builtin
author = 'Kaloyan Penev'

# The short X.Y version
version = ''
# The full version, including alpha/beta/rc tags
release = ''


# -- General configuration ---------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.ifconfig',
    'sphinx.ext.viewcode',
    'sphinx.ext.githubpages',
    'sphinx.ext.napoleon',
    'sphinx.ext.inheritance_diagram',
    'nbsphinx',
    'breathe',
    'exhale'
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = None

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
# html_theme_options = {}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Custom sidebar templates, must be a dictionary that maps document names
# to template names.
#
# The default sidebars (for documents that don't match any pattern) are
# defined by theme itself.  Builtin themes are using these templates by
# default: ``['localtoc.html', 'relations.html', 'sourcelink.html',
# 'searchbox.html']``.
#
# html_sidebars = {}


# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'PlanetaryOrbitalEvolutionduetoTidesdoc'


# -- Options for LaTeX output ------------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',

    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (
        master_doc,
        'PlanetaryOrbitalEvolutionduetoTides.tex',
        '(P)lanetary (O)rbital (E)volution due to (T)ides Documentation',
        'Kaloyan Penev',
        'manual'
    ),
]


# -- Options for manual page output ------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (
        master_doc,
        'planetaryorbitalevolutionduetotides',
        '(P)lanetary (O)rbital (E)volution due to (T)ides Documentation',
        [author],
        1
    )
]


# -- Options for Texinfo output ----------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (
        master_doc,
        'PlanetaryOrbitalEvolutionduetoTides',
        '(P)lanetary (O)rbital (E)volution due to (T)ides Documentation',
        author,
        'PlanetaryOrbitalEvolutionduetoTides',
        'One line description of project.',
        'Miscellaneous'
    ),
]


# -- Options for Epub output -------------------------------------------------

# Bibliographic Dublin Core info.
epub_title = project

# The unique identifier of the text. This can be a ISBN number
# or the project homepage.
#
# epub_identifier = ''

# A unique identification for the text.
#
# epub_uid = ''

# A list of files that should not be packed into the epub file.
epub_exclude_files = ['search.html']


# -- Extension configuration -------------------------------------------------

# -- Options for intersphinx extension ---------------------------------------

# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {'https://docs.python.org/': None}

# -- Options for todo extension ----------------------------------------------

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = True

# -- Options for autodoc extension -------------------------------------------

autodoc_default_flags = ['members',
                         'undoc-members',
                         'show-inheritance']

#Napolean extension defined names.
#pylint: disable=invalid-name
napoleon_include_private_with_doc = True
napoleon_include_special_with_doc = True
napoleon_include_init_with_doc = True
#pylint: enable=invalid-name

inheritance_graph_attrs = dict(rankdir="TB",
                               fontsize="24",
                               ratio='auto',
                               size='120')

#Call signature defined by SPHINX autodoc plugin.
#pylint: disable=too-many-arguments
#pylint: disable=unused-argument
def add_inheritance_diagram(app, what, name, obj, options, lines):
    """Add an inheritance diagram for all classes."""

    if what == 'module':
        class_list = [member[0]
                      for member in inspect.getmembers(sys.modules[name],
                                                       inspect.isclass)]
        if class_list:
            lines.insert(0, '')
            lines.insert(
                0, '.. inheritance-diagram:: '
                +
                ' '.join(class_list)
            )
            lines.insert(0, '=========================')
            lines.insert(0, 'Class Inheritance Diagram')
    elif what == 'class':
        lines.insert(0, '')
        lines.insert(0,
                     '.. inheritance-diagram:: ' + name)
#pylint: enable=too-many-arguments

def setup(app):
    """Connect handler for adding inheritance diagrams."""

    app.add_stylesheet('unlimited_width.css')
    app.connect('autodoc-process-docstring', add_inheritance_diagram)


doxygen_xml = os.path.abspath(
    os.path.join(
        os.path.dirname(__file__),
        '../../doxygen/build/xml'
    )
)

breathe_projects = {
    'C++ library': doxygen_xml
}


breathe_default_project = "C++ library"

breathe_default_members = ('members',
                           'protected-members',
                           'private-members',
                           'undoc-members')

# Setup the exhale extension
exhale_args = {
    # These arguments are required
    "containmentFolder":     "./api",
    "rootFileName":          "library_root.rst",
    "rootFileTitle":         "Library API",
    "doxygenStripFromPath":  "..",
    # Suggested optional arguments
    "createTreeView":        True
}
