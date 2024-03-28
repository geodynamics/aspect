# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'Geodynamic World Builder'
copyright = '2024, The authors of the Geodynamic World Builder'
# The full version, including alpha/beta/rc tags
release = '@WORLD_BUILDER_VERSION@'
html_title = "Manual GWB @WORLD_BUILDER_VERSION@"
version: "2"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
'myst_parser',
'sphinx.ext.todo',
'sphinx.ext.imgconverter',
'sphinx_design',
'sphinx_copybutton',
'sphinxcontrib.bibtex',
]

bibtex_default_style = 'plain'

bibtex_bibfiles = ['bibliography.bib']

myst_enable_extensions = ["colon_fence","dollarmath"]

# Breathe Configuration
breathe_default_project = "GWB"

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_book_theme'

html_logo = "@CMAKE_CURRENT_BINARY_DIR@/sphinx/_static/images/world_builder_logo_v4_b_text.png"
latex_logo = "@CMAKE_CURRENT_BINARY_DIR@/sphinx/_static/images/world_builder_logo_v4_b_text.png"

html_theme_options = {
    "navigation_with_keys":"true",
    "repository_url": "https://github.com/GeodynamicWorldBuilder/WorldBuilder/",
    "repository_branch": "main",
    "path_to_docs":"doc/sphinx/",
    "use_repository_button": True,
    "use_edit_page_button": True,
    "use_issues_button": True,
    "home_page_in_toc":False,
    "extra_footer":"Built with <a href='https://www.sphinx-doc.org/'>Sphinx</a> using a theme provided by <a href='https://ebp.jupyterbook.org/'>Executable Book Project</a>",
}


todo_include_todos = True

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['@CMAKE_CURRENT_SOURCE_DIR@/sphinx/_static/']