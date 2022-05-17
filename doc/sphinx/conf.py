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

project = 'ASPECT'
copyright = '2022'
author = 'Wolfgang Bangerth, Juliane Dannberg, Menno Fraters, Rene Gassmöller, Anne Glerum, Timo Heister, Bob Myhill, John Naliboff; with contributions by: Jacqueline Austermann, Magali Billen, Markus B{\"u}rg, Thomas Clevenger, Samuel Cox, William Durkin, Grant Euen, Thomas Geenen, Ryan Grove, Eric Heien, Ludovic Jeanniot, Louise Kellogg, Scott King, Martin Kronbichler, Marine Lasbleis, Haoyuan Li, Shangxin Liu, Hannah Mark, Elvira Mulyukova, Bart Niday, Jonathan Perry-Houts, Elbridge Gerry Puckett, Tahiry Rajaonarison, Fred Richards, Jonathan Robey, Ian Rose, Max Rudolph, Stephanie Sparks, D.~Sarah Stamps, Cedric Thieulot, Wanying Wang, Iris van Zelst, Siqi Zhang'

# The full version, including alpha/beta/rc tags
release = '2.4.0-pre'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "myst_parser",
    "sphinxcontrib.bibtex",
    "sphinxcontrib.tikz"
]
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "fieldlist",
    "substitution",
    "dollarmath",
    "amsmath",
]

tikz_proc_suite = "GhostScript"
# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = [
    "_build",
]

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_logo = "../logo/unlabeled_logo_small.png"
html_theme = 'pydata_sphinx_theme'
html_theme_options = {
    "collapse_navigation": True,
    "navigation_depth": 3,
    "show_toc_level": 3,
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/geodynamics/aspect",
            "icon": "fab fa-github-square",
        },
    ],
    "navbar_start": ["navbar-logo"],
    "footer_items": ["last-updated"],
    "use_edit_page_button": True,
}
html_sidebars = {
    "**": ["search-field", "sidebar-nav-bs", "sidebar-ethical-ads", "sidebar-cig"]
}
bibtex_bibfiles = ["references.bib"]
bibtex_default_style = "alpha"
bibtex_reference_style = "author_year"
numfig = True
html_context = {
    # "github_url": "https://github.com", # or your GitHub Enterprise interprise
    "github_user": "geodynamics",
    "github_repo": "aspect",
    "github_version": "main",
    "doc_path": "doc/sphinx",
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
html_last_updated_fmt = ""
