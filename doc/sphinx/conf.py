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
copyright = '2023'
author = 'the authors of the ASPECT manual'

# The full version, including alpha/beta/rc tags
with open('../../VERSION', 'r') as file:
    release = file.read().replace('\n','')

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "myst_parser",
    "sphinxcontrib.bibtex",
    "sphinxcontrib.tikz",
    "sphinxcontrib.cairosvgconverter"
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
html_logo = "_static/images/aspect_logo.png"
html_title = "ASPECT " + release
html_theme = 'sphinx_book_theme'
html_theme_options = {
    "collapse_navigation": True,
    "navigation_depth": 3,
    "show_toc_level": 3,
    "repository_url": "https://github.com/geodynamics/aspect/",
    "repository_branch": "main",
    "path_to_docs":"doc/sphinx/",
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/geodynamics/aspect",
            "icon": "fab fa-github-square",
        },
    ],
    "show_navbar_depth": 1,
    "use_repository_button": True,
    "use_edit_page_button": True,
    "use_issues_button": True,
    "home_page_in_toc": True,
    "logo": {
        "text": "ASPECT " + release,
    },
    "primary_sidebar_end": "navbar_end.html"
}

latex_engine = 'xelatex'

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

linkcheck_allowed_redirects = {r'https://doi.org/*': r'https://*',
                               r'https://youtu.be/*': r'https://youtube.com/*',
                               r'https://www.github.com/*': r'https://github.com/*'}

latex_elements = {
# The paper size ('letterpaper' or 'a4paper').
'papersize': 'a4paper',
'babel': '\\usepackage[english]{babel}',

# The font size ('10pt', '11pt' or '12pt').
'pointsize': '11pt',

# Additional stuff for the LaTeX preamble.
'preamble': r'''
% to allow for the degree symbol
\usepackage{gensymb}
\usepackage{siunitx}
\usepackage{etoolbox}
\pretocmd{\hyperlink}{\protect}{}{}

\usepackage{graphicx}
\setlength{\parindent}{0pt}
\setlength{\parskip}{5pt}
\usepackage{textpos}

% use the listings package for code snippets
\usepackage{listings}

\definecolor{dkgreen}{rgb}{0,0.6,0}
\definecolor{gray}{rgb}{0.5,0.5,0.5}
\definecolor{mauve}{rgb}{0.58,0,0.82}
  ''',

# Additional stuff for the LaTeX preamble.
'maketitle': r'''
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% START OF CIG MANUAL COVER TEMPLATE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This should be pasted at the start of manuals and appropriate strings entered at locations indicated with FILL.
% Be sure the TeX file includes the following packages.
% \usepackage{graphicx}
% \usepackage{times}
% \usepackage{textpos}

\definecolor{dark_grey}{gray}{0.3}
\definecolor{aspect_blue}{rgb}{0.3125,0.6875,0.9375}
\definecolor{aspect_red}{rgb}{0.522,0.0,0.0}

%LINE 1%
{
\renewcommand{\familydefault}{\sfdefault}
\newcommand{\aspect}{\textsc{ASPECT}}

\pagenumbering{gobble}
\begin{center}
\resizebox{\textwidth}{!}{\textcolor{dark_grey}{\fontfamily{\sfdefault}\selectfont
COMPUTATIONAL INFRASTRUCTURE FOR GEODYNAMICS (CIG)
}}

\hrule

%LINE 2%
\color{dark_grey}
\rule{\textwidth}{2pt}

%LINE 3%
\color{dark_grey}
% FILL: additional organizations
% e.g.: {\Large Organization 1\\Organization 2}
{\Large }
\end{center}

%COLOR AND CODENAME BLOCK%
\begin{center}
\resizebox{\textwidth}{!}{\colorbox
% FILL: color of code name text box
% e.g. blue
{aspect_red}{\fontfamily{\rmdefault}\selectfont \textcolor{white} {
% FILL: name of the code
% You may want to add \hspace to both sides of the codename to better center it, such as:
% \newcommand{\codename}{\hspace{0.1in}CodeName\hspace{0.1in}}
\hspace{0.1in}\aspect{}\hspace{0.1in}
}}}
\\[12pt]
{\Large Advanced Solver for Planetary Evolution, Convection, and Tectonics}
\end{center}

{
  \parindent 0pt
  % Place the Logo into a box of width 0 and position it at the given
  % (x,y) within this box.
  \begin{textblock*}{0in}(0.0in,0.8in)
    \begin{center}
      \vspace{1em}
        \includegraphics[width=4.5in]{{aspect_logo}.png}
      \hspace{5em}
    \end{center}
  \end{textblock*}

  % Then overlay some text into another textblock* environment that is
  % as wide as the page. We use \raggedright, so everything is shown
  % on the right side of the page. The y-coordinate of the anchor
  % (0.3in) is the same as for the picture, aligning the text nicely
  % to the image.
  \begin{textblock*}{\textwidth}(0in,0.3in)
    \vspace{1em}
    \color{dark_grey}
    \hfill{\Huge \fontfamily{\sfdefault}\selectfont User Manual \\
      \raggedleft \huge \fontfamily{\sfdefault}\selectfont Version
      % keep the following line as is so that we can replace this using a script:
      %VERSION-INFO%
      \\
      \large(generated \today)
      \\[16pt]
        {\Large
          Wolfgang Bangerth \\
          Juliane Dannberg \\
          Menno Fraters \\
          Rene Gassm{\"o}ller \\
          Anne Glerum \\
          Timo Heister \\
          Bob Myhill \\
          John Naliboff\\}
    }
  \end{textblock*}
}


%AUTHOR(S) & WEBSITE%
\null
\vfill
\color{dark_grey}
{\fontfamily{\sfdefault}\selectfont
% FILL: author list
% e.g. Author One\\Author Two\\Author Three\\
% be sure to have a newline (\\) after the final author
\large
\noindent with contributions by: \\
    Jacqueline Austermann,
    Magali Billen,
    Markus B{\"u}rg,
    Thomas Clevenger,
    Samuel Cox,
    William Durkin,
    Grant Euen,
    Thomas Geenen,
    Ryan Grove,
    Eric Heien,
    Ludovic Jeanniot,
    Louise Kellogg,
    Scott King,
    Martin Kronbichler,
    Marine Lasbleis,
    Haoyuan Li,
    Shangxin Liu,
    Hannah Mark,
    Elvira Mulyukova,
    Bart Niday,
    Jonathan Perry-Houts,
    Elbridge Gerry Puckett,
    Tahiry Rajaonarison,
    Fred Richards,
    Jonathan Robey,
    Ian Rose,
    Max Rudolph,
    Stephanie Sparks,
    D.~Sarah Stamps,
    Cedric Thieulot,
    Wanying Wang,
    Iris van Zelst,
    Siqi Zhang\\
\vspace{0.5em}
}

{\noindent
{\fontfamily{\sfdefault}\selectfont \href{https://geodynamics.org}{geodynamics.org}}
}

%LINE%
{\noindent
\color{dark_grey}
\rule{\textwidth}{2pt}
}

}

\pagebreak
\pagenumbering{arabic}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   END OF CIG MANUAL COVER TEMPLATE    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'''
}

latex_elements['maketitle'] = latex_elements['maketitle'].replace("%VERSION-INFO%", release)
