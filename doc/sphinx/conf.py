# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import re

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
    "sphinx_design",
    "sphinx_tags",
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


# Map units used by the LaTeX 'siunitx' package to text that
# can directly be used in LaTeX formulas:
_siunitx_symbols = {
    "degree": r"^\circ",
    "J": "J",
    "joule": "J",
    "K": "K",
    "kelvin": "K",
    "kg": "kg",
    "kilogram": "kg",
    "m": "m",
    "meter": "m",
    "metre": "m",
    "mol": "mol",
    "mole": "mol",
    "Pa": "Pa",
    "pascal": "Pa",
    "s": "s",
    "second": "s",
    "W": "W",
    "watt": "W",
    "hour": "h",
    "hours": "h",
    "year": "yr",
    "years": "yr",
    "yr": "yr",
}

_siunitx_prefixes = {
    "kilo": "k",
    "mega": "M",
    "giga": "G",
    "milli": "m",
    "micro": r"\mu",
    "nano": "n",
}


# Read the text of a brace-enclosed LaTeX macro, allowing for nested
# braces:
def _read_braced_argument(text, position):
    """Return the braced argument at position and the following position."""
    if position == len(text) or text[position] != "{":
        return None, position

    depth = 1
    end = position + 1
    while end < len(text) and depth:
        if text[end] == "{":
            depth += 1
        elif text[end] == "}":
            depth -= 1
        end += 1

    if depth:
        return None, position

    return text[position + 1:end - 1], end


def _format_siunitx_unit(unit):
    """Convert the subset of siunitx unit syntax used in this manual to LaTeX."""
    components = []
    denominator = False
    prefix = ""
    next_power = 1
    tokens = re.findall(r"\\[A-Za-z]+|[A-Za-z]+|[0-9]+|[*/^]", unit)
    position = 0

    def set_last_power(power):
        if components:
            components[-1][1] = -power if components[-1][1] < 0 else power

    while position < len(tokens):
        token = tokens[position]
        name = token[1:] if token.startswith("\\") else token

        if name == "per" or token == "/":
            denominator = True
        elif name == "squared":
            set_last_power(2)
        elif name == "cubed":
            set_last_power(3)
        elif name == "cubic":
            next_power = 3
        elif name in _siunitx_prefixes:
            prefix += _siunitx_prefixes[name]
        elif token == "^" and position + 1 < len(tokens):
            position += 1
            if tokens[position].isdigit():
                set_last_power(int(tokens[position]))
        elif token != "*":
            symbol = _siunitx_symbols.get(name, name)
            components.append([prefix + symbol,
                               -next_power if denominator else next_power])
            prefix = ""
            next_power = 1

        position += 1

    formatted = []
    for symbol, power in components:
        if symbol == r"^\circ" :
            factor = symbol
        else :
            factor = r"\mathrm{" + symbol + "}"
        if power != 1:
            factor += "^{" + str(power) + "}"
        formatted.append(factor)

    return r"\,".join(formatted)


def _is_inside_math(text, position):
    """Determine whether position occurs between matching dollar delimiters."""
    delimiter = None
    for match in re.finditer(r"(?<!\\)\${1,2}", text[:position]):
        value = match.group()
        if delimiter == value:
            delimiter = None
        elif delimiter is None:
            delimiter = value
    return delimiter is not None


def _convert_siunitx_commands(app, docname, source):
    """Replace legacy siunitx commands with inline LaTeX math."""
    text = source[0]
    result = []
    position = 0

    while position < len(text):
        match = re.search(r"\\(si|SI)\s*", text[position:])
        if match is None:
            result.append(text[position:])
            break

        start = position + match.start()
        command_end = position + match.end()
        first_argument, argument_end = _read_braced_argument(text, command_end)
        if first_argument is None:
            result.append(text[position:command_end])
            position = command_end
            continue

        if match.group(1) == "SI":
            while argument_end < len(text) and text[argument_end].isspace():
                argument_end += 1
            second_argument, end = _read_braced_argument(text, argument_end)
            if second_argument is None:
                result.append(text[position:argument_end])
                position = argument_end
                continue
            formula = first_argument + r"\," + _format_siunitx_unit(second_argument)
        else:
            end = argument_end
            formula = _format_siunitx_unit(first_argument)

        result.append(text[position:start])
        result.append(formula if _is_inside_math(text, start) else "$" + formula + "$")
        position = end

    source[0] = "".join(result)


def setup(app):
    app.connect("source-read", _convert_siunitx_commands)


tikz_proc_suite = "GhostScript"
# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = [
    "_build",
]

# -- Options for sphinx_tags -------------------------------------------------

# turn it on
tags_create_tags = True

# you're using MyST, so scan md files
tags_extension = ["md"]

# (optional) put generated pages here; default is "_tags"
tags_output_dir = "_tags"

# (optional) nicer titles / badge styling (needs sphinx-design)
tags_overview_title = "Page index"
tags_intro_text = "Tags:"
tags_page_title = "Tag"
tags_page_header = "Pages"
tags_index_head = "All page tags"
tags_create_badges = True
# See the color descriptions here: https://sphinx-tags.readthedocs.io/en/latest/configuration.html#badge-colors
tags_badge_colors = {"category:*": "primary", "feature:*": "success", "*":"info"}


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

# -- Options for LaTex output -------------------------------------------------

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
