(part:dev_manual:chap:contribute_to_doc:sec:important_syntax)=
Important syntax
===============

The World Builder is using [sphinx](https://www.sphinx-doc.org) to generate the documentation. Sphinx accepts multiple variants of extended markdown, but in this project we are using [myst](https://myst-parser.readthedocs.io). Below we list some common commands to use.

:::{note}
Everything in this file will look like garbage in a regular markdown viewer, like if you’re viewing this on GitHub. Viewing it on readthedocs will render everything properly.
:::

(title)=
# Title
Every page starts with a header for this file:
::::{code-block} md
(part:dev_manual:chap:contribute_to_doc:sec:important_syntax)=
Important syntax
=====================================

::::

Note that the first line contains a reference to this file that can be used to link to this file from other places of the documentation.

(text)=
# Text

Put single asterisks (*) around text you *want to italicize*, but double asterisks around text you **want to bold**. If you want to `highlight code` you put single Grave Accent (\`) around it. You can escape chacters with the backlash (\), like so \`not highlighted\`, \*not italic\*, \*\*not bold\*\*.

:::{code-block} md
Put single asterisks (*) around text you *want to italicize*, 
but double asterisks around text you **want to bold**. 
If you want to `highlight code` you put single Grave Accent (\`) around it.
You can escape chacters with the backlash (\\), like so \`not highlighted\`, \*not italic\*, \*\*not bold\*\*.
:::

(headers)=
# Headers

Headers start with a # for level 1, ## for level 2, etc. 

(header2)=
## Header 2

(header3)=
### Header 3

:::{code-block} md
(header1)=
# Header 1

(header2)=
## Header 2

(header3)=
### Header 3
:::


# Links
Links are written in the following way:
:::{code-block} md
[link within this page](title)
:::
[link within this page](title)

:::{code-block} md
[link to other page](part:dev_manual:chap:contribute_to_doc:sec:index)
:::

[link to other page](part:dev_manual:chap:contribute_to_doc:sec:index)

:::{code-block} md
[link to external website](https://github.com/GeodynamicWorldBuilder/WorldBuilder)
:::
[link to external website](https://github.com/GeodynamicWorldBuilder/WorldBuilder)

# Admonitions

Adding notes, warnings or todos boxes can be done with the following syntax:

::::{code-block} md
:::{todo}
This is a `TODO` admonition.
:::
::::
:::{todo}
This is a `TODO` admonition.
:::

::::{code-block} md
:::{note}
This is a `note` admonition.
:::
::::
:::{note}
This is a `note` admonition.
:::

::::{code-block} md
:::{tip}
This is a `tip` admonition.
:::
::::
:::{tip}
This is a `tip` admonition.
:::

::::{code-block} md
:::{seealso}
This is a `seealso` admonition.
:::
::::
:::{seealso}
This is a `seealso` admonition.
:::

::::{code-block} md
:::{important}
This is an `important` admonition.
:::
::::
:::{important}
This is an `important` admonition.
:::

::::{code-block} md
:::{warning}
This is a `warning` admonition.
:::
::::
:::{warning}
This is a `warning` admonition.
:::

::::{code-block} md
:::{attention}
This is an `attention` admonition.
:::
::::
:::{attention}
This is an `attention` admonition.
:::

::::{code-block} md
:::{danger}
This is a `danger` admonition.
:::
::::
:::{danger}
This is a `danger` admonition.
:::

::::{code-block} md
:::{error}
This is an `error` admonition.
:::
::::
:::{error}
This is an `error` admonition.
:::

# Figures

Adding figures looks similar to notes and todo's:

::::{code-block} md
:::{figure} ../../../../doc/sphinx/_static/images/user_manual/basic_starter_tutorial/BST_19.png
:name: BST_19_syntax
:alt: Basic Starter Tutorial section 19 highres result. 
:align: center

Basic Starter Tutorial section 19 high resolution result. This has 8 times the resolution than the grid file above. Note that you can see the Earth's curvature! 
:::
::::

:::{figure} ../../../../doc/sphinx/_static/images/user_manual/basic_starter_tutorial/BST_19.png
:name: BST_19_syntax
:alt: Basic Starter Tutorial section 19 highres result. 
:align: center

Basic Starter Tutorial section 19 high resolution result. This has 8 times the resolution than the grid file above. Note that you can see the Earth's curvature! 
:::

# Code blocks

Creating code blocks is similar to making notes or todos. Instead of `{note}` you need to use `{code-block}` and then the name of the language which the code should be colored by. For example for json, you would use:

::::{code-block} md
:::{code-block} json
:lineno-start: 1
:emphasize-lines: 3,4
 {
    "version":1.0,
    "coordinate system": {"model":"cartesian"}, 
 }
:::
::::

The lines starting with a colon (:) specify certain parameters which can be set, such as showing line numbers and what number they need to start at (lineno-start: 1) and what lines to mark (emphasize-lines: 3,4).

:::{code-block} json
:lineno-start: 1
:emphasize-lines: 3,4
 {
    "version":1.0,
    "coordinate system": {"model":"cartesian"}, 
    "Features":[]
 }
:::

If you need to wrap a code block in a code block, you just use more : on the outside:

:::::{code-block} md
::::{code-block} md
:::{code-block} json
:lineno-start: 1
:emphasize-lines: 3,4
 {
    "version":1.0,
    "coordinate system": {"model":"cartesian"}, 
 }
:::
::::
:::::

Including code from files can be done in the following way:
::::{code-block} md
:::{literalinclude} ../../_static/gwb_input_files/BST_02_minimal_box.wb
:language: json
:lineno-start: 1
:::
::::

:::{literalinclude} ../../_static/gwb_input_files/BST_02_minimal_box.wb
:language: json
:lineno-start: 1
:::

# Dropdown
:::::{code-block}

::::{dropdown} Click me to open a dropdown box
:name: this:is:my:label

This text is only shown when you click on the dropdown button.
::::
:::::

::::{dropdown} Click me to open a dropdown box
:name: this:is:my:label

This text is only shown when you click on the dropdown button.
::::

# Grids

Girds work through multiple layers, like a code-block in a code block:

:::::{code-block}
::::{grid} 3
:::{grid-item-card} BST_19_spherical_models.wb
:link: ../../_static/gwb_input_files/BST_19_spherical_models.wb
:::
:::{grid-item-card} BST_19_spherical_models.grid
:link: ../../_static/gwb_input_files/BST_19_spherical_models.grid
:::
:::{grid-item-card} Paraview Spherical state file 
:link: ../../_static/paraview_state_files/BST_spherical.pvsm
:::
::::
:::::

::::{grid} 3
:::{grid-item-card} BST_19_spherical_models.wb
:link: ../../_static/gwb_input_files/BST_19_spherical_models.wb
:::
:::{grid-item-card} BST_19_spherical_models.grid
:link: ../../_static/gwb_input_files/BST_19_spherical_models.grid
:::
:::{grid-item-card} Paraview Spherical state file 
:link: ../../_static/paraview_state_files/BST_spherical.pvsm
:::
::::


# Tabs

:::::{code-block}

::::{tab-set}
:::{tab-item} tab 1
:sync: tab1
 tab 1
:::

:::{tab-item} tab 2
:sync: tab2

tab 2
:::
::::
:::::

::::{tab-set}
:::{tab-item} tab 1
:sync: tab1
 tab 1
:::

:::{tab-item} tab 2
:sync: tab2

tab 2
:::
::::

# Citations
:::{code-block} md
Traditional citation {cite}`Fraters_Thieulot_etal_2019`.
Citation as noun {cite:t}`Fraters_Thieulot_etal_2019`.
Use a comma for multiple, like {cite:t}`Fraters_Thieulot_etal_2019,Kronbichler_Heister_etal_2012,Heister_Dannberg_etal_2017`.
:::

becomes

Traditional citation {cite}`Fraters_Thieulot_etal_2019`.
Citation as noun {cite:t}`Fraters_Thieulot_etal_2019`.
Use a comma for multiple, like ({cite:t}`Fraters_Thieulot_etal_2019,Kronbichler_Heister_etal_2012,Heister_Dannberg_etal_2017`.

If the paper you wish to cite is not already in [references.bib](https://github.com/geodynamics/aspect/blob/main/doc/sphinx/bibliography.bib), you will need to add its bibtex entry. Please format the citation tag as:

- `firstauthor_year` for one author
- `firstauthor+secondauthor_year` for two authors
- `firstauthor+secondauthor_year` for three authors
- `firstauthor+secondauthor_etal_year` for more than two authors

where `firstauthor` and `secondauthor` refer to the *last* names of the authors.

# Math

Inline math works just like in latex. Use dollar signs to put stuff like $\eta_r(z)$ in a sentence. That's:

:::{code-block} md
Use dollar signs to put stuff like $\eta_r(z)$ in a sentence.
:::

Be aware that not all latex packages are known by Sphinx, so some commands will need to be reformatted (i.e., \num{4e5} needs to be rewritten as 4\times 10^{5}.)

For math blocks:

::::{code-block} md
:::{math}
:label: eqn:one:two
F = m a
:::
::::

becomes

:::{math}
:label: eqn:one:two
F = m a
:::
::::{code-block} md
Refer to the labeled equation like {math:numref}`eqn:one:two`.
::::
Refer to the labeled equation like {math:numref}`eqn:one:two`.
(The label is not necessary if the equation doesn't have one, in which case it would be formatted like):
::::md
:::{math}
F = ma
:::
::::

The split environment is built into the math directive; which allows you to use one align operator (`&`).

::::{code-block} md
:::{math}
:label: eq:aligned
  -\nabla \cdot \left[2\eta
   \left(\varepsilon(\mathbf u) - \frac{1}{3}(\nabla \cdot \mathbf u)\mathbf 1\right)\right]
   + \nabla p' &=  -\bar \alpha \bar\rho T' \mathbf g \qquad \textrm{in $\Omega$},  \\
  \nabla \cdot (\bar\rho \mathbf u) &= 0  \qquad  \textrm{in $\Omega$}.
:::
::::

becomes:
:::{math}
:label: eq:aligned
  -\nabla \cdot \left[2\eta
   \left(\varepsilon(\mathbf u) - \frac{1}{3}(\nabla \cdot \mathbf u)\mathbf 1\right)\right]
   + \nabla p' &=  -\bar \alpha \bar\rho T' \mathbf g \qquad \textrm{in $\Omega$},  \\
  \nabla \cdot (\bar\rho \mathbf u) &= 0  \qquad \textrm{in $\Omega$}.
:::

If you require more than one alignment character you need to start and end an additional `aligned` environment. Be sure to add the same number of align operators (`&`) to each line and use `\\` to denote a new line.

::::{code-block} md
:::{math}
:label: eq:aligned
\begin{aligned}
  -\nabla \cdot \left[2\eta
   \left(\varepsilon(\mathbf u) - \frac{1}{3}(\nabla \cdot \mathbf u)\mathbf 1\right)\right]
   + \nabla p' &=  -\bar \alpha \bar\rho T' \mathbf g & \qquad  & \textrm{in $\Omega$},  \\
  \nabla \cdot (\bar\rho \mathbf u) &= 0  & \qquad  & \textrm{in $\Omega$}.
\end{aligned}
:::
::::

becomes:

:::{math}
:label: eq:multi-aligned
\begin{aligned}
  -\nabla \cdot \left[2\eta
   \left(\varepsilon(\mathbf u) - \frac{1}{3}(\nabla \cdot \mathbf u)\mathbf 1\right)\right]
   + \nabla p' &=  -\bar \alpha \bar\rho T' \mathbf g & \qquad  & \textrm{in $\Omega$},  \\
  \nabla \cdot (\bar\rho \mathbf u) &= 0  & \qquad  & \textrm{in $\Omega$}.
\end{aligned}
:::

:::{todo}
Replace equations with more appropriate equations for the World Builder. Maybe the MKenzie slab temperature equations?
:::


# Footnotes

:::{code-block} md
Here's how you put a footnote[^footnote1] in a sentence.
:::

Here's how you put a footnote[^footnote1] in a sentence (see the bottom of this file).



:::{code-block} md
[^footnote1]: And here's what the footnote says (put this at the very bottom of the document)
:::

[^footnote1]: And here's what the footnote says (put this at the very bottom of the document)



## More info

- [Myst-parser website](https://myst-parser.readthedocs.io/en/latest/syntax/syntax.html)
- [Jupyter book website cheatsheet](https://jupyterbook.org/en/stable/reference/cheatsheet.html)
- [ASPECT myst quick reference](https://aspect-documentation.readthedocs.io/en/latest/user/extending/write-a-cookbook/quickref.html)
