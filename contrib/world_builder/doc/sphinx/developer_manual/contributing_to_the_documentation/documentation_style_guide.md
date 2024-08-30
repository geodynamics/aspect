(part:dev_manual:chap:contribute_to_doc:sec:doc_style_guide)=
Documentation style guide
=========================

Not only are there many ways to style the documents, there are also many ways with Sphinx to reach the same output. To attempt to unify the style of the documentation and the documentation code, this page contains some general styling recommendations.

# Page titles

Page titles should always use equal signs below it (`=`) to create the heading. This creates a clear visual distinction between normal headers and the page title.

# Page label

The page label should always be a path, separated by colons (`:`) in the following style: `part:dev_manual:chap:contribute_to_doc:sec:doc_style_guide`.

# Including advanced elements

When possible, the use of colons (`:`) should be preferred over the use of grave accents (\`) when creating advanced elements such as todo notes, figures, code blocks or tables. This means that:

::::{code-block} md
:::{note}
My Note
:::
::::

is preferred over the use of 

::::{code-block} md
```{note}
My Note
````
::::

# Including World Builder files

When including a World Builder file, it is important to be able to focus on the parts which are important, while also enabling the reader to see the surrounding context and the whole file. Furthermore, the user should have the capability to run the models themselves and get the same output as shown. To standardize this the following structure is strongly encouraged. 

1. The World Builder file and a corresponding grid file need to be stored separately, either in the cookbook, or in the _static directory.
2. The use of a tab-set with a tab with important lines and a tab with the full World Builder file, highlighting important lines.
3. Adding at least a link to the World Builder file and a grid file to visualize it, preferably also a link to a `.pvtu` file which allows to get the same output as the figure, if present.
4. Add a figure of the output.

Below is an example code block.

:::::::{code-block} md

::::::{tab-set}

:::::{tab-item} Important lines
:sync: Partial

```{literalinclude} ../../_static/gwb_input_files/BST_05_limit_depth.wb
:language: json
:lineno-start: 5
:lines: 5-15
:emphasize-lines: 3,8
```
::::{grid} 3
:::{grid-item-card} BST_05_limit_depth.wb
:link: ../../_static/gwb_input_files/BST_05_limit_depth.wb
:::
:::{grid-item-card} BST_05_limit_depth.grid
:link: ../../_static/gwb_input_files/BST_05_limit_depth.grid
:::
:::{grid-item-card} Paraview v1 state file 
:link: ../../_static/paraview_state_files/BST_v1.pvsm
:::
::::
:::::

:::::{tab-item} Full file
:sync: Full


```{literalinclude} ../../_static/gwb_input_files/BST_05_limit_depth.wb
:language: json
:lineno-start: 1
:emphasize-lines: 7,12
```

::::{grid} 3
:::{grid-item-card} BST_05_limit_depth.wb
:link: ../../_static/gwb_input_files/BST_05_limit_depth.wb
:::
:::{grid-item-card} BST_05_limit_depth.grid
:link: ../../_static/gwb_input_files/BST_05_limit_depth.grid
:::
:::{grid-item-card} Paraview v1 state file 
:link: ../../_static/paraview_state_files/BST_v1.pvsm
:::
::::
:::::

::::::


```{figure} ../../../../doc/sphinx/_static/images/user_manual/basic_starter_tutorial/BST_05.png
:name: BST_05_example
:alt: Basic Starter Tutorial section 5. 
:align: center

Image description
```
:::::::


which looks like this when rendered:


::::::{tab-set}

:::::{tab-item} Important lines
:sync: Partial

```{literalinclude} ../../_static/gwb_input_files/BST_05_limit_depth.wb
:language: json
:lineno-start: 5
:lines: 5-15
:emphasize-lines: 3,8
```
::::{grid} 3
:::{grid-item-card} BST_05_limit_depth.wb
:link: ../../_static/gwb_input_files/BST_05_limit_depth.wb
:::
:::{grid-item-card} BST_05_limit_depth.grid
:link: ../../_static/gwb_input_files/BST_05_limit_depth.grid
:::
:::{grid-item-card} Paraview v1 state file 
:link: ../../_static/paraview_state_files/BST_v1.pvsm
:::
::::
:::::

:::::{tab-item} Full file
:sync: Full


```{literalinclude} ../../_static/gwb_input_files/BST_05_limit_depth.wb
:language: json
:lineno-start: 1
:emphasize-lines: 7,12
```

::::{grid} 3
:::{grid-item-card} BST_05_limit_depth.wb
:link: ../../_static/gwb_input_files/BST_05_limit_depth.wb
:::
:::{grid-item-card} BST_05_limit_depth.grid
:link: ../../_static/gwb_input_files/BST_05_limit_depth.grid
:::
:::{grid-item-card} Paraview v1 state file 
:link: ../../_static/paraview_state_files/BST_v1.pvsm
:::
::::
:::::

::::::


```{figure} ../../../../doc/sphinx/_static/images/user_manual/basic_starter_tutorial/BST_05.png
:name: BST_05_example
:alt: Basic Starter Tutorial section 5. 
:align: center

Image descriptionta
```

# Units

The easiest way of adding superscripts needed for units such as kg/m<sup>3</sup> is to use HTML tags, as in `<sup>superscript</sup>`.
