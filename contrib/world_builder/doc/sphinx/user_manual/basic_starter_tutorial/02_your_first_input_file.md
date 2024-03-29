(part:user_manual:chap:basic_starter_tutorial:sec:02_first_input_file)=
Your first input file
=====================

GWB input files are plain-text files which generally have the extension `.wb`. At minimum, the World Builder input file should contain a object (indicated by `{}` in JSON) containing a version number (JSON key: `version`) and a list (indicated by `[]` in JSON) of features (JSON key: `features`). This results in the following minimal input file:  



```{literalinclude} ../../_static/gwb_input_files/BST_02_minimal_box.wb
:language: json
:lineno-start: 1
```

::::{grid} 3
:::{grid-item-card} BST_1_minimal_box.wb
:link: ../../_static/gwb_input_files/BST_02_minimal_box.wb
:::
:::{grid-item-card} BST_1_minimal_box.grid
:link: ../../_static/gwb_input_files/BST_02_minimal_box.grid
:::
:::{grid-item-card} Paraview v1 state file 
:link: ../../_static/paraview_state_files/BST_v1.pvsm
:::
::::

```{note}
You can copy the text in the box by selecting it or by hovering over the box and clicking the copy symbol in the top right corner. The file can also be downloaded by clicking on the `BST_1_minimal_box.wb` button below the textbox. 

If you want to inspect the results yourself, you can also download the corresponding `.grid` file and use the `gwb-grid` application (`gwb-grid BST_1_minimal_box.wb BST_1_minimal_box.grid`) to create a vtk output (`BST_1_minimal_box.vtk`) and view it in, for example, Paraview.
```

Congratulations on creating your first World Builder input file! If you visualize the result (shown below), you will notice an adiabatic temperature profile and that the compositional field is zero everywhere. This is the background, or our canvas if we want to stay in the painting analogy. Now, before grabbing our brushes, we should first discuss the shape of our canvas.


```{figure} ../../../../doc/sphinx/_static/images/user_manual/basic_starter_tutorial/BST_02.png
:name: BST_02
:alt: Basic Starter Tutorial section 16 highres result. 
:align: center

Basic Starter Tutorial section 2. The top part of the figure (blank) shows where the composition has been assigned. Since no compositions have been assigned yet, nothing is plotted when we view the composition solution. The bottom part shows the temperature solution that has been assigned. This shows the adiabatic temperature gradient.
Unless otherwise stated, all figures are generated with double the resolution given as default in the above grid file. Since all tutorials are automatically tested, running at full resolution for every pull request would take too much computational time. The resolution shown in the figures is commented out under the subsection `shown grid properties`.
```
