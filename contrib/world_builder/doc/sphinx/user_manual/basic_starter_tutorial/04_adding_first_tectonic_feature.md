(part:user_manual:chap:basic_starter_tutorial:sec:04_first_tectonic_feature)=
Adding your first tectonic feature
========================

Now we are finally ready to get out our brushes and start coloring in the world. We are going to start by adding an oceanic plate to the feature list. 

Each feature is an object and is enclosed in curly braces (`{}`). You will need to specify the following keys:

1. The `model` key always defines which feature is being set. In this case it is "oceanic plate". 
2. The `name` key should contain a descriptive name. The name can be anything and is used to help keep the input file readable. Here we name it "Overriding Plate".
3. The `coordinates` key is a list of 2D coordinates: x,y in Cartesian and long,lat in spherical. A list in JSON is indicated by square brackets (`[]`) and the items are separated by commas.
4. (Optional) The `temperature models` key is a list of temperature models. Each temperature model is an object. For now we will just choose the simplest one: `uniform` where we set the temperature to `293`K.
5. (Optional) The `compositional models` key is a list of compositional models. Each compositional model is an object. Like with the temperature model, we will just chose the simplest one: `uniform` and set the composition in compositional field to 0. We will go into more details about how this works later in the tutorial.

There are some other options that will be covered later. The file below shows the result. To focus your attention, by default only the lines of interest are shown but you can always view the full file by clicking on the `Full file` tab.


::::::{tab-set}

:::::{tab-item} Important lines
:sync: Partial

```{literalinclude} ../../_static/gwb_input_files/BST_04_overriding_plate.wb
:language: json
:lineno-start: 4
:lines: 4-12
```
::::{grid} 3
:::{grid-item-card} BST_04_overriding_plate.wb
:link: ../../_static/gwb_input_files/BST_04_overriding_plate.wb
:::
:::{grid-item-card} BST_04_overriding_plate.grid
:link: ../../_static/gwb_input_files/BST_04_overriding_plate.grid
:::
:::{grid-item-card} Paraview v1 state file 
:link: ../../_static/paraview_state_files/BST_v1.pvsm
:::
::::
:::::

:::::{tab-item} Full file
:sync: Full


```{literalinclude} ../../_static/gwb_input_files/BST_04_overriding_plate.wb
:language: json
:lineno-start: 1
```

::::{grid} 3
:::{grid-item-card} BST_04_overriding_plate.wb
:link: ../../_static/gwb_input_files/BST_04_overriding_plate.wb
:::
:::{grid-item-card} BST_04_overriding_plate.grid
:link: ../../_static/gwb_input_files/BST_04_overriding_plate.grid
:::
:::{grid-item-card} Paraview v1 state file 
:link: ../../_static/paraview_state_files/BST_v1.pvsm
:::
::::
:::::

::::::


```{figure} ../../../../doc/sphinx/_static/images/user_manual/basic_starter_tutorial/BST_04.png
:name: BST_04
:alt: Basic Starter Tutorial section 16 highres result. 
:align: center

Basic Starter Tutorial section 04. The top part of the figure shows where the composition has been assigned as an object (green box). Here, the composition has been assigned the value 0. The bottom figure shows the temperature as viewed slightly from below. This shows that the overriding oceanic plate has a temperature of 293K that extends from the surface to the bottom of the model.
```
