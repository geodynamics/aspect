(part:user_manual:chap:basic_starter_tutorial:sec:BST_10_adding_basic_subducting_plate)=
Adding a basic subducting plate
============================

Now that we have a solid understanding of how area features like `oceanic plate`, `continental plate` and `mantle layer` work, we can start to take a look at [line features](part:user_manual:chap:concepts:sec:line_features). In this case we are going to take a look at the `subducting plate`, but note that there is also the `fault` feature. On a fundamental level, they both function in a very similar way. The main difference is that the slab has a slab surface and a up and down direction from that surface while the fault has a center surface and is symmetric on both sides. The temperature and compositional models available for each may also differ.

There are five properties these line features objects require:

1. `model`: The name of the model, e.g., `subducting plate` or `fault`
2. `name`: Like with area features, a descriptive name to improve readability of the file.
3. `coordinates`: The location of where the feature intersects with the surface. For a subducting plate, this would be the location of the trench. This is what the line in line feature refers. These coordinates represent a line in a map view.
4. `dip point`: This is often the most confusing concept to new users but it is really quite simple. If the slab or fault makes an angle, we need to decide which side of the line (or trench for a slab) the slab or fault is dipping. There are many ways of doing this but in GWB it is done through simply defining a point on one side of the line. The slab or fault will then dip in that direction. It is often good to choose a point far away from the trench especially if the trench makes a lot of curves.
5. `segments`: Defines the downwards part of the subducting plate (i.e., the slab) or fault. Each segment requires at least a `length` as a number, a `thickness` which is an array of one or two numbers and an `angle` which is also an array of one or two numbers. We will first show an example to get our first slab in the model before explaining further. It will hopefully be easier to understand when you see it in action in this and the next few sections.

For the example in this section, we will make a trench at the interface of two area features where the slab is dipping with a constant angle of 60 degrees in the direction of the overriding plate. We make this slab (segment) 300 km long and give it a constant thickness of 100 km. We set a dip point at a location on the side of the overriding plate with respect to the trench. In this case the origin (0,0) will work fine and the location is shown as the turquoise sphere in the figure below. 

```{note}
The result will not change based on the exact location of the dip point, as long as it is on one side of all the lines formed by connecting the line features coordinate points.
```

::::::{tab-set}

:::::{tab-item} Important lines
:sync: Partial

```{literalinclude} ../../_static/gwb_input_files/BST_10_adding_basic_subducting_plate.wb
:language: json
:lineno-start: 43
:lines: 43-49
:emphasize-lines: 2,4
```
::::{grid} 3
:::{grid-item-card} BST_10_adding_basic_subducting_plate.wb
:link: ../../_static/gwb_input_files/BST_10_adding_basic_subducting_plate.wb
:::
:::{grid-item-card} BST_10_adding_basic_subducting_plate.grid
:link: ../../_static/gwb_input_files/BST_10_adding_basic_subducting_plate.grid
:::
:::{grid-item-card} Paraview v3 state file 
:link: ../../_static/paraview_state_files/BST_v3.pvsm
:::
::::
:::::

:::::{tab-item} Full file
:sync: Full


```{literalinclude} ../../_static/gwb_input_files/BST_10_adding_basic_subducting_plate.wb
:language: json
:lineno-start: 1
:emphasize-lines: 44,46
```

::::{grid} 3
:::{grid-item-card} BST_10_adding_basic_subducting_plate.wb
:link: ../../_static/gwb_input_files/BST_10_adding_basic_subducting_plate.wb
:::
:::{grid-item-card} BST_10_adding_basic_subducting_plate.grid
:link: ../../_static/gwb_input_files/BST_10_adding_basic_subducting_plate.grid
:::
:::{grid-item-card} Paraview v3 state file 
:link: ../../_static/paraview_state_files/BST_v3.pvsm
:::
::::
:::::

::::::


```{todo}
The explanation can be significantly improved by adding conceptual figures
```
```{figure} ../../../../doc/sphinx/_static/images/user_manual/basic_starter_tutorial/BST_10.png
:name: BST_10
:alt: Basic Starter Tutorial section 10. 
:align: center

Basic Starter Tutorial section 10. The top part of the figure shows where the composition has been assigned as an object. Currently it shows composition 0 as green, composition 1 as yellow, composition 2 as purple and composition 3 as blue. Composition 4 is not shown to be able to see the slab. The front half of the overriding plate (composition 1) has also been removed to be able to better view the slab. The turquoise sphere is the location of the dip point. The bottom part shows the temperature as seen slightly from below where only temperatures between 300K and 1600K are shown. 
```
