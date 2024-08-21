(part:user_manual:chap:basic_starter_tutorial:sec:12_subducting_plate_temp)=
Subducting plate temperatures
=============================

Before exploring the segments more, let's briefly discuss temperature models for the subducting plates. The temperature structure of a slab can be quite complicated. 

A model which was often used as a temperature structure is from or based on the McKenzie temperature structure {cite:p}`McKenzie_1970`. In the World Builder this is implemented as the subducting plate temperature model called the `plate model` ({cite:t}`turcotte_schubert_2014`). If you want to use this, you will need to provide a reference density and a plate velocity like: 
```{code-block} json
---
lineno-start: 1
---
"temperature models":[
  {"model":"plate model", "density":3300, "plate velocity":0.02}
]
```

This is still perfectly valid and usable but the approach has some fundamental issues including that the temperature structure set-up is not mass conservative. This is why the `mass conserving` slab temperature model was developed (Todo: insert reference when paper is done). This temperature model also brings a lot of other large and small improvements which we can't go over in this starter tutorial. Although this temperature is generally recommended over the older McKenzie model, it is also is a bit more complicated in usage. So to keep this tutorial simple, we will use the `plate model` here, but we strongly recommend looking into the `mass conserving` temperature model for subduction models.

```{todo}
Create a cookbook for the mass conserving temperature structure and point to it from here.
```


::::::{tab-set}

:::::{tab-item} Important lines
:sync: Partial

```{literalinclude} ../../_static/gwb_input_files/BST_12_subducting_plate_temperatures.wb
:language: json
:lineno-start: 43
:lines: 43-50
:emphasize-lines: 5
```
::::{grid} 3
:::{grid-item-card} BST_12_subducting_plate_temperatures.wb
:link: ../../_static/gwb_input_files/BST_12_subducting_plate_temperatures.wb
:::
:::{grid-item-card} BST_12_subducting_plate_temperatures.grid
:link: ../../_static/gwb_input_files/BST_12_subducting_plate_temperatures.grid
:::
:::{grid-item-card} Paraview v3 state file 
:link: ../../_static/paraview_state_files/BST_v3.pvsm
:::
::::
:::::

:::::{tab-item} Full file
:sync: Full


```{literalinclude} ../../_static/gwb_input_files/BST_12_subducting_plate_temperatures.wb
:language: json
:lineno-start: 1
:emphasize-lines: 47
```

::::{grid} 3
:::{grid-item-card} BST_12_subducting_plate_temperatures.wb
:link: ../../_static/gwb_input_files/BST_12_subducting_plate_temperatures.wb
:::
:::{grid-item-card} BST_12_subducting_plate_temperatures.grid
:link: ../../_static/gwb_input_files/BST_12_subducting_plate_temperatures.grid
:::
:::{grid-item-card} Paraview v3 state file 
:link: ../../_static/paraview_state_files/BST_v3.pvsm
:::
::::
:::::

::::::


```{figure} ../../../../doc/sphinx/_static/images/user_manual/basic_starter_tutorial/BST_12.png
:name: BST_12
:alt: Basic Starter Tutorial section 12. 
:align: center

Basic Starter Tutorial section 12. The top part of the figure shows where the composition has been assigned as an object. Currently it shows composition 0 as green, composition 1 as yellow, composition 2 as purple and composition 3 as blue. Composition 4 is not shown to be able to see the slab. The front half of the overriding plate (composition 1) has also been removed to be able to better view the slab. The bottom part shows the temperature as seen slightly from below where only temperatures between 300K and 1600K are shown. 
```
