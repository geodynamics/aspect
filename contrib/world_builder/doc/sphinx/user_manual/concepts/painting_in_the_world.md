(part:user_manual:chap:concepts:sec:painting_the_world)=
Painting in the world
=====================

In the last section we introduced the concept of a feature. In this section we are going to discuss how features are placed inside the world and how they can interact with each other. 

The first thing to understand about how GWB features are primarily defined by the user as 2D features on a map. That means that the user provides 2D coordinates of where the feature is located. Only after that, the user provides information to actually make the feature a 3D object. The second thing to understand is that features are provided in a list. The GWB goes through that list from top to bottom and adds each feature to the world in order. If there is any overlap, the user can decide what to do (how to blend the the old and new value), but by default each feature overwrites the old value with the new value. You can think of it like paining, where each features add a layer of paint and the painter can decide to completely paint over an area, or create an interesting blend!

Starting with a topographic map
-------------------------------

In the previous section, we referred to that a user should be able to take a map of a region and start to make a model. That is exactly what we are going to do here.

Below you can see three pictures. In the leftmost one we are starting with a map of a region we are interested in making a model of. For this example we have made a simplified map of the eastern Caribbean plate. For simplicity's sake, we are only going to focus on a small part, as indicated by the thick black line in the middle figure. Next we are going to define important points on that map, as indicated by P1 to P7. The way we define these points is through looking at the different tectonic units in the area. P1, P2, P7, P5 and P6 define the area to the left of the trench, which represents the Caribbean plate, and is indicated by the green area in the right most figure. P2, P3, P4, P5 and P7 define the area to the right of the trench, which represents the Atlantic part of the South American plate, and is indicated by the blue color in the right most figure. Finally, P2, P7 and P5, indicate the trench and is indicated by the purple color in the right most figure.


::::{grid} 3

:::{grid-item-card}

**Starting with a "geologic map"**

![Starting with a geologic map](../../_static/images/user_manual/map_to_top_view_v2_plain.svg)
:::
:::{grid-item-card}

**Now we select a region and define points**

![Starting with a geologic map](../../_static/images/user_manual/map_to_top_view_v2_plain_text_frame_sds.svg)
:::
:::{grid-item-card}

**Coloring in the rest of the region**

![Starting with a geologic map](../../_static/images/user_manual/map_to_top_view_v2_plain_text_frame_sds_orp_sdp.svg)
:::
::::



Using the topographic map to make a 3D model
--------------------------------------------

Now that we have colored in our map, it is time to talk about how this map helps us making a 3D setup. To do that we are going to use the figures below. The right side of each figure shows a map view as shown in the previous subsection. The right part of the figure show a 3D rendering of the setup. We are going to call the colors we are adding compositions. You can think of this different types of rock.

:::{card}

**We start with a box without compositions**

![Starting with a box without compositions](../../_static/images/user_manual/gwb_box_building_plain_text_frame.svg)

:::

To show how we stack features, we are first we are going to add a composition for the mantle. That is we are going to use a feature called a `mantle layer`. We give the feature the coordinates P1, P2, P4 and P6, and tell it to start a the surface and have a thickness which is enough to fill the whole model. In the figure below we have given the area which is now defined as the mantle composition a red color.

:::{card}

**Next we add a mantle composition**

![Starting with a box without compositions](../../_static/images/user_manual/gwb_box_building_plain_text_frame_mtl.svg)
:::

Next we are going to add the overriding Caribbean plate. We know that it is make of oceanic lithosphere, so we are going to use a feature called `oceanic plate`. Like before we provide it with the coordinates (P1, P2, P7, P5 and P6) and we give it a thickness of 90km. As you can see, a part of the area which was red before is now *painted over* by a green color, which represent a different composition.

:::{card}

**Now we add an overriding plate**


![adding an overriding plate](../../_static/images/user_manual/gwb_box_building_plain_text_frame_mtl_orp.svg)
:::

We do the same thing for the Atlantic part of the South American plate and color it blue like before.

:::{card}

**Then we add an oceanic plate**


![adding an oceanic plate](../../_static/images/user_manual/gwb_box_building_plain_text_frame_mtl_orp_sdp.svg)
:::

Now we get to the trench. The important thing about the trench is that is has subducting lithosphere, so we are going to use a feature called `subducting plate`. We provide the location of the trench to the subducting plate (P2, P7 and P5). Next we need to provide some other information like the angle of the plate with the surface and thickness. In practice the slab can have different segments with different angles and thicknesses, but for now we are just going to use a single slab segment with a constant angle and thickness. Like before, we are going to color in the area and overwrite the existing colors in those places. This provides everything we need to make a 3D model! 

:::{card}

**Finally, we add a slab, the "Painting" of the model is now complete**


![Painting in model completed by adding a slab](../../_static/images/user_manual/gwb_box_building_plain_text_frame_mtl_orp_sdp_sds.svg)
:::
