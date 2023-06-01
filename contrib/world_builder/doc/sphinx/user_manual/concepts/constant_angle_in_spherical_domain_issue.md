(part:user_manual:chap:concepts:sec:const_angle_spherical)=
Constant angle in spherical domain issue
========================================

There is one last concept we need to discuss before we can move on to the {ref}`part:user_manual:chap:basic_starter_tutorial:sec:index`, which is the problem of defining what a constant angle in a sphere is. This does not matter for Cartesian models, but is important fore line features in spherical models.

The issue arises when you want to define a line with a certain angle to the surface in a sphere. The problem here is that there are two ways of looking at this. The first way is to see a constant angle with the point it starts at. This works fine, except that for large distances it will burst through the surface of the planet out into space. This can be see in the left figure below. The other option is to look at a constant angle as an angle which remains constant with the surface above. This results in a logarithmic spiral as seen in the right figure below.

::::{grid} 2

:::{grid-item-card}

**Constant angle with the starting position**

![spherical subduction line](../../../manual/images/spherical_approximations_line.png)
:::
:::{grid-item-card}

**Constant angle with the surface above the current position**

![spherical subduction spiral](../../../manual/images/spherical_approximations_spiral.png)
:::
::::

Currently only the left figure, and something in between the left and right figure has been implemented in the GWB. The in between option is that for each line segment it can adjust the angle to be the correct angle with respect to the surface. Implementing the way it is done in the right figure, is possible if there is enough interest. 

In the GWB this is defined in the spherical coordinate system through something called a `depth method`. The left figure method is called `starting point`, because the angle is set and kept relative to the starting point. The in between option is called `begin segment`, because the angle is relative to the beginning of each segment. The right figure method is called `contiuous`, and is not implemented.
