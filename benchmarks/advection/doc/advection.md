(sec:benchmarks:advection)=
# Advection stabilization benchmarks

The underlying PDEs of the temperature and compositional field are typically
advection-dominated and as such, require a stabilization scheme, see
{ref}`sec:advection-stabilization` for an introduction for the methods
implemented in ASPECT.

We have several benchmarks to test the robustness, quality of solutions (size
of overshoots, smearing of sharp interfaces). Here, we give a short summary of
the benchmarks implemented:

-   Dropping box (`benchmarks/drop_*.prm`): This is a simple 2d box with a
    prescribed, constant, vertical velocity. An initial condition creates a
    square box with a high temperature, which is advected vertically.
    See {numref}`fig:benchmark-drop`.

-   Rotating Shapes: `benchmarks/rotate_shape_*.prm`: A collection of shapes
    in a 2d box rotated by 360 degrees by a prescribed velocity.
    See {numref}`fig:benchmark-rotate-shape`.

Both benchmarks have the identical setup in the temperature and a
compositional field. The only difference is that the temperature equation
contains a (small) physical diffusion term.

**[Description of benchmark files](../README.md)**

```{figure-md} fig:benchmark-drop
<img src="drop.png" />

 Dropping box benchmark at final time. Left: entropy viscosity. Right: SUPG.
```

```{figure-md} fig:benchmark-rotate-shape
<img src="rotate_shape.png" />

 Rotating shapes benchmark at final time: Left: reference. Middle: Entropy viscosity. Right: SUPG.
```

:::{toctree}
../README.md
:::
