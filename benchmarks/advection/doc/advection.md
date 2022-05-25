#### Advection stabilization benchmarks

The underlying PDEs of the temperature and compositional field are typically
advection-dominated and as such, require a stabilization scheme, see
[\[sec:advection-stabilization\]][1] for an introduction for the methods
implemented in .

We have several benchmarks to test the robustness, quality of solutions (size
of overshoots, smearing of sharp interfaces). Here, we give a short summary of
the benchmarks implemented:

-   Dropping box (`benchmarks/drop_*.prm`): This is a simple 2d box with a
    prescribed, constant, vertical velocity. An initial condition creates a
    square box with a high temperature, which is advected vertically. See
    Figure&nbsp;[1][].

-   Rotating Shapes: `benchmarks/rotate_shape_*.prm`: A collection of shapes
    in a 2d box rotated by 360 degrees by a prescribed velocity. See
    Figure&nbsp;[2][].

Both benchmarks have the identical setup in the temperature and a
compositional field. The only difference is that the temperature equation
contains a (small) physical diffusion term.

```{figure-md} fig:benchmark-drop
<img src="drop.png" alt="Figure" />

Dropping box benchmark at final time. Left: entropy viscosity. Right: SUPG.
```

```{figure-md} fig:benchmark-rotate-shape
<img src="rotate_shape.png" alt="Figure" />

Rotating shapes benchmark at final time: Left: reference. Middle: Entropy viscosity. Right: SUPG.
```

  [1]: #sec:advection-stabilization
  [1]: #fig:benchmark-drop
  [2]: #fig:benchmark-rotate-shape
