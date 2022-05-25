#### Reading in compositional initial composition files generated with geomIO

*This section was contributed by Juliane Dannberg*

Many geophysical setups require initial conditions with several different
materials and complex geometries. Hence, sometimes it would be easier to
generate the initial geometries of the materials as a drawing instead of by
writing code. The MATLAB-based library geomIO
(<https://bitbucket.org/geomio/geomio>, (Bauville and Baumann 2019)) provides
a convenient tool to convert a drawing generated with the vector graphics
editor Inkscape (<https://inkscape.org/en/>) to a data file that can be read
into . Here, we will demonstrate how this can be done for a 2D setup for a
model with one compositional field, but geomIO also has the capability to
create 3D volumes based on a series of 2D vector drawings using any number of
different materials. Similarly, initial conditions defined in this way can
also be used with particles instead of compositional fields.

To obtain the developer version of geomIO, you can clone the bitbucket
repository by executing the command

     git clone https://bitbucket.org/geomio/geomio.git

or you can download geomIO [here][]. You will then need to add the geomIO
source folders to your MATLAB path by running the file located in
`/path/to/geomio/installation/InstallGeomIO.m`. An extensive documentation for
how to use geomIO can be found [here][1]. Among other things, it explains [how
to generate drawings in Inkscape][] that can be read in by geomIO, which
involves assigning new attributes to paths in Inkscape&rsquo;s XML editor. In
particular, a new property &lsquo;phase&rsquo; has to be added to each path,
and set to a value corresponding to the index of the material that should be
present in this region in the initial condition of the geodynamic model.

We will here use a drawing of a jellyfish located in
[cookbooks/geomio/doc/jellyfish.svg][], where different phases have already
been assigned to each path (Figure&nbsp;[1][]).

```{figure-md} fig:jelly-picture
<embed src="cookbooks/geomio/doc/jellyfish.pdf" style="width:20.0%" />

Vector drawing of a jellyfish.
```

After geomIO is initialized in MATLAB, we [run geomIO as described in the
documentation][], loading the default options and then specifying all the
option we want to change, such as the path to the input file, or the
resolution:

``` matlab
```

You can view all of the options available by typing `opt` in MATLAB.

In the next step we create the grid that is used for the coordinates in the
`ascii data` initial conditions file and assign a phase to each grid point:

``` matlab
```

You can plot the `Phase` variable in MATLAB to see if the drawing was read in
and all phases are assigned correctly (Figure&nbsp;[2][]).

```{figure-md} fig:jelly-plot
<img src="jelly.png" style="width:45.0%" alt="Figure" />

Plot of the <code>Phase</code> variable in MATLAB.
```

Finally, we want to write output in a format that can be read in by &rsquo;s
`ascii data` compositional initial conditions plugin. We write the data into
the file `jelly.txt`:

``` matlab
```

To read in the file we just created (a copy is located in &rsquo;s data
directory), we set up a model with a box geometry with the same extents we
specified for the drawing in px and one compositional field. We choose the
`ascii data` compositional initial conditions and specify that we want to read
in our jellyfish. The relevant parts of the input file are listed below:

``` prmfile
```

If we look at the output in `ParaView`, we can see our jellyfish, with the
mesh refined at the boundaries between the different phases
(Figure&nbsp;[3][]).

```{figure-md} fig:jelly-paraview
<embed src="cookbooks/geomio/doc/jelly-paraview.pdf" style="width:55.0%" />

<figcaption aria-hidden="true"> <em>model output of the jellyfish and corresponding mesh in ParaView.
```

For a geophysical setup, the MATLAB code could be extended to write out the
phases into several different columns of the ASCII data file (corresponding to
different compositional fields). This initial conditions file could then be
used in with a material model such as the `multicomponent` model, assigning
each phase different material properties.

An animation of a model using the jellyfish as initial condition and assigning
it a higher viscosity can be found here:
<https://www.youtube.com/watch?v=YzNTubNG83Q>.

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-bauvillegeomio" class="csl-entry">

Bauville, A, and TS Baumann. 2019. &ldquo;geomIO: An Open-Source MATLAB
Toolbox to Create the Initial Configuration of 2d/3d Thermo-Mechanical
Simulations from 2d Vector Drawings.&rdquo; *Geochemistry, Geophysics,
Geosystems*. <https://doi.org/10.1029/2018GC008057>.

</div>

</div>

  [here]: https://bitbucket.org/geomio/geomio/downloads
  [1]: http://geomio-doc.bitbucket.org/
  [how to generate drawings in Inkscape]: http://geomio-doc.bitbucket.org/tuto2D.html#drawing
  [cookbooks/geomio/doc/jellyfish.svg]: cookbooks/geomio/doc/jellyfish.svg
  [1]: #fig:jelly-picture
  [run geomIO as described in the documentation]: http://geomio-doc.bitbucket.org/tuto2D.html#assigning-phase-to-markers
  [2]: #fig:jelly-plot
  [3]: #fig:jelly-paraview
