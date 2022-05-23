(sec:run-aspect:visualize:stat-data)=
# Visualizing statistical data

In addition to the graphical output discussed above,
ASPECT produces a statistics file that collects
information produced during each time step. For the remainder of this section,
let us assume that we have run ASPECT with the
input file discussed in {ref}`sec:cookbooks:convection-box`, simulating convection in a
box. After running ASPECT, you will find a file
called `statistics` in the output directory that, at the time of writing this,
looked like this:

``` ksh
# 1: Time step number
# 2: Time (seconds)
# 3: Number of mesh cells
# 4: Number of Stokes degrees of freedom
# 5: Number of temperature degrees of freedom
# 6: Iterations for temperature solver
# 7: Iterations for Stokes solver
# 8: Velocity iterations in Stokes preconditioner
# 9: Schur complement iterations in Stokes preconditioner
# 10: Time step size (seconds)
# 11: RMS velocity (m/s)
# 12: Max. velocity (m/s)
# 13: Minimal temperature (K)
# 14: Average temperature (K)
# 15: Maximal temperature (K)
# 16: Average nondimensional temperature (K)
# 17: Outward heat flux through boundary with indicator 0 ("left") (W)
# 18: Outward heat flux through boundary with indicator 1 ("right") (W)
# 19: Outward heat flux through boundary with indicator 2 ("bottom") (W)
# 20: Outward heat flux through boundary with indicator 3 ("top") (W)
# 21: Visualization file name
 0 0.0000e+00 256 2467 1089  0 29 30 29 1.2268e-02 1.79026783e+00 2.54322608e+00
 1 1.2268e-02 256 2467 1089 32 29 30 30 3.7388e-03 5.89844152e+00 8.35160076e+00
 2 1.6007e-02 256 2467 1089 20 28 29 29 2.0239e-03 1.09071922e+01 1.54298908e+01
 3 1.8031e-02 256 2467 1089 15 27 28 28 1.3644e-03 1.61759153e+01 2.28931189e+01
 4 1.9395e-02 256 2467 1089 13 26 27 27 1.0284e-03 2.14465789e+01 3.03731397e+01
 5 2.0424e-02 256 2467 1089 11 25 26 26 8.2812e-04 2.66110761e+01 3.77180480e+01
```

In other words, it first lists what the individual columns mean with a hash
mark at the beginning of the line and then has one line for each time step in
which the individual columns list what has been explained above.[^footnote1]

This file is easy to visualize. For example, one can import it as a whitespace
separated file into a spreadsheet such as Microsoft Excel or
OpenOffice/LibreOffice Calc and then generate graphs of one column against
another. Or, maybe simpler, there is a multitude of simple graphing programs
that do not need the overhead of a full fledged spreadsheet engine and simply
plot graphs. One that is particularly simple to use and available on every
major platform is `Gnuplot`. It is extensively documented at
<http://www.gnuplot.info/>.

`Gnuplot` is a command line program in which you enter commands that plot data
or modify the way data is plotted. When you call it, you will first get a
screen that looks like this:

```
/home/user/aspect/output gnuplot

        G N U P L O T
        Version 4.6 patchlevel 0    last modified 2012-03-04
        Build System: Linux x86_64

        Copyright (C) 1986-1993, 1998, 2004, 2007-2012
        Thomas Williams, Colin Kelley and many others

        gnuplot home:     http://www.gnuplot.info
        faq, bugs, etc:   type "help FAQ"
        immediate help:   type "help"  (plot window: hit 'h')

Terminal type set to 'qt'
gnuplot>
```

At the prompt on the last line, you can then enter commands. Given the
description of the individual columns given above, let us first try to plot
the heat flux through boundary 2 (the bottom boundary of the box), i.e.,
column 19, as a function of time (column 2). This can be achieved using the
following command:

``` gnuplot
plot "statistics" using 2:19
```

The left panel of {numref}`fig:viz-gnuplot` shows what `Gnuplot` will display in its
output window. There are many things one can configure in these plots (see the
`Gnuplot` manual referenced above). For example, let us assume that we want to
add labels to the $x$- and $y$-axes, use not just points but lines and points
for the curves, restrict the time axis to the range $[0,0.2]$ and the heat
flux axis to $[-10:10]$, plot not only the flux through the bottom but also
through the top boundary (column 20) and finally add a key to the figure, then
the following commands achieve this:

``` gnuplot
set xlabel "Time"
  set ylabel "Heat flux"
  set style data linespoints
  plot [0:0.2][-10:10] "statistics" using 2:19 title "Bottom boundary", \
                       "statistics" using 2:20 title "Top boundary"
```

If a line gets too long, you can continue it by ending it in a backslash as
above. This is rarely used on the command line but useful when writing the
commands above into a script file, see below. We have done it here to get the
entire command into the width of the page.

```{figure-md} fig:viz-gnuplot
<img src="images/viz-gnuplot.*" alt="Figure" width="80%"/>

Visualizing the statistics file obtained from the example in {ref}`sec:cookbooks:convection-box` using Gnuplot: Output using simple commands.
```

For those who are lazy, `Gnuplot` allows to abbreviate things in many
different ways. For example, one can abbreviate most commands. Furthermore,
one does not need to repeat the name of an input file if it is the same as the
previous one in a plot command. Thus, instead of the commands above, the
following abbreviated form would have achieved the same effect:

``` gnuplot
se xl "Time"
  se yl "Heat flux"
  se sty da lp
  pl [:0.2][-10:10] "statistics" us 2:19 t "Bottom boundary", "" us 2:20 t "Top boundary"
```

This is of course unreadable at first but becomes useful once you become more
familiar with the commands offered by this program.

Once you have gotten the commands that create the plot you want right, you
probably want to save it into a file. `Gnuplot` can write output in many
different formats. For inclusion in publications, either `eps` or `png` are
the most common. In the latter case, the commands to achieve this are

``` gnuplot
set terminal png
  set output "heatflux.png"
  replot
```

The last command will simply generate the same plot again but this time into
the given file. The result is a graphics file similar to the one shown in
{numref}`fig:convection-box-stats`.

:::{note}
After setting output to a file, *all* following plot commands will want to write to this file.
Thus, if you want to create more plots after the one just created, you need to reset output back to
the screen. On Linux, this is done using the command `set terminal X11`. You can then continue
experimenting with plots and when you have the next plot ready, switch back to output to a file.
:::

What makes `Gnuplot` so useful is that it doesn't just allow entering
all these commands at the prompt. Rather, one can write them all into a file,
say `plot-heatflux.gnuplot`, and then, on the command line, call

``` ksh
gnuplot plot-heatflux.gnuplot
```

to generate the `heatflux.png` file. This comes in handy if one wants to
create the same plot for multiple simulations while playing with parameters of
the physical setup. It is also a very useful tool if one wants to generate the
same kind of plot again later with a different data set, for example when a
reviewer requested additional computations to be made for a paper or if one
realizes that one has forgotten or misspelled an axis label in a plot.[^footnote2]

`Gnuplot` has many many more features we have not even touched upon. For
example, it is equally happy to produce three-dimensional graphics, and it
also has statistics modules that can do things like curve fits, statistical
regression, and many more operations on the data you provide in the columns of
an input file. We will not try to cover them here but instead refer to the
manual at <http://www.gnuplot.info/>. You can also get a good amount of
information by typing `help` at the prompt, or a command like `help plot` to
get help on the `plot` command.

[^footnote1]: With input files that ask for initial adaptive refinement, the first time step may appear twice because we solve on a mesh that
is globally refined and we then start the entire computation over again on a once adaptively refined mesh (see the parameters
in {ref}`parameters:Mesh_20refinement` for how to do that).
[^footnote2]: In my own work, I usually save the ASPECT input file, the statistics output file and the Gnuplot script along with the
actual figure I want to include in a paper. This way, it is easy to either re-run an entire simulation, or just tweak the graphic
at a later time. Speaking from experience, you will not believe how often one wants to tweak a figure long after it was first
created. In such situations it is outstandingly helpful if one still has both the actual data as well as the script that generated
the graphic.
