# Visualizing results

Among the postprocessors that can be selected in the input parameter file (see
Sections&nbsp;[4.2][] and {ref}`parameters:Postprocess/Visualization`52]) are
some that can produce files in a format that can later be used to generate a
graphical visualization of the solution variables $\mathbf u, p$ and $T$ at
select time steps, or of quantities derived from these variables (for the
latter, see {ref}`sec:viz-postpostprocessors`53]).

By default, the files that are generated are in VTU format, i.e., the
XML-based, compressed format defined by the VTK library, see
<http://public.kitware.com/VTK/>. This file format has become a broadly
accepted pseudo-standard that many visualization program support, including
two of the visualization programs used most widely in computational science:
Visit (see <https://visit.llnl.gov/>) and ParaView (see
<http://www.paraview.org/>). The VTU format has a number of advantages beyond
being widely distributed:

-   It allows for compression, keeping files relatively small even for sizable
    computations.

-   It is a structured XML format, allowing other programs to read it without
    too much trouble.

-   It has a degree of support for parallel computations where every processor
    would only write that part of the data to a file that this processor in
    fact owns, avoiding the need to communicate all data to a single processor
    that then generates a single file. This requires a master file for each
    time step that then contains a reference to the individual files that
    together make up the output of a single time step. Unfortunately, there
    doesn't appear to be a standard for these master records; however,
    both ParaView and Visit have defined a format that each of these programs
    understand and that requires placing a file with ending `.pvtu` or
    `.visit` into the same directory as the output files from each processor.
    {ref}`4.2][] gives an example of what can be found in the output
    directory.

<div class="center">

</div>

:::{toctree}
visit.md
stat-data.md
large-data-issues.md
:::
