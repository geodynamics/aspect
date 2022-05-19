# Visualization the graphical output using *Visit*


In the following, let us discuss the process of visualizing a 2d computation
using Visit. The steps necessary for other visualization programs will
obviously differ but are, in principle, similar.

To this end, let us consider a simulation of convection in a box-shaped, 2d
region (see the "cookbooks" section, {ref}`5][], and in
particular {ref}`5.2.1][] for the input file for this particular
model). We can run the program with 4 processors using

      mpirun -np 4 ./aspect cookbooks/convection-box/convection-box.prm

Letting the program run for a while will result in several output files as
discussed in {ref}`4.2][] above.

In order to visualize one time step, follow these steps:[15]

-   *Selecting input files:* As mentioned above, in parallel computations we
    usually generate one output file per processor in each time step for which
    visualization data is produced (see, however, {ref}`4.5.3][]). To
    tell Visit which files together make up one time step,
    ASPECT creates a
    `output/solution/solution-XXXXX.visit` file in the output directory. To
    open it, start Visit, click on the "Open" button in the
    "Sources" area of its main window (see
    {numref}`fig:visit-1:a`54]) and select the file you want.
    Alternatively, you can also select files using the "File $>$
    Open" menu item, or hit the corresponding keyboard short-cut. After
    adding an input source, the "Sources" area of the main window
    should list the selected file name. More easily, you can also just open
    `output/solution.visit` which references *all* output files for all time
    steps. If you open this, Visit will display a slider that allows you to
    select which time step you want to visualize, along with forward,
    backward, and play buttons that allow you to move between time steps.

-   *Selecting what to plot:* ASPECT outputs
    all sorts of quantities that characterize the solution, such as
    temperature, pressure, velocity, and many others on demand (see
    {ref}`parameters:Postprocess/Visualization`52]). Once an
    input file has been opened, you will want to add graphical representations
    of some of this data to the still empty canvas. To this end, click on the
    "Add" button of the "Plots" area. The resulting
    menu provides a number of different kinds of plots. The most important for
    our purpose are: (i) "Pseudocolor" allows the visualization of
    a scalar field (e.g., temperature, pressure, density) by using a color
    field. (ii) "Vector" displays a vector-valued field (e.g.,
    velocity) using arrows. (iii) "Mesh" displays the mesh. The
    "Contour," "Streamline" and "Volume"
    options are also frequently useful, in particular in 3d.

    Let us choose the "Pseudocolor" item and select the
    temperature field as the quantity to plot. Your main window should now
    look as shown in {numref}`fig:visit-1:b`55]. Then hit the
    "Draw" button to make Visit generate data for the selected
    plots. This will yield a picture such as shown in
    {numref}`fig:visit-2:a`56] in the display window of Visit.

-   *Overlaying data:* Visit can overlay multiple plots in the same view. To
    this end, add another plot to the view using again the "Add"
    button to obtain the menu of possible plots, then the "Draw"
    button to actually draw things. For example, if we add velocity vectors
    and the mesh, the main window looks as in
    {numref}`fig:visit-1:c`57] and the main view as in
    {numref}`fig:visit-2:b`58].

-   *Adjusting how data is displayed:* Without going into too much detail, if
    you double click onto the name of a plot in the "Plots"
    window, you get a dialog in which many of the properties of this plot can
    be adjusted. Further details can be changed by using
    "Operators" on a plot.

-   *Making the output prettier:* As can be seen in
    {numref}`fig:visit-2`59], Visit by default puts a lot of clutter
    around the figure - the name of the user, the name of the input
    file, color bars, axes labels and ticks, etc. This may be useful to
    explore data in the beginning but does not yield good pictures for
    presentations or publications. To reduce the amount of information
    displayed, go to the "Controls $>$ Annotations" menu item to
    get a dialog in which all of these displays can be selectively switched on
    and off.

-   *Saving figures:* To save a visualization into a file that can then be
    included into presentations and publications, go to the menu item
    "File $>$ Save window." This will create successively numbered
    files in the directory from which Visit was started each time a view is
    saved. Things like the format used for these files can be chosen using the
    "File $>$ Set save options" menu item. We have found that one
    can often get better looking pictures by selecting the
    "Screenshot" method in this dialog.

More information on all of these topics can be found in the Visit
documentation, see <https://visit.llnl.gov/>. We have also recorded video
lectures demonstrating this process interactively at
<http://www.youtube.com/watch?v=3ChnUxqtt08> for Visit, and at
<http://www.youtube.com/watch?v=w-65jufR-bc> for ParaView.
