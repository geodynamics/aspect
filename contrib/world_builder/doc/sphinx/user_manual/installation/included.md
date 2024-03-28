(part:user_manual:chap:installation:sec:included)=
What is in the Geodynamic World Builder package?
================================================

The Geodynamic World Builder contains more than just an initial conditions builder library.
When it is downloaded, the folder contains the main library, two programs, and a tester package. Each one is described in more detail in the next few sections.

The World Builder library
-------------------------

This is the main code of the Geodynamic World Builder.
In this library is all the code related to setting up and querying the World Builder world.
It also contains the wrapper code to interface with C, Python, and Fortran programs.

The World Builder app
---------------------

This is a program which can be used to query the World Builder from the command line, by providing it a World Builder file, and then a data file.
This data file should contain in the header information on the dimension you want to use and the number of compositions, and in the main part, the required information (for example, for a 3d case: x, y, and z; position, depth, and gravity).
 It then outputs a file with these properties, and the temperature and compositional values behind them.
 For more information on how to use the World Builder app see {ref}`part:user_manual:chap:how_to_use_the_apps:sec:gwb-dat_app`.

The World Builder Visualization program
---------------------------------------

This program helps with visualizing the World Builder file by producing vtu files which can be opened with visualization programs like Paraview.
It requires a World Builder file and a grid file.
A grid file is a file which contains information about what part of the World Builder domain should be visualized and with what resolution.
For more information on how to use the World Builder visualizer see {ref}`part:user_manual:chap:how_to_use_the_apps:sec:gwb-grid_app`.

The tester
----------

It is important for every software to be properly tested.
The World Builder is no exception.
We currently have two types of tests implemented.
The first, and currently the most important one, is the unit tester.
This tester allows testing of individual functions of the World Builder library in relative isolation.
The second type of tester is an integration tester, which works through the World Builder app.
This tester tests whether the whole library works as expected by providing a World Builder file and data points to obtain temperature and composition.
The tester package is run every time on proposed new code before that code is added to the main World Builder repository, and all tests have to pass before the code can be merged.
This happens through GitHub actions (see <https://github.com/GeodynamicWorldBuilder/WorldBuilder/actions>) and AppVeyor (see <https://ci.appveyor.com/project/MFraters/worldbuilder>).

Having tests alone is not good enough to make sure that the World Builder actually does what it is supposed to do.
The tester should theoretically cover all the possible use cases.
In practice we test the coverage by counting the number of 'relevant' lines and how many of these lines are touched when running the tester.
This is counted and reported by the program Gcov (<https://gcc.gnu.org/onlinedocs/gcc/Gcov.html>).

This approach is not perfect and has two main problems.
The first problem is that a 100% coverage is practically not achievable since the code might have fail assertions in places which should never be reached.
The second problem is that even though a line of code is touched by the tester, it may not mean that all possible cases in that line are tested.
Think for example of an inline if statement, or an assertion macro lines.
These lines count as being touched by the tester but only one case may actually be tested.

As long as these limits are kept in mind, there is no problem in using this kind of coverage to test the tester quality.
With this limitation, we try to keep the code above 95% coverage.
At the time of writing, the coverage is above 98%.
The coverage is measured and reported every time when new code is proposed.
This happens through Coveralls (see <https://coveralls.io/github/GeodynamicWorldBuilder/WorldBuilder>).
