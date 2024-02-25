(part:user_manual:chap:installation:sec:stand-alone-install)=
Stand-alone installation with all apps
======================================

1. Make a directory to install to (e.g., `mkdir world_builder`).
2. Enter that directory (e.g., `cd world_builder`).
3. Clone the git repository from GitHub (e.g., `git clone git@github.com:GeodynamicWorldBuilder/WorldBuilder.git`). It is strongly recommended to make sure you have a working GitHub account first, with correctly setup ssh keys.
4. Enter the new World Builder directory (e.g., `cd WorldBuilder`).
5. Make a build directory and enter it.
(For steps 6-10, select a tab):

::::{tab-set}
:::{tab-item} For Windows with Visual Studio
6. Run CMake by entering: `cmake MAKE_FILE_GENERATOR="Visual Studio 15 2017 Win64"..`, or the version of Visual Studio you have installed, and make sure CMake finds all the dependencies.
7. For production runs, set build type to release by entering `-DCMAKE_BUILD_TYPE=Release`.
8. Run make with the amount of threads you want to use (e.g., use 8 processes: `make -j 8`).
9. If you want the Geodynamic World Builder to be installed on your system, run `cmake -build . -target install -j 8`
10. Run the tests to make sure everything is installed correctly (`cmake -build . -target run_tests -j 8`).
:::

:::{tab-item} For all other configurations
6. Run CMake by entering: `cmake ..` and make sure CMake finds all the dependencies.
7. For production runs, set build type to release by entering `make release`.
8. Run make with the amount of threads you want to use (e.g., use 8 processes: `make -j 8`).
9. If you want the Geodynamic World Builder to be installed on your system, run `sudo make install -j 4`
10. Run the tests to make sure everything is installed correctly (`ctest`).
:::
::::


Now the library, tester, and the two programs (described in {ref}`part:user_manual:chap:installation:sec:included`) are ready for use.

:::{note}
By default, the Geodynamic World Builder is configured in `debug` mode, i.e., it includes several assertions to make sure that the input and computed properties (temperature or composition) are reasonable. This makes the program much slower in the `debug ` mode compared to the `release` mode (10 times or more depending on the problem type). It is therefore recommended that you first test the feasibility of the generated output using a small problem (lower resolution or lesser features) in the `debug` mode and then use the `release` mode to run the full-scale problem for faster computation.
:::
