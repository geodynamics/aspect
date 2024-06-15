
# Using candi to compile dependencies

In its default configuration `candi` downloads and compiles a
deal.II configuration that is able to run
ASPECT, but it also contains a number of packages
that are not required. We strive to keep
the development version of ASPECT compatible
with the latest release of deal.II and the
current deal.II development version at any
time, and we usually support several older versions of <span
class="smallcaps">p4est</span> and Trilinos.

1.  *Obtaining candi:* Download `candi` by running

            git clone https://github.com/dealii/candi

    in a directory of your choice.

2.  *Obtaining a suitable candi configuration file:* As
    mentioned above the default configuration of candi includes
    a number of packages that are not necessary for ASPECT and
    some optional ASPECT dependencies are not enabled by default.
    In addition there are some configuration options that make ASPECT faster.
    We provide a candi configuration file that is optimized for
    ASPECT at

            https://github.com/geodynamics/aspect/tree/main/contrib/install/local.cfg

    While not mandatory, we recommend to download this file and
    place it inside the `candi` directory (you should then have two
    configuration files in that directory, named `candi.cfg` and `local.cfg`).
    When both `candi.cfg` and `local.cfg` are present in the directory, running
    candi (see below) will by default use the configuration options in `local.cfg`.

    If you want to work without the `local.cfg`, be aware that you should enable
      - p4est
      - trilinos
      - hdf5 (optional)
      - netcdf (optional)
      - sundials
      - deal.II
    and that you should consider enabling `NATIVE_OPTIMIZATIONS`. The packages above
    marked `(optional)` are not required to be able to run ASPECT, but some features
    like HDF5 output will not be available if the corresponding package is not installed.
    Therefore, we recommend installing all packages as listed above (or use the
    provided `local.cfg`).

2.  *Installing deal.II and its dependencies:*
    Execute `candi` by running

            cd candi
            ./candi.sh -p INSTALL_PATH

    (here we assume you replace `INSTALL_PATH` by the path were you want to
    install all dependencies and deal.II,
    typically a directory inside `$HOME/bin` or a similar place). This step
    might take a long time, but can be parallelized by adding `-jN`, where `N`
    is the number of CPU cores available on your computer. Further
    configuration options and parameters are listed at
    <https://github.com/dealii/candi>. In case you encounter problems during
    this step, please read the error message, and consult our wiki
    (<https://github.com/geodynamics/aspect/wiki>) for common installation
    problems, before asking on the forum
    (<https://community.geodynamics.org/c/aspect>).

3.  You may now want to configure your environment to make it aware of the
    newly installed packages. This can be achieved by adding the line
    `source INSTALL_PATH/configuration/enable.sh` to the file responsible for
    setting up your shell environment[^footnote1] (again we assume you replace
    `INSTALL_PATH` by the patch chosen in the previous step). Then close the
    terminal and open it again to activate the change.

4.  *Testing your installation:* Test that your installation works by
    compiling the `step-32` example that you can find in
    `$DEAL_II_DIR/examples/step-32`. Prepare and compile by running
    `cmake . && make` and run with `mpirun -n 2 ./step-32`.

Congratulations, you are now set up for compiling
ASPECT itself.

[^footnote1]: For bash this would be the file `˜/.bashrc.`
