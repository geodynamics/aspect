
# Using candi to compile dependencies

#### Using candi to compile dependencies

In its default configuration `candi` downloads and compiles a <span
class="smallcaps">deal.II</span> configuration that is able to run <span
ASPECT, but it also contains a number of packages
that are not required (and that can be safely disabled if problems occur
during the installation). We require at least the packages <span
class="smallcaps">p4est</span>, TRILINOS, and
finally DEAL.II.

At the time of this writing `candi` will install <span
class="smallcaps">p4est</span> 2.2, TRILINOS
12.18.1, and DEAL.II 9.3.0. We strive to keep
the development version of ASPECT compatible
with the latest release of DEAL.II and the
current DEAL.II development version at any
time, and we usually support several older versions of <span
class="smallcaps">p4est</span> and TRILINOS.

1.  *Obtaining candi:* Download `candi` by running

            git clone https://github.com/dealii/candi

    in a directory of your choice.

2.  *Installing DEAL.II and its dependencies:*
    Execute `candi` by running

            cd candi
            ./candi.sh -p INSTALL_PATH

    (here we assume you replace `INSTALL_PATH` by the path were you want to
    install all dependencies and DEAL.II,
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
    setting up your shell environment[6] (again we assume you replace
    `INSTALL_PATH` by the patch chosen in the previous step). Then close the
    terminal and open it again to activate the change.

4.  *Testing your installation:* Test that your installation works by
    compiling the `step-32` example that you can find in
    `$DEAL_II_DIR/examples/step-32`. Prepare and compile by running
    `cmake . && make` and run with `mpirun -n 2 ./step-32`.

Congratulations, you are now set up for compiling <span
ASPECT itself.
