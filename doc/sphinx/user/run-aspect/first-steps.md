# First steps

Before trying to set up a model to answer your particular research questions,
we advise you to get familiar with ASPECT and
its functionalities by following these steps:

1.  Watch the CIG ASPECT tutorials
    (<https://www.youtube.com/playlist?list=PLdy04DoEepEyeS_HZwa0Ws0kW5Rs2wsQ6>)
    that will show you how to run ASPECT and
    construct new setups yourself.

2.  Go through the cookbooks in this manual, see {ref}`cha:cookbooks`.

3.  Go through the benchmarks in this manual, see
    {ref}`cha:benchmarks`.

4.  If you want to use some existing functionality that is not discussed in
    these resources, search in the extensive tests directory. For example, to
    search for an initial temperature condition called "spherical
    gaussian perturbation" while in the
    ASPECT directory, type:

          grep 'spherical gaussian perturbation' tests/*.prm

    This command will show you all the test input files that use this initial
    temperature condition. You can also look up any of the parameters used in
    the input files in this manual.

5.  Have a look at the ASPECT GitHub
    repository. Here you can see the planned developments
    (<https://github.com/geodynamics/aspect/projects/2>), current issues that
    others have reported (<https://github.com/geodynamics/aspect/issues>), and
    what is currently being worked on
    (<https://github.com/geodynamics/aspect/pulls>).

6.  Have a look at our discussion forum when your model behaves unexpectedly
    or you need functionality that does not exist yet. The
    ASPECT community can tell you whether they
    experienced something similar or are already working on the topic.

7.  If you experience unexpected behavior that you expect is a bug and this
    problem has not been reported as an issue on GitHub, please create a new
    issue so that everybody is aware of the potential problem and can think of
    a fix. When creating a new issue, it is very useful if you can provide a
    minimum working example, i.e. a small test setup that demonstrates the
    issue and does not require modifications to the code. You can for example
    modify one of the existing test input files, which typically take less
    than a minute to run using only a few cores. The test input file and an
    image illustrating the problem can be attached to the issue.
