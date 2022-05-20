(sec:extending:running-tests)=
# Running tests

In order to run the tests, it is necessary to have either Diff or Numdiff to
compare the results to the known good case. Diff is installed by default on
most Linux systems, and Numdiff is usually available as a package so this is
not a severe limitation. While it is possible to use Diff, Numdiff is
preferred due to being able to more accurately identify whether a variation in
numerical output is significant. The test suite is run using the `ctest`
program that comes with `cmake`, and should therefore be available on all
systems that have compiled ASPECT.

After running `cmake` and then compiling
ASPECT, you can run the test suite by using the
command `ctest` in your build directory. By default, this will only run a
small subset of all tests given that both setting up all tests (several
hundred) and running them takes a non-trivial amount of time. To set up the
full test suite, you can run

``` ksh
make setup_tests
```

in the build directory. To run the entire set of tests, then execute

``` ksh
ctest
```

Unless you have a very fast machine with lots of processors, running the
entire test suite will take hours, though it can be made substantially faster
if you use

``` ksh
ctest -j <N>
```

where `<N>` is the number of tests you want `ctest` to run in parallel; you
may want to choose `<N>` equal to or slightly smaller than the number of
processors you have. Alternatively, you can run only a subset of all tests by
saying

``` ksh
ctest -R <regex>
```

where `<regex>` is a regular expression and the only tests that will be run
are those whose names match the expression.

When `ctest` runs a test, it will ultimately output results of the form

``` ksh
build> ctest -R additional_outputs
Test project /home/fac/f/bangerth/p/deal.II/1/projects/build
    Start 1: additional_outputs
1/3 Test #1: additional_outputs ...............   Passed    2.03 sec
    Start 2: additional_outputs_02
2/3 Test #2: additional_outputs_02 ............   Passed    1.84 sec
    Start 3: additional_outputs_03
3/3 Test #3: additional_outputs_03 ............   Passed    1.91 sec

100% tests passed, 0 tests failed out of 3

Total Test time (real) =   5.88 sec
```

While the small default subset of tests should work on almost all platforms,
you will find that some of the tests fail on your machine when you run the
entire test suite. This is because success or failure of a test is determined
by looking at whether its output matches the one saved at the time when the
test was written to the last digit, both as far as numerical output in
floating point precision is concerned (e.g., for heat fluxes or other things
we compute via postprocessors) as well as for integers such as the number of
iterations that is printed in the screen output.[8] Unfortunately, systems
almost always differ by compiler version, processor type and version, system
libraries, etc, that can all lead to small changes in output &ndash; generally
(and hopefully!) not large enough to produce *qualitatively* different
results, but *quantitatively* large enough to change the number of iterations
necessary to reach a specific tolerance, or to change the computed heat flux
by one part in a million. This leads to `ctest` reporting that a test failed,
when in reality it produced output that is qualitatively correct.

Given that some tests are expected to fail on any given system raises the
question why it makes sense to have tests at all? The answer is that there is
*one* system on which all tests are supposed to succeed: This system is a
machine that churns through all tests every time someone proposes a change to
the ASPECT code base via the
ASPECT GitHub page.[9] Upon completion of the test
suite, both the general summary (pass/fail) and a full verbose log will
available from the GitHub page. Because the official test setup is set up in a
Docker container, it is simple to replicate the results on a local machine. To
this end, follow the instructions in
{ref}`sec:\[subsec:docker_container`23] to set up Docker, and then run
the following command in any terminal (replace `ASPECT_SOURCE_DIR` with the
path to your ASPECT directory):

``` ksh
docker run -v ASPECT_SOURCE_DIR:/home/dealii/aspect \
    --name=aspect-tester --rm -it \
    tjhei/dealii:v9.2.0-full-v9.2.0-r2-gcc5 \
    bash /home/dealii/aspect/cmake/compile_and_update_tests.sh
```

This command executes the shell script `cmake/compile_and_update_tests.sh`
*inside* the docker container that contains the official
ASPECT test system. Note that by mounting your
ASPECT folder into the container you are actually
updating the reference test results on the host system (i.e. your computer).
