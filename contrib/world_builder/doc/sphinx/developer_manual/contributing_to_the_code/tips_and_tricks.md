(part:dev_manual:chap:contribute_to_code:sec:tips_tricks)=
Tips and tricks
===============

```{todo}
Give some tick and tricks, like start copying plugins which already so almost what you want.
```

# Tricks to Run Tests

## Trick 1: Different method to run all tests

This command would compile and run all the tests:

    make test

Not so much information is given for the failing tests, so the next step is needed for more information:

    ctest -V

by appending two "V"s instead of one, more information could be given:

    ctest -VV

## Trick 2: Run one test with more information

This command would run one test and print more information to the screen with the "-R" option:

    ctest -V -R name_of_test

Notice that the "-R" option means using the regular expression, so it will run tests that match the string of "name_of_test" and might 
end up running multiple tests with matching names.
This is useful to check for failing tests. At the end of the message, the differences between the generated file of the test and the standard file would be shown. One can also use the "-VV" for more information.

In case the test fails, then there are usually paths to the files that don't match the standard output files. With this information, one can look into the contents of that file and decide what to do next.

## Trick 3: Run the unit_tests separately

In trick 1, the unit_tests are grouped as one big test. Again, not so much information is given about the individual unit tests.
Once the tests are compiled in step 1, a separate executable could be found in the build folder which provides more information regarding the unit tests:

    ./build/bin/unit_tests

you can actually list all test cases:

    ./build/bin/unit_test --list-test-cases

and then run an individual test case with:

    ./build/bin/unit_test --test-case="name of test case"

or even exclude some test cases:

    ./build/bin/unit_test --test-case-exclude="name of test case"

for more info, check with:

    ./bin/unit_test --help

## Trick 4: Interact with visualization when generating new tests

The tests in the World Builder contain a .wb file for the setup and a .dat file for the mesh point. 
For example, the mass_conserving_slab_temperature test contains two files:

    tests/gwb-dat/mass_conserving_slab_temperature.wb
    tests/gwb-dat/mass_conserving_slab_temperature.dat

And a folder:

     tests/gwb-dat/mass_conserving_slab_temperature


The first few lines of the contents in the .dat file are:

```{code-block} python
---
lineno-start: 1
---
    # Profiles across a slab dipping at 30 degrees 
    # Now define parameters:
    # dim = 2
    # compositions = 0
    # x z d T 
    2490e3 0 5e3
    2495e3 0 5e3
    2460e3 0 25e3
    2475e3 0 35e3
```

Here, the header contains information on dimensions, compositions, etc, and the data are coordinates in 2D. In this case, the 'z' coordinate is actually not used, and the 'd' coordinate is used for assigning a depth of the points. And it's the 'T' field that is going to be outputted and checked for.

For more information, please refer to the manual page for the gwb-dat app:
{ref}`part:user_manual:chap:how_to_use_the_apps:sec:gwb-dat_app`

A trick is to interact with the gwb-grid application and pick up appropriate representative points from the visualization:

    ./build/bin/gwb-grid  tests/gwb-dat/some_test.wb  some_grid.grid

This command is composed of the "gwb-grid" executable in the "bin" sub-folder, a wb file for the test, and a grid for the gwb-grid to generate output. There are already a few examples of grids in the tests/gwb-grid folder one could modify from. The output of this command is a vtu file that could be easily visualized in software like Paraview.

For more information, please refer to the manual page for the gwb-grid app:
{ref}`part:user_manual:chap:how_to_use_the_apps:sec:gwb-grid_app`

Once the appropriate points are picked up, their coordinates need to be entered in the dat file with the aforementioned layout.
This would lead to the change of outputs of the tests in the "screen-output". At last, apply the trick mentioned in trick 2, the standard "tests/gwb-dat/some_test/screen-output" file would be fixed.
