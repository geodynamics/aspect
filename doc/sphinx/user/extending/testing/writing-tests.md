(sec:extending:writing-tests)=
# Writing tests

To write a test for a new feature, copy one of the existing parameter files in
the `tests/` folder in the ASPECT source
directory, or simply any other parameter file, modify it to use the new
feature, check that the new feature does what it is supposed to do, and then
just add the parameter file to the tests directory. You will then need to add
another folder to that directory that is named exactly like the parameter
file, and add the model output files that prove that the feature is working
(usually, these are the log file and the statistics file, and you will have to
rename `log.txt` to `screen-output` for historical reasons). The test and
output files should be as small and quick to run as possible. If you need to
include graphical output to test your feature, you will have to use the
gnuplot output format, so that the tester can compare the actual numbers (in
the vtu format, the output files are compressed, and can not be compared using
Numdiff). An easy way to create all of the files you need is to copy the
folder of an existing test and rename it to the name of your parameter file.

To actually run the test, you have to go to your
ASPECT build directory and run

``` ksh
make setup_tests
```

so that your new test is added to the test suite. Then you can run it by
executing

``` ksh
ctest -R name_of_your_test -V
```

and you will get an output telling you if the test has passed or failed (and
why it failed). If you have just copied the output files of a different test
in the `tests/` directory to make your test, you of course expect your test to
fail. In this case, the output you see should contain a line that starts with
`******* Check` and then just shows two paths. Those two paths are the one
where the output files of the test are located (the ones you just created by
running the test) and where the reference output of the test is located (the
one you created by copying an existing test). So you can copy this whole line
and replace `******* Check` by `cp` to copy the output you just created over
the reference output. Of course, you should only do that after you have made
sure that these output files show that the feature you want to test is working
as expected.

When you make a new test part of a pull request on GitHub, then as explained
above that will lead to a run of all tests &ndash; including your new one
&ndash; on a "reference machine." The reference machine that runs
the tests may of course produce slightly different results than the machine on
which a pull request was developed and from which the output was taken. If
this has been confirmed to be the source of a failed test run, a file that
contains the differences between the test output you submitted as part of your
pull request and the "reference" tester output will be available
from GitHub (You will have to click on the link labeled "Details"
next to the line that tells you if tests have failed; for the jenkins tester
that will bring you to a new page, where you have to go to the
"Artifacts" tab in the top right corner. Depending on the tester,
the file might be called `changes-gcc.diff` or `changes-test-results.diff`).
To use this file to update your test output, you will have to download it and
put it into your top-level ASPECT directory.
There you can apply the .diff file using:

``` ksh
git apply changes-test-results.diff
```

This will update your test output so that it matches the results from the
official tester.

On the other hand, if a change leads to even a single *existing* test failing
on that system, then we know that some more investigation as to the causes is
necessary.
