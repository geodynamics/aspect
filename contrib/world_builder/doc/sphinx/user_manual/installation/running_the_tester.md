(part:user_manual:chap:installation:sec:running_the_tester)=
Running the tester
==================

Once you have installed the GWB, it is always a good idea to run the tester to see if the installation works. There are a few ways to call the tester: `make test`, `ninja test` and `ctest`. They can all run the same tests suite, but `ctest` is the most universal and versatile of the options. Running CTest may give you output like below:

```
Test project /path/toworld-builder/WorldBuilder/build-rl
      Start  1: unit_tests
 1/74 Test  #1: unit_tests ........................................   Passed    2.99 sec
      Start  2: testing_help
 2/74 Test  #2: testing_help ......................................   Passed    0.01 sec
      Start  3: testing_no_file
 3/74 Test  #3: testing_no_file ...................................   Passed    0.01 sec
      Start  4: testing_one_file
 4/74 Test  #4: testing_one_file ..................................   Passed    0.01 sec
      Start  5: testing_too_many_arguments
 5/74 Test  #5: testing_too_many_arguments ........................   Passed    0.01 sec
      Start  6: app_continental_plate_2d
 6/74 Test  #6: app_continental_plate_2d ..........................   Passed    0.22 sec
      Start  7: app_continental_plate_3d
 7/74 Test  #7: app_continental_plate_3d ..........................   Passed    0.21 sec
 ...
 ```

 If everything went fine, it will report success. If not, it will show something like this at the end:

```
70/74 Test #70: compile_simple_fortran_example ....................   Passed    0.28 sec
      Start 71: run_simple_fortran_example
71/74 Test #71: run_simple_fortran_example ........................   Passed    0.06 sec
      Start 72: run_python_version_test
72/74 Test #72: run_python_version_test ...........................   Passed    0.00 sec
      Start 73: run_simple_python_test
73/74 Test #73: run_simple_python_test ............................   Passed    0.02 sec
      Start 74: run_simple_python_example
74/74 Test #74: run_simple_python_example .........................***Failed    0.02 sec

99% tests passed, 1 tests failed out of 74

Total Test time (real) =  11.74 sec

The following tests FAILED:
         74 - run_simple_python_example (Failed)
Errors while running CTest
Output from these tests are in: /path/to/world-builder/WorldBuilder/build-rl/Testing/Temporary/LastTest.log
Use "--rerun-failed --output-on-failure" to re-run the failed cases verbosely.
```

Tests can fail for a number of reasons without there being a major issue but you will need to look at them to know whether there is an issue or not. One example is that on different computers, there may be very small differences in the outcome. That is, if the differences between the reference output and the generated output are **very** small, the installation is usually fine. Below we will explain how to find out how big the differences are. One way to mitigate this problem is to install a program called `numdiff`. `Numdiff` is a program which can recognize numbers and only reports a difference if it is significant. The GWB will prefer to use `numdiff` over the normal `diff`, if `numdiff` is installed. 

An other example is the problem shown above, which we will now investigate. The easiest way is to use `ctest --rerun-failed --output-on-failure`. This will only rerun the tests which failed and provide output when it fails.

```
Test project /path/toworld-builder/WorldBuilder/build-rl
    Start 74: run_simple_python_example
1/1 Test #74: run_simple_python_example ........***Failed    0.02 sec
CMake Error at /path/toworld-builder/WorldBuilder/tests/python/run_python_tests.cmake:32 (message):
  Failed: Test program /usr/bin/python3.10 exited != 0.

   Test args where: /path/toworld-builder/WorldBuilder/tests/python/example.py

  Traceback (most recent call last):

    File "/path/toworld-builder/WorldBuilder/tests/python/example.py", line 1, in <module>
      from gwb import WorldBuilderWrapper

  ModuleNotFoundError: No module named 'gwb'

   
   The test outpup was: 
   




0% tests passed, 1 tests failed out of 1

Total Test time (real) =   0.02 sec

The following tests FAILED:
         74 - run_simple_python_example (Failed)
Errors while running CTest
```

The reason for this particular error is that although you compiled the GWB Python module, you didn't install it to a place Python can find it. If you are not planning to use the GWB Python module, you can safely ignore this error.
