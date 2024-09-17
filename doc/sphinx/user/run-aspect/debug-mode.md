(sec:run-aspect:debug-mode)=
# Debug or optimized mode

ASPECT utilizes a deal.II feature called *debug mode*. By default,
ASPECT is configured in debug mode, i.e., it calls a
version of the deal.II library that contains
lots of checks for the correctness of function arguments, the consistency of
the internal state of data structure, and more. ASPECT also contains
many similar consistency checks enabled only in debug mode. Finally,
you can only expect a good debugging experience (for example using *gdb*)
because we only enable debug symbols in *debug mode* (at least by default).

If you program deal.II, for example to extend
ASPECT, it has been our experience over the years
that, by number, most programming errors are of the kind where one forgets to
initialize a vector, one accesses data that has not been updated, one tries to
write into a vector that has ghost elements, etc. If not caught, the result of
these bugs is that parts of the program use invalid data (data written into
ghost elements is not communicated to other processors), that operations
simply make no sense (adding vectors of different length), that memory is
corrupted (writing past the end of an array) or, in rare and fortunate cases,
that the program simply crashes.

Debug mode is designed to catch most of these errors: It enables some 7,300
assertions (as of late 2011) in deal.II where
we check for errors like the above and, if the condition is violated, abort
the program with a detailed message that shows the failed check, the location
in the source code, and a stacktrace how the program got there. The downside
of debug mode is, of course, that it makes the program much slower -
depending on application by a factor of 4-10. An example of the speedup
one can get is shown in {ref}`sec:cookbooks:convection-box`.

ASPECT by default uses debug mode because most
users will want to play with the source code, and because it is also a way to
verify that the compilation process worked correctly. If you have verified
that the program runs correctly with your input parameters, for example by
letting it run for the first 10 time steps, then you can switch to optimized
mode by compiling ASPECT with the command[^footnote1]

     make release

and then compile using

     make

To switch back to debug mode type:

     make debug

:::{note}
It goes without saying that if you make significant modifications to the program, you
should do the first runs in debug mode to verify that your program still works as expected.
Or, if you encounter bugs or a strange behavior when using ASPECT, make sure you run in
debug mode first.

On the other hand, if you are running any large computations that take
significant resources (number of processors and/or time), make sure
you run in *optimized mode*.
:::

Finally, ASPECT recently learned to compile both a debug build and a
release build simultaneously. This feature can be enabled by setting
`CMAKE_BUILD_TYPE` to `DebugRelease` or by executing `make debugrelease`.
After compilation, your build directory will contain two binaries, `aspect`
and `aspect-release` corresponding to a debug build and a release build, respectively.

:::{note}
The header that ASPECT prints to the screen at the beginning of each
computation shows you which mode was selected.
:::

More details about optimized versus debug mode are also discussed in
{ref}`sec:gmg`.

[^footnote1]: The make targets described here modify the CMake variable
`CMAKE_BUILD_TYPE` to either `Debug` or `Release`. You can of course also
modify these directly when configuring your project with CMake.
