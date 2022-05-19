# Debug or optimized mode

### Debug or optimized mode

ASPECT utilizes a <span
class="smallcaps">deal.II</span> feature called *debug mode*. By default,
ASPECT uses debug mode, i.e., it calls a
version of the DEAL.II library that contain
lots of checks for the correctness of function arguments, the consistency of
the internal state of data structure, etc. If you program with <span
class="smallcaps">deal.II</span>, for example to extend 
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
assertions (as of late 2011) in DEAL.II where
we check for errors like the above and, if the condition is violated, abort
the program with a detailed message that shows the failed check, the location
in the source code, and a stacktrace how the program got there. The downside
of debug mode is, of course, that it makes the program much slower -
depending on application by a factor of 4-10. An example of the speedup
one can get is shown in {ref}`5.2.1][].

ASPECT by default uses debug mode because most
users will want to play with the source code, and because it is also a way to
verify that the compilation process worked correctly. If you have verified
that the program runs correctly with your input parameters, for example by
letting it run for the first 10 time steps, then you can switch to optimized
mode by compiling ASPECT with the command[14]

     make release

and then compile using

     make

To switch back to debug mode type:

     make debug

<div class="center">

</div>
