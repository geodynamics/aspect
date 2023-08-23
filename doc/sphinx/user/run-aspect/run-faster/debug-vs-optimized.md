(sec:run-aspect:run-faster:debug-vs-optimized)=
# Debug vs. optimized mode

Both deal.II and
ASPECT by default have a great deal of internal
checking to make sure that the code's state is valid. For example, if
you write a new postprocessing plugin (see
{ref}`sec:extending:idea-of-plugins`) in which you need to access the solution
vector, then deal.II's `Vector` class
will make sure that you are only accessing elements of the vector that
actually exist and are available on the current machine if this is a parallel
computation. We do so because it turns out that by far the most bugs one
introduces in programs are of the kind where one tries to do something that
obviously doesn't make sense (such as accessing vector element 101 when
it only has 100 elements). These kinds of bugs are more frequent than
implementing a wrong algorithm, but they are fortunately easy to find if you
have a sufficient number of assertions in your code. The downside is that
assertions cost run time.

As mentioned above, the default is to have all of these assertions in the code
to catch those places where we may otherwise silently access invalid memory
locations. However, once you have a plugin running and verified that your
input file runs without problems, you can switch off all of these checks by
switching from debug to optimized mode. This means re-compiling
ASPECT and linking against a version of the <span
class="smallcaps">deal.II</span> library without all of these internal checks.
Because this is the first thing you will likely want to do, we have already
discussed how to do all of this in {ref}`sec:run-aspect:debug-mode`.
