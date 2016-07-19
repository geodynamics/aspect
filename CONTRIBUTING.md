# Contributing to ASPECT
ASPECT is a community project that lives by the participation of its
members — i.e., including you! It is our goal to build an inclusive
and participatory community so we are happy that you are interested in
participating! 

## Getting started with git and GitHub
If you are new to using git or the GitHub platform you may find it
helpful to review [video lecture
32.8](http://www.math.tamu.edu/~bangerth/videos.676.32.8.html), which
should be enough to help you get started with using GitHub or possibly
contributing to ASPECT itself. 

## Asking and answering questions about ASPECT
The ASPECT community maintains an active mailing list hosted by CIG
[here](https://lists.geodynamics.org/cgi-bin/mailman/listinfo/aspect-devel). The
mailing list is for questions about ASPECT at all levels.

## Bug reports
It is a great help to the community if you report any bugs in the
code that you may find. We keep track of all open issues with ASPECT
[here](https://github.com/geodynamics/aspect/issues). If you
can, please try to include a minimal failing example that can help us
reproduce the problem.

## Making changes to ASPECT
To make a change to ASPECT you should create a *fork* (through GitHub)
of the code base and a separate *branch* (sometimes called a feature
branch). You can propose that your branch be combined with the rest of
code by opening a *pull request*. This will give a chance for others
to review your code. While this seems very formal, keeping all of the
code review in one place makes it easier to coordinate changes to the
code (and there are usually several people making changes to the code
at once). This process is described at length in the deal.II video [lecture 32.8](http://www.math.tamu.edu/~bangerth/videos.676.32.8.html). Please do not hesitate to ask questions about the workflow on the mailing list if you are not sure what to do.

Since ASPECT is a fairly large project with lots of contributors we
use a set of [coding
conventions](https://www.dealii.org/developer/doxygen/deal.II/CodingConventions.html)
equal to those used by <a href="http://www.dealii.org">deal.II</a>
upon which ASPECT is based, so as to keep the style of the source code
consistent. This convention essentially consists of using `astyle` for
indentation, camel case for classes, and lower case letters with
underscores for everything else. If you are new to the project then we
will work with you to ensure your contributions are formatted with
this style, so please do not think of it as a road block if you would
like to contribute some code.

ASPECT is licensed under the GNU General Public License; while you
will retain copyright on your contributions, all changes to the code
must be provided under this common license.
