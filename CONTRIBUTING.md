# Contributing to ASPECT
ASPECT is a community project that lives by the participation of its
members — i.e., including you! It is our goal to build an inclusive
and participatory community so we are happy that you are interested in
participating! 

## Getting started with git and GitHub
If you are new to using git or the GitHub platform you may find it
helpful to review [video lecture
32.8](http://www.math.colostate.edu/~bangerth/videos.676.32.8.html), which
should be enough to help you get started with using GitHub or possibly
contributing to ASPECT itself. Alternatively, GitHub provides a helpful
guide on the process of contributing to an open-source project
[here](https://opensource.guide/how-to-contribute/).

## Asking and answering questions about ASPECT
For questions about ASPECT on all levels, please use the 
[ASPECT forum](https://community.geodynamics.org/c/aspect). 
Archived discussions from the inactive aspect-devel mailing list can be downloaded 
at [aspect-devel archives](http://lists.geodynamics.org/pipermail/aspect-devel).

## Bug reports
It is a great help to the community if you report any bugs that you
may find. We keep track of all open issues related to ASPECT
[here](https://github.com/geodynamics/aspect/issues). 

Please follow these simple instructions before opening a new bug report:

- Do a quick search in the list of open and closed issues for a duplicate of
  your issue.
- Do a google search in the archived mailing list discussions for a
  duplicate of your issue by searching for 

  `your search term site:http://lists.geodynamics.org/pipermail/aspect-devel/`
- If you did not find an answer, open a new
  [issue](https://github.com/geodynamics/aspect/issues/new) and explain your
  problem in as much detail as possible.
- Attach as much as possible of the following information to your issue:
  - a minimal parameter file that reproduces the issue,
  - the `log.txt` file that was created during the model run,
  - the error message you saw on your screen,
  - any information that helps us understand why you think this is a bug, and
    how to reproduce it.

## Making ASPECT better
ASPECT is a community project, and we are encouraging all kinds of
contributions. Obvious candidates are implementations of new plugins as discussed in
the manual, since they are typically self-contained and do not
require much knowledge of the details of the remaining code. Other much
appreciated contributions are new examples (cookbooks, tests, or benchmarks),
extended documentation (every paragraph helps), and in particular fixing typos
or updating outdated documentation. Obviously, we also encourage
contributions to the core functionality in any form! If you consider making a
larger contribution to the core functionality, please open a new
[issue](https://github.com/geodynamics/aspect/issues/new) first, to discuss
your idea with one of the maintainers. This allows us to give you early
feedback and prevents you from spending much time on a project that might already be
planned, or that conflicts with other plans.

To make a change to ASPECT you should:
- Create a
[fork](https://guides.github.com/activities/forking/#fork) (through GitHub) of
the code base.
- Create a separate
[branch](https://guides.github.com/introduction/flow/) (sometimes called a
feature branch) on which you do your modifications.
- You can propose that your branch be merged into the ASPECT
code by opening a [pull request](https://guides.github.com/introduction/flow/).
This will give others a chance to review your code. 

We follow the philosophy that no pull request (independent of the author) is
merged without a review from one other member of the community, and approval of
one of the maintainers. This applies to maintainers as well as to first-time
contributors. We know that a review can be a daunting process, but pledge to
keep all comments friendly and supportive, as you can see in the list of [open
pull requests](https://github.com/geodynamics/aspect/pulls)! We are as
interested in making ASPECT better as you are!

While this seems very
formal, keeping all of the code review in one place makes it easier to
coordinate changes to the code (and there are usually several people making
changes to the code at once). This process is described at length in the
deal.II video [lecture
32.8](http://www.math.colostate.edu/~bangerth/videos.676.32.8.html).  Please do
not hesitate to ask questions about the workflow on the mailing list if you are
not sure what to do.

Since ASPECT is a fairly large project with lots of contributors we
use a set of [coding
conventions](https://www.dealii.org/developer/doxygen/deal.II/CodingConventions.html)
equal to those used by <a href="http://www.dealii.org">deal.II</a>
upon which ASPECT is based, so as to keep the style of the source code
consistent. This convention essentially consists of using
[astyle](http://astyle.sourceforge.net/astyle.html) v.2.04 with a
[style file](doc/astyle.rc) for indentation, CamelCase for classes and
namespaces, and lower_case_names_with_underscores for everything else. If you
are new to the project then we will work with you to ensure your contributions
are formatted with this style, so please do not think of it as a road block if
you would like to contribute some code.


## Acknowledgment of contributions

The ASPECT community is grateful for every contribution! But, beyond
this gratitude, there are also several *formal*
ways in which your contribution will be acknowledged by the ASPECT community:
- Every commit that is merged into the ASPECT repository makes you part of
  the growing group of ASPECT
  [contributors](https://github.com/geodynamics/aspect/graphs/contributors).
- Our biweekly mailing list newsletter mentions all raised issues and proposed
  and merged pull requests, including an acknowledgment of the author of the
  issue/pull request.
- If you contributed significant functionality to the source code we will ask
  you to provide an entry into our
  [changelog](http://aspect.geodynamics.org/doc/doxygen/changes_current.html)
  that shows your name. To add such an entry, place
  a file into the [doc/modules/changes/](doc/modules/changes/) folder.
  This changelog is used to generate our release notes, and remains available
  for all previous releases of ASPECT.
- If you contributed a significant part of the manual (such as a new cookbook,
  benchmark, or subsection), you will be listed as one of the contributing
  authors of the manual.
- The Principal Developers of ASPECT come together on a regular basis and discuss
  whether others should be invited to join the
  group of Principal Developers. Criteria that *Principal Developers*
  for this decision include:

  - A profound understanding of ASPECT's structure and vision;
  - A proven willingness to further the project's goals and help other users;
  - Significant contributions to ASPECT (not necessarily only source code,
    also mailing list advice, documentation, benchmarks, tutorials);
  - Regular and active contributions to ASPECT for more than one year,
    not restricted to user meetings.

  The group of current Principal Developers is listed in the [AUTHORS](AUTHORS.md)
  file in the main repository.

## License
ASPECT is published under the [GPL v2 or newer](LICENSE); while you
will retain copyright on your contributions, all changes to the code
must be provided under this common license.
