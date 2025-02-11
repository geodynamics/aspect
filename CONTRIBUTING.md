# Contributing to ASPECT

ASPECT is a community project that lives by the participation of its
members — i.e., including you! It is our goal to build an inclusive
and participatory community so we are happy that you are interested in
participating!

## Asking and answering questions about ASPECT

For questions about ASPECT on all levels, please use the
[ASPECT forum](https://community.geodynamics.org/c/aspect).
Archived discussions from the inactive aspect-devel mailing list can be downloaded
at [aspect-devel archives](https://geodynamics.org/pipermail/aspect-devel/).

Please note that no one in the forum is paid or obliged to answer
your questions. Everyone answering forum posts is volunteering
their time to help make our community a friendly and helpful
place to inhabit. Therefore, we ask that you keep your inquiries
polite, take some effort to make them easily readable, and
include all useful information as described below. If there is
something you can test on your own, please test it first,
before asking a question on the forum.

There is no guarantee that we can or will help with your problem, but
you are welcome to contact us and we will try to jointly improve
ASPECT for your application case. Depending on the level of
effort provided a thank you note or an acknowledgment in
a paper is always appreciated. Keep in mind that for significant
and scientifically creative contributions by a community member
(but usually only then) a co-authorship on a publication
is appropriate.

## Reporting bugs and asking questions

It is a great help to the community if you report any bugs that you
may find. We keep track of all open issues related to ASPECT
[here](https://github.com/geodynamics/aspect/issues).

Please follow these instructions before opening a new bug report or ask
a question:

- Make sure you have run your model in
  [DEBUG mode](https://aspect-documentation.readthedocs.io/en/latest/user/run-aspect/debug-mode.html).
  DEBUG mode reports many errors in greater detail and may already solve your
  question.
- Search in the
  [list of open and closed issues](https://github.com/geodynamics/aspect/issues?q=is%3Aissue)
  for a duplicate of your question.
- Search in the [ASPECT forum](https://community.geodynamics.org/c/aspect) for
  a duplicate of your question.
- If you did not find an answer in the previous searches, open a new question:
  - If you suspect you have found a bug, open a new
    [issue](https://github.com/geodynamics/aspect/issues/new) and explain your
    problem as described below.
  - If you are not sure how to set up a model, ask a question in the
    [ASPECT forum](https://community.geodynamics.org/c/aspect) as described below.
  - If you are not sure what to do, you can post a question in the
    [ASPECT forum](https://community.geodynamics.org/c/aspect).
- In either case, attach the following information:
  - a parameter file with a simplified and small model that reproduces the
    issue,
  - the `log.txt` file that was created during the model run,
  - one or several screenshots of the full error message you saw on your
    screen,
  - any information that helps us understand why you think this is a bug
    (screenshots, data series, comparisons to reference results, etc.),
  - instructions for how to reproduce the problem.

Without providing the information above, we will be less likely able to help
you, it may take significantly longer until you receive a reply, and we will
just ask you for this information anyway.

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

### Getting started with git and GitHub

If you are new to using git or the GitHub platform you may find it
helpful to review [video lecture
32.8](http://www.math.colostate.edu/~bangerth/videos.676.32.8.html), which
should be enough to help you get started with using GitHub or possibly
contributing to ASPECT itself. Alternatively, GitHub provides a helpful
guide on the process of contributing to an open-source project
[here](https://opensource.guide/how-to-contribute/).

### Opening pull requests

To make a change to ASPECT you should:

- Create a
[fork](https://guides.github.com/activities/forking) (through GitHub) of
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
not hesitate to ask questions about the workflow on the forum if you are
not sure what to do.

### Coding conventions

Since ASPECT is a fairly large project with lots of contributors we
use a set of
[coding conventions](https://www.dealii.org/developer/doxygen/deal.II/CodingConventions.html)
equal to those used by [deal.II](http://www.dealii.org)
upon which ASPECT is based, so as to keep the style of the source code
consistent. This convention essentially consists of using
[astyle](http://astyle.sourceforge.net/astyle.html) v.2.04 with a
[style file](https://github.com/geodynamics/aspect/blob/main/contrib/utilities/astyle.rc)
for indentation, CamelCase for classes and
namespaces, and lower_case_names_with_underscores for everything else. If you
are new to the project then we will work with you to ensure your contributions
are formatted with this style, so please do not think of it as a road block if
you would like to contribute some code.

### Installing astyle

To correctly indent the code in ASPECT you can use `make indent` or `ninja indent`
depending on the build system you are using. The indenting script needs version 2.04 astyle.
If you installed deal.II through candi, the correct astyle may already on your system and in your path.
If the indenting script complains that it could not find (the correct version of) astyle,
you can find it [here](https://sourceforge.net/projects/astyle/files/astyle/astyle%202.04/).

An easy way to install it is through using the following command in Linux (do not do this in the aspect directory):
`mkdir astyle && cd astyle && wget 'https://sourceforge.net/projects/astyle/files/astyle/astyle 2.04/astyle_2.04_linux.tar.gz' && tar -zxvf astyle_2.04_linux.tar.gz && cd astyle/build/gcc && make && sudo make install`.
This will create a new directory called astyle, download, unpack, compile and install it.
When you add the bin directory to to your path, the indent command should find astyle.

### Changelog entries

If your new pull request creates a change that is noticeable to ASPECT users,
please add an entry to our
[changelog](https://aspect.geodynamics.org/doc/doxygen/changes_current.html). To
add such an entry, make a copy of one of the files in the
[doc/modules/changes/](https://github.com/geodynamics/aspect/tree/main/doc/modules/changes)
folder and modify it to explain
your change. The name of the file consists of the date of the change
(approximately) and the name of the author. Start the file with a tag
(added/changed/fixed/removed) that explains the nature of your change and
include your name and date in the file in exactly the format of the other
files, this ensures you will get credit for your work.

## Acknowledgment of contributions

The ASPECT community is grateful for every contribution! But, beyond
this gratitude, there are also several *formal*
ways in which your contribution will be acknowledged by the ASPECT community:

- Every commit that is merged into the ASPECT repository makes you part of
  the growing group of ASPECT
  [contributors](https://github.com/geodynamics/aspect/graphs/contributors).
- Our biweekly forum newsletter mentions all raised issues and proposed
  and merged pull requests, including an acknowledgment of the author of the
  issue/pull request.
- For every release the most significant entries of our
  [changelog](https://aspect.geodynamics.org/doc/doxygen/changes_current.html)
  are selected to generate our release announcements. Additionally, all entries
  remain available for all previous releases of ASPECT inside the
  [repository](https://github.com/geodynamics/aspect/tree/main/doc/modules)
  and on our website (under "Development"), which documents the functionality
  that you contributed to ASPECT.
- If you contributed a significant part of the manual (such as a new cookbook,
  benchmark, or subsection), you will be listed as one of the contributing
  authors of the manual.
- The Principal Developers of ASPECT come together on a regular basis and discuss
  whether others should be invited to join the
  group of *Principal Developers*. Criteria
  for this decision include:

  - A profound understanding of ASPECT's structure and vision;
  - A proven willingness to further the project's goals and help other users;
  - Significant contributions to ASPECT (not necessarily only source code,
    also forum advice, documentation, benchmarks, tutorials);
  - Regular and active contributions to ASPECT for more than one year,
    not restricted to user meetings.

  The group of current Principal Developers is listed in the
  [AUTHORS.md](https://github.com/geodynamics/aspect/blob/main/AUTHORS.md)
  file in the main repository.

## License

ASPECT is published under the
[GPL v2 or newer](https://github.com/geodynamics/aspect/blob/main/LICENSE);
while you will retain copyright on your contributions, all changes to the code
must be provided under this common license.
