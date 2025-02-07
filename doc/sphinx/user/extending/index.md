(cha:extending)=
# Extending and contributing

After you have familiarized yourself with
ASPECT using the examples of
{ref}`cha:cookbooks` you will invariably want to set up your
own models. During this process you might experience that not all of your
ideas are already possible with existing functionality, and you will need to
make changes to the source code.

ASPECT is designed to be an extensible code. In
particular, it uses a plugin architecture and a set of signals through which
it is relatively easy to replace or extend certain components of the program.
Examples of things that are simple to extend are the material description, the
model geometry, the gravity field, the initial conditions, the boundary
conditions, the functions that postprocess the solution, and the behavior of
the adaptive mesh refinement. This list may also have grown since this section
was written. Changing the core functionality, i.e., the basic equations
{math:numref}`eq:stokes-1`&ndash;{math:numref}`eq:temperature`, and how they are solved is
arguably more involved. We will discuss this in {ref}`sec:extending:solver`.

There are several ways to add new functionality in plugins, and we want to
highlight advantages and disadvantages of each of them:

1.  Modify existing files: The simplest way to start modifying
    ASPECT is to modify one of the existing source
    files and then recompile the program as described in
    {ref}`sec:install:local-installation:compiling`. This process does not require any
    additional setup, and is therefore ideal for learning how to make simple
    modifications. However, it comes with several severe disadvantages. If you
    modify files the history of your local copy of
    ASPECT diverges from the official development
    version. You will therefore run into conflicts if you want to update your
    version later, for example, because there are new features or bug fixes
    available in the development version. Also these modifications make your
    results less reproducible. If you used your results in a publication, you
    could no longer say *which* version of
    ASPECT was used to produce these results, because
    you modified it yourself. Therefore, we discourage this form of
    modification for productive use (it can still be helpful for teaching).

2.  Create a feature branch: If you are familiar with the version control
    system `git` that we use to organize the development of
    ASPECT (an excellent tutorial is available at:
    <http://swcarpentry.github.io/git-novice/>) you might think of creating a
    separate branch inside your ASPECT
    repository and making your changes in this branch. This way you keep the
    history of your local modifications separate from the changes made to the
    main version. You can also uniquely describe the
    ASPECT version you used for a set of models, and
    you can upload your branch to make your changes reproducible. This
    approach is also the ideal starting point if you intend to contribute your
    changes back, as it already is the first step of our guide to contributing
    back (see also {ref}`sec:extending:contributing`). However, for projects with
    functionality that is not intended to be merged into the main version
    (e.g. because it is too specific to be of general use), we have found that
    this approach is not ideal, as you will still run into conflicts when you
    want to update your ASPECT version, and you
    need to merge the main version into your branch, or rebase the branch
    every time you want to update. Thus, while ideal for contributing to
    ASPECT, we do not recommend this approach for
    keeping model-specific functionality around.

3.  Create a shared library that contains your changes: The main benefit of
    the plugin architecture described in the paragraph above is that if you
    want to extend ASPECT for your own
    purposes, you can do this in a separate set of files that describe your
    situation, rather than by modifying the
    ASPECT source files themselves. This is
    advantageous, because (i) it makes it possible for you to update
    ASPECT itself to a newer version without losing
    the functionality you added (because you did not make any changes to the
    ASPECT files themselves), (ii) because it
    makes it possible to keep unrelated changes separate in your own set of
    files, in a place where they are simple to find, and (iii) because it
    makes it much easier for you to share your modifications and additions
    with others, you can for example include them as supplementary material in
    your publications. Of course you can (and should) also use version control
    on your separate set of files to keep track of which version of files was
    used for a given set of models. Two examples for keeping a separate shared
    library for model specific changes are discussed in
    {ref}`sec:cookbooks:prescribed_velocity`, and in
    {ref}`sec:cookbooks:inner_core_convection`. We will discuss
    the concept of plugins in {ref}`sec:extending:idea-of-plugins`, and how to write a plugin
    in {ref}`sec:extending:write-a-plugin`.

Since ASPECT is written in C++ using the
deal.II library, you will have to be proficient in
C++. You will also likely have to familiarize yourself with this library for
which there is an extensive amount of documentation:

-   The manual at
    <https://www.dealii.org/developer/doxygen/deal.II/index.html> that
    describes in detail what every class, function and variable in
    deal.II does.

-   A collection of modules at
    <https://www.dealii.org/developer/doxygen/deal.II/topics.html> that give
    an overview of whole groups of classes and functions and how they work
    together to achieve their goal.

-   The deal.II tutorial at
    <https://dealii.org/developer/doxygen/deal.II/Tutorial.html> that
    provides a step-by-step introduction to the library using a sequence of
    several dozen programs that introduce gradually more complex topics. In
    particular, you will learn deal.II's
    way of *dimension independent programming* that allows you to write the
    program once, test it in 2d, and run the exact same code in 3d without
    having to debug it a second time.

-   The step-31 and step-32 tutorial programs at
    <https://www.dealii.org/developer/doxygen/deal.II/step_31.html> and
    <https://www.dealii.org/developer/doxygen/deal.II/step_32.html> from which
    ASPECT directly descends.

-   An overview of many general approaches to numerical methods, but also a
    discussion of deal.II and tools we use in
    programming, debugging and visualizing data are given in Wolfgang
    Bangerth's video lectures. These are linked from the
    deal.II website at <https://www.dealii.org/> and
    directly available at
    <http://www.math.colostate.edu/~bangerth/videos.html>.

-   The deal.II Frequently Asked Questions at
    <https://github.com/dealii/dealii/wiki/Frequently-Asked-Questions> that
    also have extensive sections on developing code with
    deal.II as well as on debugging. It also answers
    a number of questions we frequently get about the use of C++ in deal.II.

-   Several other parts of the deal.II website
    at <https://www.dealii.org/> also have information that may be relevant if
    you dive deeper into developing code. If you have questions, the mailing
    lists at <https://www.dealii.org/mail.html> are also of general help.

-   A general overview of deal.II is also
    provided in the paper {cite:t}`bangerth:etal:2007`.

As described in {ref}`sec:run-aspect:debug-mode` you should always compile
and run ASPECT in *debug mode* when you are
making changes to the source code, as it will capture the vast majority of
bugs everyone invariably introduces in the code.

When you write new functionality and run the code for the first time, you will
almost invariably first have to deal with a number of assertions that point
out problems in your code. While this may be annoying at first, remember that
these are actual bugs in your code that have to be fixed anyway and that are
much easier to find if the program aborts than if you have to go by their more
indirect results such as wrong answers. The Frequently Asked Questions at
<https://github.com/dealii/dealii/wiki/Frequently-Asked-Questions> contain a
section on how to debug deal.II programs.


:::{toctree}
---
maxdepth: 1
---
contributing.md
idea-of-plugins.md
write-a-plugin.md
write-a-cookbook/index.md
plugin-types/index.md
compatibility.md
signals.md
extending-solver.md
testing/index.md
benchmarking-run-time.md
future-plans.md
release-tasklist-link.md
:::
