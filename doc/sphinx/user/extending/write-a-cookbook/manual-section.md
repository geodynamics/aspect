(sec:extending:write-manual-section)=
# Adding a section to the manual

Then you have to decide if the cookbook you want to contribute is a *Simple
setup* (that explains how to use one specific feature, but does not try to
reproduce any earth-like setting, see
{ref}`sec:cookbooks:simple-setups`), a *Geophysical setup* (that
teaches how to setup a specific type of geodynamic model like a global
convection model, a subduction zone or a mid-ocean ridge, see
{ref}`sec:cookbooks:geophysical-setups`) or a *Benchmark* (see
{ref}`cha:benchmarks`). Depending on that choice, you will add a new entry
to the toctree of the corresponding .md file in the
[doc/sphinx/user](https://github.com/geodynamics/aspect/tree/main/doc/sphinx/user)
directory: either
[cookbooks/simple-setups.md](https://github.com/geodynamics/aspect/blob/main/doc/sphinx/user/cookbooks/simple-setups.md),
[cookbooks/geophysical-setups.md](https://github.com/geodynamics/aspect/blob/main/doc/sphinx/user/cookbooks/geophysical-setups.md),
or [benchmarks/index.md](https://github.com/geodynamics/aspect/blob/main/doc/sphinx/user/benchmarks/index.md).
Your toctree entry should follow the structure
`cookbooks/YOUR-COOKBOOK-NAME/doc/YOUR-COOKBOOK-NAME.md` or
`benchmarks/YOUR-COOKBOOK-NAME/doc/YOUR-COOKBOOK-NAME.md` depending on which top-level
directory your cookbook is going in (see the other toctree entries for examples).

Now that it has a place to live in the manual, you will need to actually create
the documentation, in the form of a .md file in the
`aspect/cookbooks/YOUR-COOKBOOK-NAME/doc/` directory (replace `/cookbooks/` with
`/benchmarks/` if applicable). The format in which cookbook documentation (and
the rest of the manual) needs to be written is
[MyST](https://myst-parser.readthedocs.io/en/latest/#), a flavor of markdown.
For your convenience, we provide a {ref}`sec:myst-quickref` containing all the
proper formatting for tables, figures, math, citations, etc. You can also refer
to the documentation of an existing cookbook, such as
[2D compressible convection with a reference profile and material properties from BurnMan](https://github.com/geodynamics/aspect/blob/main/cookbooks/burnman/doc/burnman.md)
or [Simple convection in a quarter of a 2d annulus](https://github.com/geodynamics/aspect/blob/main/cookbooks/shell_simple_2d/doc/shell_simple_2d.md),
for examples of how to format your own documentation.
Your documentation should follow this general structure:

-   Start with a short description of what feature the cookbook introduces or
    what the model setup is meant to accomplish, including the relevant
    physics. Specifically, this paragraph should also address the question of
    what motivates the model. If the setup comes from a publication, make sure
    to mention that and include the reference (see the {ref}`sec:myst-quickref`
    for proper citation formatting).

-   If the model uses a new plugin, describe the new feature this plugin
    introduces and how this is implemented in the code. Ideally, this
    paragraph includes essential code snippets from the plugin file that
    complement and illustrate the description in the text. Place the code
    snippet in the same directory as your .md file and see
    {ref}`sec:quickref:external-code-blocks` for how to reference it.

-   Explain what the important input parameters in this setup are, what values
    you set them to and why. This paragraph should give an overview of your
    model setup, including the initial conditions, boundary conditions,
    geometry, etc., and anything that is special about the setup. Ideally,
    this description includes snippets from the input file; include them the
    same way you included the code snippets above.

-   Show the model results in form of figures and/or plots, accompanied by an
    explanation of what happens in the model. This can also include a link to
    an animation of the model you made and uploaded somewhere, for example on
    YouTube. When creating figures or animations, you should think about the
    color scale that you use. Some colormaps &ndash; like the rainbow color
    palette that is still the default in some visualization tools &ndash; can
    obscure features present in the data and introduce artifacts because the
    rainbow color scale is not perceptually uniform. For more background on
    this topic, start here
    <https://matplotlib.org/2.0.2/users/colormaps.html>. To state some of their
    recommendations here, in most cases it is best to choose a perceptually
    uniform color palette. For representing information that has ordering,
    they recommend sequential color palettes that change in lightness/color
    incrementally like "viridis," "inferno,"
    "plasma" and "magma." For representing data that
    deviates around zero, they recommend diverging color palettes where two
    different colors change in lightness and meet at an unsaturated color in
    the middle such as "BrBG" and "RdBu." If you use a
    recent version of ParaView or VisIt, these color palettes are included
    with the preset color maps under the names given above, and you may want
    to choose one of these options rather than the default.

-   Finally, mention some ways the users could modify or extend the cookbook,
    such as parameters to vary to get new and interesting results, or to
    better understand the numerical methods or the physical processes
    occurring in the model. These can just be suggestions, or you can also
    extend on these ideas by adding subsections that illustrate how these
    modifications influence the model results.

And that's it, you have just created your first cookbook! Make a
[pull request](https://docs.github.com/en/get-started/quickstart/github-flow)
to contribute it to the main repository! You can find more information on how to
do that on [our github page](https://github.com/geodynamics/aspect/blob/main/CONTRIBUTING.md).
Once the pull request is submitted, you can preview what your page will look
like in rendered form by scrolling to the checks near the bottom of the open
pull request and clicking `Details` next to `docs/readthedocs.org:aspect-documentation`.
The documentation takes a few minutes to build; once done, you will be able to
browse a preview of what the manual will look like once your additions are
merged. Navigate to your newly-added page and use this opportunity to ensure
all equations, citations, figures, etc. are formatted correctly and appear as you would like.

We also run automated checks of every cookbook and benchmark as part
of the pull request review. The script that does this might require
modifications for your new cookbook if your plugin is in a different
directory, changes are necessary to run your `.prm` files, or some
files can not be run directly. The rules are `GNU Make` rules
specified at the bottom of
[cookbooks/check.mk](https://github.com/geodynamics/aspect/blob/main/cookbooks/check.mk)
and
[benchmarks/check.mk](https://github.com/geodynamics/aspect/blob/main/benchmarks/check.mk),
respectively. The default rule is
```
%/: dummy
    +@$(def); make_lib $@
    @$(def); run_all_prms $@
```
which compiles the library in the directory and runs all prm
files. Assuming your cookbook is in a folder
`free_surface_with_crust/`, you can replace this rule by adding
```
free_surface_with_crust/: dummy
    +@$(def); make_lib $@/plugin
    @$(def); run_prm $@ test.prm
```
at the end of the document, which will compile a plugin found in
`free_surface_with_crust/plugin/` and then run `test.prm` in the
folder `free_surface_with_crust/`. Feel free to ask for help in your
pull request.

:::{note}

The `make` program has a peculiar format where the indented rows in
the snippets above consist of exactly one tab character, rather than
four spaces. This does not translate to what is shown above, but
`make` will generate difficult-to-read error messages if you use
spaces instead of the required tab.
:::

Finally, you will get bonus points if you also create a test (see
{ref}`sec:extending:writing-tests`) that only runs the first time step (or a
lower resolution version) of your cookbook.
