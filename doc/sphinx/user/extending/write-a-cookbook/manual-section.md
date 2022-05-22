# Adding a section to the manual

Next, you have to decide if the cookbook you want to contribute is a *Simple
setup* (that explains how to use one specific feature, but does not try to
reproduce any earth-like setting, see
{ref}`sec:cookbooks-simple`), a *Geophysical setup* (that
teaches how to setup a specific type of geodynamic model like a global
convection model, a subduction zone or a mid-ocean ridge, see
{ref}`sec:cookbooks-geophysical`) or a *Benchmark* (see
{ref}`sec:cookbooks-benchmarks`]). Depending on that choice,
you will then start a new `\subsubsection` in the [manual.tex][] file at the
end of the corresponding subsection (*Simple setups*, *Geophysical setups* or
*Benchmarks*). This is where your description of the model will go.

In addition to the text in the manual, you also have to create a subfolder in
the [doc/manual/cookbooks][] directory. This is the place where all figures
and input file/code snippets that accompany the description go into.

Note also one special case: If your setup is a benchmark, you will have to put
your input file into the [benchmarks][] folder rather than into the
[cookbooks][] folder, and you have to create the subfolder for your figures
and code snippets in the [doc/manual/cookbooks/benchmarks][] directory.

To give you some guidelines on how to write the section in the manual, you can
follow this general structure:

-   Start with a short description of what feature the cookbook introduces or
    what the model setup is meant to accomplish, including the relevant
    physics. Specifically, this paragraph should also address the question of
    what motivates the model. If the setup comes from a publication, make sure
    to mention that and include the reference.

-   If the model uses a new plugin, describe the new feature this plugin
    introduces and how this is implemented in the code. Ideally, this
    paragraph includes essential code snippets from the plugin file that
    complement and illustrate the description in the text. Place the code
    snippet in the corresponding subfolder you created in the
    [doc/manual/cookbooks][] directory and use the command

          \lstinline{\lstinputlisting[language=C++]{cookbooks/subfolder_name/code_snippet.cc}!

    to insert the code in the manual.tex file.

-   Explain what the important input parameters in this setup are, what values
    you set them to and why. This paragraph should give an overview of your
    model setup, including the initial conditions, boundary conditions,
    geometry, etc., and anything that is special about the setup. Ideally,
    this description includes snippets from the input file. You can place
    these snippets in the subfolder you created in the
    [doc/manual/cookbooks][] directory and include them in the `manual.tex`
    file using a command like

          \lstinputlisting[language=prmfile]{cookbooks/subfolder_name/doc/input_snippet.prm.out}

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

And that's it, you have just created your first cookbook! Make a [pull
request](https://docs.github.com/en/get-started/quickstart/github-flow) to contribute it to the main repository! You can find more
information on how to do that on [our github page](https://github.com/geodynamics/aspect/blob/main/CONTRIBUTING.md).

You will get bonus points if you also create a test (see
{ref}`sec:1.8.2`) that only runs the first time step (or a lower
resolution version) of your cookbook.

:::{admonition} TODO
:class: error

This section needs substantial updating by C. Mills to update the process on how to add a section to the manual.
A number of relative links also need to be added.
:::
