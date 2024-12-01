(sec:run-aspect:parameters-overview)=
# Input parameter files

What ASPECT computes is driven by two things:

-   The models implemented in ASPECT. This
    includes the geometries, the material laws, or the initial conditions
    currently supported. Which of these models are currently implemented is
    discussed below; {ref}`cha:extending` discusses in great
    detail the process of implementing additional models.

-   The run-time parameters of the selected model. For example, you could select a model that prescribes
    constant coefficients throughout the domain from all the material models
    currently implemented; you could then select appropriate values for all of
    these constants. Both of these selections happen from a parameter file
    that is read at run time and whose name is specified on the command line.
    (See also {ref}`sec:run-aspect:overview`.)

In this section, let us give an overview of the structure of these input
parameter files. Specific parameters, their default values, and allowed values
for these parameters are documented on our
[parameter website](https://aspect.geodynamics.org/doc/parameter_view/parameters.xml),
and the {ref}`parameters` chapter of this documentation.

When writing or modifying parameter
files we strongly recommend to use [Visual Studio Code](https://code.visualstudio.com/)
with the [ASPECT extension](https://marketplace.visualstudio.com/items?itemName=zhikui.vscode-aspect).
This extension enables syntax highlighting,
auto-completion, and includes the documentation for all ASPECT parameters, thus greatly simplifying
the process of writing ASPECT parameter files. Take note that by default the extension will make
use of the input parameters of the latest ASPECT release. If you are using a different ASPECT
version you can change the settings of the extension.

:::{toctree}
structure.md
categories.md
muparser-format.md
compatibility.md
:::
