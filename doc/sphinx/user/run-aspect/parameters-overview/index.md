(sec:run-aspect:parameters-overview)=
# Input parameter files

What ASPECT computes is driven by two things:

-   The models implemented in ASPECT. This
    includes the geometries, the material laws, or the initial conditions
    currently supported. Which of these models are currently implemented is
    discussed below; {ref}`sec:extending` discusses in great
    detail the process of implementing additional models.

-   The run-time parameters of the selected model. For example, you could select a model that prescribes
    constant coefficients throughout the domain from all the material models
    currently implemented; you could then select appropriate values for all of
    these constants. Both of these selections happen from a parameter file
    that is read at run time and whose name is specified on the command line.
    (See also {ref}`4.2`.)

In this section, let us give an overview of what can be selected in the
parameter file. Specific parameters, their default values, and allowed values
for these parameters are documented in {ref}`sec:parameters`.
An index with page numbers for all run-time parameters can be found on
page&nbsp;.



:::{toctree}
structure.md
categories.md
muparser-format.md
compatibility.md
:::
