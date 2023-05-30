(sec:run-aspect:2d-vs-3d)=
# Selecting between 2d and 3d

ASPECT can solve both two- and
three-dimensional problems.[^footnote1] You select which one you want by putting a
line like the following into the parameter file (see
<https://aspect.geodynamics.org/doc/parameter_view/parameters.xml>):

```{literalinclude} ../../../manual/cookbooks/overview/doc/dim.part.prm
```

Internally, dealing with the dimension builds on a feature in deal.II, upon which
ASPECT is based, that is called
*dimension-independent programming*. In essence, what this does is that you
write your code only once in a way so that the space dimension is a variable
(or, in fact, a template parameter) and you can compile the code for either 2d
or 3d. The advantage is that codes can be tested and debugged in 2d where
simulations are relatively cheap, and the same code can then be re-compiled
and executed in 3d where simulations would otherwise be prohibitively
expensive for finding bugs; it is also a useful feature when scoping out
whether certain parameter settings will have the desired effect by testing
them in 2d first, before running them in 3d. This feature is discussed in
detail in the [deal.II tutorial program
step-4](https://www.dealii.org/developer/doxygen/deal.II/step_4.html). Like there, all the functions and classes in
ASPECT are compiled for both 2d and 3d. Which
dimension is actually called internally depends on what you have set in the
input file, but in either case, the machine code generated for 2d and 3d
results from the same source code and should, thus, contain the same set of
features and bugs. Running in 2d and 3d should therefore yield comparable
results. Be prepared to wait much longer for computations to finish in the
latter case, however.

[^footnote1]: For a description of what exactly we mean when we consider two-dimensional models, see {ref}`sec:methods:2d-models`
