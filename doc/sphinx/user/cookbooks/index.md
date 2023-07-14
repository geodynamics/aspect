(cha:cookbooks)=
# Cookbooks
In this section, let us present a number of &ldquo;cookbooks&rdquo; &ndash;
examples of how to use <span class="smallcaps">ASPECT</span> in typical or
less typical ways. As discussed in {ref}`cha:run-aspect`, ASPECT is driven by
run-time parameter files, and so setting up a particular situation primarily
comes down to creating a parameter file that has the entries that correctly describe your model. Thus, the
subsections below will discuss in detail what parameters to set and to what
values. Note that parameter files need not specify *all* parameters &ndash; of
which there is a bewildering number &ndash; but only those that are relevant
to the particular situation we would like to model. All parameters not listed
explicitly in the input file are simply left at their default value (the
default values are also documented in <https://aspect.geodynamics.org/doc/parameter_view/parameters.xml>).

Of course, there are situations where what you want to do is not covered by
the models already implemented. Specifically, you may want to try a different
geometry, a different material or gravity model, or different boundary
conditions. In such cases, you will need to implement these extensions in the
actual source code. {ref}`cha:extending` provides information
on how to do that.

The remainder of this section shows a number of applications of <span
class="smallcaps">ASPECT</span>. They are grouped into three categories:
Simple setups of examples that show thermal convection ({ref}`sec:cookbooks:simple-setups`),
setups that try to model geophysical situations ({ref}`sec:cookbooks:geophysical-setups`) and
setups that are used to benchmark <span class="smallcaps">ASPECT</span> to
ensure correctness or to test accuracy of our solvers ({ref}`cha:benchmarks`).
Before we get there, however, we will review how one usually approaches
setting up computations in {ref}`sec:cookbooks:overview`.

:::{note}
The input files discussed in the following sections can generally be found in the [cookbooks/](https://github.com/geodynamics/aspect/tree/main/cookbooks) directory of your ASPECT installation.
:::

:::{toctree}
cookbooks-overview.md
simple-setups.md
geophysical-setups.md
teaching-setups.md
:::
