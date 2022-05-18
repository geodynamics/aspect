(cha:cookbooks)=
# Cookbooks
In this section, let us present a number of &ldquo;cookbooks&rdquo; &ndash;
examples of how to use <span class="smallcaps">ASPECT</span> in typical or
less typical ways. As discussed in Sections&nbsp;[\[sec:running\]][1] and
[\[sec:parameters\]][2], <span class="smallcaps">ASPECT</span> is driven by
run-time parameter files, and so setting up a particular situation primarily
comes down to creating a parameter file that has the right entries. Thus, the
subsections below will discuss in detail what parameters to set and to what
values. Note that parameter files need not specify *all* parameters &ndash; of
which there is a bewildering number &ndash; but only those that are relevant
to the particular situation we would like to model. All parameters not listed
explicitly in the input file are simply left at their default value (the
default values are also documented in Section&nbsp;[\[sec:parameters\]][2]).

Of course, there are situations where what you want to do is not covered by
the models already implemented. Specifically, you may want to try a different
geometry, a different material or gravity model, or different boundary
conditions. In such cases, you will need to implement these extensions in the
actual source code. Section&nbsp;[\[sec:extending\]][3] provides information
on how to do that.

The remainder of this section shows a number of applications of <span
class="smallcaps">ASPECT</span>. They are grouped into three categories:
Simple setups of examples that show thermal convection (Section&nbsp;[1.2][]),
setups that try to model geophysical situations (Section&nbsp;[1.3][]) and
setups that are used to benchmark <span class="smallcaps">ASPECT</span> to
ensure correctness or to test accuracy of our solvers (Section&nbsp;[1.4][]).
Before we get there, however, we will review how one usually approaches
setting up computations in Section&nbsp;[1.1][].

<div class="center">

</div>

:::{toctree}
cookbooks-overview.md
simple-setups/index.md
geophysical-setups/index.md
:::
