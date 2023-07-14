# How to write a cookbook

ASPECT has a number of cookbooks (see {ref}`cha:cookbooks` and
{ref}`cha:benchmarks`) that introduce certain features of the
code to new users or explain how to set up a certain type of application
model. If you have a model setup that fits into one of those categories and
are willing to share it and write some explanation about it, we are always
happy about that! We also keep a list of cookbooks we think would be great
additions to ASPECT as an [issue on github](https://github.com/geodynamics/aspect/issues/2110).

All cookbooks consist of an input file for the model run, which is located in
the [cookbooks](https://github.com/geodynamics/aspect/tree/main/cookbooks) folder,
a page in the online documentation describing the setup, and
&ndash; if additional plugins are required to run the model &ndash; the
corresponding .cc file(s) located in a subdirectory of the cookbooks folder
corresponding to the individual cookbook.

:::{toctree}
---
maxdepth: 1
---
parameter-file.md
plugins.md
manual-section.md
quickref.md
:::
