(sec:install:local-installation:documentation)=
# Building the documentation

If you intend to modify
ASPECT for your own experiments, you may want to also
generate documentation about the source code and the user manual.
This step is optional and not necessary if you do not intend to
modify the source code or contribute to the project.
The source code documentation can be build using the command

      cd doc; make aspect.tag

which assumes that you have the `doxygen` documentation generation tool
installed. Most Linux distributions have packages for `doxygen`. The result
will be the file `doc/doxygen/index.html` that is the starting point for
exploring the documentation. You can see an online version of the documentation
at <https://aspect.geodynamics.org/doc/doxygen/index.html>.

The user manual is created using the Sphinx documentation system and can be
built using the command

      cd doc/sphinx; make html

This will create the file `doc/sphinx/_build/html/index.html` that is the
starting point for exploring the user manual. The Sphinx documentation system
requires Python and a number of packages that are documented in the file
[environment.yml](https://github.com/geodynamics/aspect/blob/main/doc/sphinx/environment.yml).
You can create a Python environment with all
required packages using the command

      conda env create -f doc/sphinx/environment.yml

which assumes that you have the Anaconda python installer available.
