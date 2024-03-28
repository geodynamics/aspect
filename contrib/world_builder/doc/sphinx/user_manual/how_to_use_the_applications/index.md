(part:user_manual:chap:how_to_use_the_apps:sec:index)=
How to use the applications
===========================

The Geodynamic World Builder is in its core a code library which can be used by other codes. To use it, you need to write a program. The GWB apps provided are small programs that allow users to directly use the GWB. Two apps currently come with the GWB: `gwb-dat` and `gwb-grid`. The first application allows users to provide a World Builder file and a `.dat` file, which is a space-separated value file with some extra options, to query individual points of interest in the World Builder world. The second application allows users to use a `.grid` file to create a grid output (`.vtu`) of a whole domain which can then be visualized in a program like Paraview.

```{toctree}
:caption: User manual
:hidden:

gwb-dat_app
gwb-grid_app
```
