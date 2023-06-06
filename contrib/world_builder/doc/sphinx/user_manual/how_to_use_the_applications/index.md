(part:user_manual:chap:how_to_use_the_apps:sec:index)=
How to use the applications
===========================

The Geodynamic World Builder is in its core a code library, which can be used by other codes. That means that to use it you need to write a program which uses it. The apps the gwb provide do exactly this. They are small programs which allow users to directly use the GWB without writing a program themselves. There are currently two apps which come with the GWB: `gwb-dat` and `gwb-grid`. The first application allows users to provide a world builder file and a `.dat`, which is a space separated value file with some extra options, to query individual points of in the world builder world. The second application allows users to use a `.grid` file to create a grid output (`.vtu`) of a whole domain, which then can be visualized in a program like Paraview.

```{toctree}
:caption: User manual
:hidden:

gwb-dat_app
gwb-grid_app
```
