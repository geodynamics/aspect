(part:user_manual:chap:introduction:sec:citing_the_world_builder)=
Referencing the World Builder
=============================

## Citation
The developers of the Geodynamic World Builder request that you cite the following publication that describes the code and cite the specific version of the software used.

:::::{tab-set}
:class: full-width

::::{tab-item} v0.5.0 Bibtex
:name: v0.5.0_bibtex

:::{include} ../../../../CITATION
::::

:::::

## Data Availability
We strongly recommend making your data available for reproducibilty and replicability. Consider depositing your data e.g., code and .wb, in an approved repository, which will assign an identifier (e.g., a DOI) and enable citation of your data. For visualization of your model, we strongly recommend also including your grid (.grid), paraview state file (.pvsm), and image (e.g., .png) files. See [geodynamics.org software publishing guidance](https://geodynamics.org/software/software-bp/software-publishing).
Then add the following to your data availability statement:

:::{code-block} latex

The code modifications, parameter, and visualization files used for the models in this study 
are available at DOI (Authors X, Y, Z) under CC BY-NC-SA 4.0.

GWB version 0.5 (\cite{se-10-1785-2019,gwb-doi-v0.5.0}) used to build these models is freely 
available under the LGPL v2.1 license and is being actively developed on GitHub and can be 
accessed via https://github.com/GeodynamicWorldBuilder/WorldBuilder.
:::
