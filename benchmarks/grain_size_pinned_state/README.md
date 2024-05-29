# Grain Size Pinned State

This folder contains information about benchmarks that test the implementation of
grain size evolution following the paper of

Mulyukova, E., & Bercovici, D. (2018). Collapse of passive margins by
lithospheric damage and plunging grain size. Earth and Planetary Science
Letters, 484, 341-352.

## Grain size growth

[tests/grain_size_growth_pinned.prm](https://github.com/geodynamics/aspect/blob/main/tests/grain_size_growth_pinned.prm) is a test that ensures the grain growth term is implemented correctly. See the parameter file for a detailed description.

## Equilibrium grain size

[tests/grain_size_strain_pinned_state.prm](https://github.com/geodynamics/aspect/blob/main/tests/grain_size_strain_pinned_state.prm) is a test that ensures the correct equilibrium grain size is reached for simple flow conditions. See the parameter file and the notebook `equilibrium_grain_size_pinned_state.ipynb` in this folder for more information.

## Dynamic grain size evolution

[tests/grain_size_strain_pinned.prm](https://github.com/geodynamics/aspect/blob/main/tests/grain_size_strain_pinned.prm) is a test that ensures the correct dynamic grain size change is computed for simple flow conditions. See the parameter file and the notebook `dynamic_grain_size_pinned_state.ipynb` in this folder for more information.

## Grain size plunge at passive margins

[grain_size_plunge.prm](grain_size_plunge.prm) is a benchmark case that reproduces the model in Fig. 4 of
{cite:t}`mulyukova:bercovici:2018`.
It models the grain size evolution under constant stress in a passive margin.
No analytical expression for the time-evolution of the grain size exists, but
plotting the evolution using the script [plot.plt](plot.plt) shows a grain size and viscosity
evolution that is very close to Fig.4 of the publication, see {numref}`fig:grainsize:plunge`.

```{figure-md} fig:grainsize:plunge
<img src="comparison_plot.*" width="100%" />

 Comparison of the grain size evolution and viscosity evolution of the ASPECT benchmark compared to
 Fig.4 of Mulyukova Bercovici 2018.
```
