(sec:cookbooks:CPO_induced_anisotropic_viscosity)=
# CPO induced anisotropic viscosity

*This section was contributed by Yijun Wang and Ágnes Király.*

This cookbook explains how to use the CPO-induced anisotropic viscosity material model to set up a shear box model.

## Motivation


## Model setup

We prescribe simple shear in a 3d Cartesian box/cube
with dimensions of $1 \times 1 \times 1 $ $[m^3]$. The shear strain rate is set to
$\dot{\epsilon}_{xz} = -5\times 10^{-6} [s^{-1}]$. The particle with olivine/enstatite volume fraction of 0.7/0.3 is placed
right at the center of the cubic box so it stays stationary. The DRex implementation
keep tracks of rotations of crystal grains within the particle under macroscopic deformation.

```{figure-md} fig:model_setup_3D_box
<img src="model_setup_3D_box.png" style="width:80.0%" />

Description
```

The model computes how crystal grains rotate and align under simple shear and the pole figures visualize this rotation.


## Model results


## Extending the model


## References:
- Fraters, M. R. T., & Billen, M. I. (2021).
On the implementation and usability of crystal preferred orientation evolution in geodynamic modeling. Geochemistry, Geophysics, Geosystems, 22(10), e2021GC009846.

- Kaminski, E., Ribe, N. M., & Browaeys, J. T. (2004).
D-Rex, a program for calculation of seismic anisotropy due to crystal lattice preferred orientation in the convective upper mantle. Geophysical Journal International, 158(2), 744–752.

- Karato, S., Jung, H., Katayama, I., & Skemer, P. (2008).
Geodynamic significance of seismic anisotropy of the upper mantle: New insights from laboratory studies. Annu. Rev. Earth Planet. Sci., 36, 59–95.
