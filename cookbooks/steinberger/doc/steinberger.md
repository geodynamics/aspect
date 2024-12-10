(sec:cookbooks:steinberger)=
# Convection using a pressure--temperature look-up table and the rheology of Steinberger and Calderwood (2006)

*This section was contributed by Juliane Dannberg and Ren&eacute;
Gassm&ouml;ller.*

In this cookbook we will go one step further from the last one and set up a
fully compressible mantle convection model using the projected density
approximation (where the density is interpolated onto the finite element grid
to compute the density gradients in the mass conservation equation rather than
approximating these gradients using a reference profile or
temperature/pressure derivatives of the density, see {cite:t}`gassmoller:etal:2020`).
To compute the material properties, we read in a look-up table of
material properties in dependence of temperature and pressure, originally
computed using a mineral physics software (in this case, Perple_X, {cite}`connolly:2005`).
The table is based on the thermodynamic database of {cite:t}`stixrude:lithgow-bertelloni:2011`
and a pyrolitic composition {cite}`ringwood:irifune:1988`.
Compared to a 1D profile, a temperature-pressure look-up table has the
advantage that material properties are accurate not only around one reference
adiabat, but also for strongly deviating pressures and temperatures. This is
particularly important at phase transitions, because their depth depends on
the temperature and pressure.

This cookbook also demonstrates how to read in a viscosity profile from a data
file. Specifically, we use the profile and lateral viscosity variations due to
temperature from {cite:t}`steinberger:calderwood:2006`, which are based on mineral
physics constraints and surface observations.

In addition, this cookbook shows the use of periodic boundary conditions.

## Geometry and periodic boundaries.

The model setup is a quarter spherical shell with periodic side boundaries.
The inner and outer radius are 3481&nbsp;km and 6371&nbsp;km, respectively, so
that the mantle is 2900&nbsp;km deep. In the same section of the input file,
we also need to specify that the model should have periodic boundaries in
angular ($\phi$) direction:

```{literalinclude} geometry.part.prm
```

Both the top and bottom boundaries allow for free slip. Because the model has
periodic side boundary conditions and free slip boundaries at top and bottom,
the amount of rigid-body rotation in $\phi$ direction is not constrained. In
other words: There is no unique solution. ASPECT can remove this nullspace from the
model (see {ref}`sec:methods:nullspace-removal`). Here, we do this by setting
the net rotation to zero:

```{literalinclude} nullspace.part.prm
```

The temperature is fixed to 273&nbsp;K at the top and 3773&nbsp;K at the
bottom boundary. The initial temperature model consists of an adiabatic
profile, thermal boundary layers at the surface and the core-mantle boundary,
and a small harmonic perturbation to initiate convection. The gravity profile
in the model is based on PREM.

## The equation of state

To use material properties from a temperature--pressure look-up table, we
use the Steinberger material model. We have to specify the path to the
directory where all the data files we want to use for this model are stored.
This includes the files for the viscosity profile, the lateral viscosity
variations due to temperature, and all material files containing look-up
tables computed by mineral physics software. In addition we have to specify
the names of these files. In our case, we only have one of these look-up
tables, because we only have one composition: pyrolite. But in principle, the
material model can use several compositions with one look-up table for each.
For intermediate composition values, material properties will then be averaged
based on the mass/volume fractions of the individual compositions.

In addition, there are a few options we can select about how these look-up
tables should be used: We can decide between interpolating between data points
in the lookup table based on the pressure and temperature at the point we need
the material properties for, or we can simply take the value from the table
that is closest. In our case, we choose the bilinear interpolation because it
is more accurate. Second, we can decide how latent heat should be computed:
from the thermal expansivity and specific heat, or from the enthalpy (all
three properties should be columns in the look-up table). In some cases the
look-up table contains the effective thermal expansivity and specific heat.
Using these effective properties automatically includes the latent heat
release and consumption at phase transformations in the adiabatic heating term
and the left-hand side term (change in thermal energy over time) of the energy
equation. In a case like that, we simply want to use these values without
using additional latent heat terms because latent heat is already included
automatically when using the properties from the look-up table. If the look-up
table contains thermal expansivity and specific heat without the effect of
phase transitions, then ASPECT can compute latent heat effects based on the pressure
and temperature derivatives of the specific enthalpy (using the approach of
{cite:t}`nakagawa:etal:2009`). In our case, we simply do not include latent heat at
all in our model. So the look-up table is computed without latent heat
effects, and we set the "Latent heat" parameter to `false`.

```{literalinclude} lookup.part.prm
```

In an actual research application, it would be appropriate (and consistent
with the projected density approximation, or any other compressible
approximation) to compute latent heat instead of neglecting it as we do. This
often leads to numerical instabilities that one typically addresses by
ensuring that either the resolution is fine enough so that each phase
transitions is resolved by several mesh cells, or the energy equations needs
to be solved for entropy instead of pressure (which is an option available in
ASPECT; in this case, the look-up table needs to be given in terms of entropy and
pressure).

## The look-up table format

The format of these look-up tables is described in the documentation of the
[aspect::MaterialModel::MaterialUtilities::Lookup::MaterialLookup] class.
Two different formats are currently supported: Perple_X and HeFESTo. The
format needs to be selected in the input file, and each format has a specific
header and needs to be structured in a specific way. The paragraph below
explains how to structure a Perple_X file. This file format is the default, it
is also the more flexible format and it is what is used in this cookbook.
Since the only requirements for the format are the header and the order of
some of the columns, files created with other mineral physics software can
also be converted to this format.

The Perple_X header contains the following in the first 13 lines:

1.  The Perple_X version,

2.  the name of the data table,

3.  the dimensions of the data table (for example, for a table with one
    dimension being pressure, the other temperature, this would be 2),

4.  the variable in the first dimension (this either needs to be `T(K)` for
    temperature, or `P(bar)` for pressure),

5.  the minimum value of this variable,

6.  the increments this variable will be increased with in the table,

7.  the number of different values of this variable the table contains,

8.  the second variable,

9.  the minimum value of this second variable,

10. the increments this second variable will be increased with in the table,

11. the number of different values of this second variable the table contains,

12. the number of material properties in the table, and finally,

13. the names of the columns.

The first two columns need to be the pressure and temperature (in any order).
The other required column names are: `rho,kg/m3` (for the density),
`alpha,1/K` (for the thermal expansivity), `cp,J/K/kg` (for the specific
heat), `vp,km/s` (for the P-wave velocity), `vs,km/s` (for the S-wave
velocity)), `h,J/kg` (for the specific enthalpy). Optionally, the file can
contain columns with the name `phase` (to read in the name of the dominant
phase), and columns named `vol_fraction_` and the name of a phase after the
second underscore (to read in volume fractions of different phases). As an
example, the header of the table used in this cookbook is given below:

    |6.6.6
    PYR-Ringwood88_2.tab
               2
    T(K)
       400.00360000000001
       19.999960000020000
             181
    P(bar)
       15001.334999999999
       5114.9322988556905
             262
               8
    T(K)           P(bar)         rho,kg/m3      alpha,1/K      cp,J/K/kg      vp,km/s        vs,km/s        h,J/kg

Below this header, the table contains the actual data values, using one column
for each of the property names given in the last line of the header. It is
also useful to know that ASPECT does not actually read in the values of the pressures
and temperatures in the first two columns, but instead uses the minimum,
increment, and number of values parameters given in the header, assuming a
uniform step size. The first column is always assumed to be the inner loop
(i.e., it needs to increase first while the second column stays constant).

## The rheology.

The rheology of this model consists of two parts: The viscosity profile, and
the lateral variations due to temperature. For each of these, we need to read
in a data file. In this example, we use files that are based on
{cite:t}`steinberger:calderwood:2006` for both. The viscosity profile is based on mineral
physics and surface constrains, and the lateral viscosity variations use an
Arrhenius law with a depth-dependent activation enthalpy. For more details and
a derivation, see {cite:t}`steinberger:calderwood:2006`.

Other rheology models can be used by reading in different files. The
formatting of these files is the following: The radial viscosity file contains
two columns, where the first is the viscosity in Pa&nbsp;s, and the second is
the depth in km (note that this is an exception to the usual ASPECT convention of
using SI units). The lateral viscosity file also contains two columns, the
first being the activation enthalpy divided by the gas constant and the
nondimensional stress exponent (which is 1 for diffusion creep/in the lower
mantle, and 3.5 for dislocation creep/in the upper mantle and transition zone
in the model of {cite:t}`steinberger:calderwood:2006`). The second column is
depth, again in km. Both parts are combined to compute the viscosity in the
following way:
```{math}
\eta = \eta_\text{rad} \exp{ \left( -\frac{V_\text{lat} \Delta T}{T T_\text{ref}} \right)} ,
```
where $\eta_\text{rad}$ is the value from the radial viscosity file,
$V_\text{lat}$ is the value from the lateral viscosity file, $T$ is
temperature, $T_\text{ref}$ is the reference temperature profile, and
$\Delta T$ is the deviation from the reference temperature profile.

This reference profile can be chosen in several different ways: On the one hand,
it can be chosen as the laterally averaged temperature (and in this case, a
number of depth slices for this lateral averaging can be specified as well).
This is the original formulation of {cite:t}`steinberger:calderwood:2006`, and the
default of the material model. On the other hand, the adiabatic temperature
profile can be chosen as the reference. However, the radial profile needs to
be adapted based on how this reference temperature is chosen. If the reference
profile uses the laterally averaged temperature, then the radial profile needs
to include a high viscosity in the lithosphere (where it is cold), and a low
viscosity near the core-mantle boundary (where it is warm). If the reference
profile is the adiabatic profile, then the temperature will deviate from this
reference in the top and bottom thermal boundary layers already, leading to
changes in viscosity. So in this case, the radial profile should not include
these boundary layers (because otherwise we would compute their effect twice).
This option allows the viscosity in the boundary layers to develop based on
the temperature in the model, which is why we choose it for this cookbook.

The default data directory already contains several radial viscosity files.
The file
[data/material-model/steinberger/radial-visc.txt](https://www.github.com/geodynamics/aspect/blob/main/data/material-model/steinberger/radial-visc.txt) is the original
{cite:t}`steinberger:calderwood:2006` profile (with an
interpolation between the original discrete layers) and for use with the
laterally averaged temperature. The file
[data/material-model/steinberger/radial-visc-simple.txt](https://www.github.com/geodynamics/aspect/blob/main/data/material-model/steinberger/radial-visc-simple.txt) is for use with the
adiabatic profile. To illustrate the difference, the content of both files is
plotted in {numref}`fig:steinberger-viscosity1` and {numref}`fig:steinberger-viscosity2`.


```{figure-md} fig:steinberger-viscosity1
<img src="radial-visc.svg" style="width:48.0%" />

Viscosity profile based on the original {cite:t}`steinberger:calderwood:2006` formulation, intended for use with a temperature-dependence of viscosity based on the laterally averaged temperature.
```

```{figure-md} fig:steinberger-viscosity2
<img src="radial-visc-simple.svg" style="width:48.0%" />

Modified viscosity profile without boundary layers, intended for use with a temperature dependence of viscosity based on an adiabatic temperature profile.
```

In addition, there are two different ways these viscosity profiles are
represented in {cite:t}`steinberger:calderwood:2006`. They start out with
non-optimized normalized viscosity profiles (as in their Figure 4), and then
optimize them based on observational data. During this optimization, the different
layers that make up the profiles (usually lithosphere, upper mantle, transition
zone, and lower mantle) are shifted left or right relative to each other.
The figures that show these optimized profiles in the paper represent them as
piece-wise constant values within 22 layers, and the profiles shown in Figures
{numref}`fig:steinberger-viscosity1` and {numref}`fig:steinberger-viscosity2`
are interpolations between these 22 layers. However, we can also represent
the profile as a variation of the initial normalized viscosity profile, with
the four different layers being shifted relative to each other based on the
optimization. This viscosity profile is given in
[data/material-model/steinberger/radial-visc-continuous.txt](https://www.github.com/geodynamics/aspect/blob/main/data/material-model/steinberger/radial-visc-continuous.txt).
The image below illustrates the difference compared to
[data/material-model/steinberger/radial-visc-simple.txt](https://www.github.com/geodynamics/aspect/blob/main/data/material-model/steinberger/radial-visc-simple.txt), the profile based
on interpolation of the discrete layers.
The profile with the continuously varying viscosity within 4 layers is the more faithful
representation of the results of {cite:t}`steinberger:calderwood:2006`, and we therefore
recommend its use in models that use their viscosity profile.

```{figure-md} fig:steinberger-viscosity3
<img src="comparison.svg" style="width:48.0%" />

Comparison of the viscosity profiles without boundary layers, showing the interpolation of the piece-wise constant profile [data/material-model/steinberger/radial-visc-simple.txt](https://www.github.com/geodynamics/aspect/blob/main/data/material-model/steinberger/radial-visc-simple.txt) in green and the profile that is continuous within 4 layers, but features jumps between these layers [data/material-model/steinberger/radial-visc-continuous.txt](https://www.github.com/geodynamics/aspect/blob/main/data/material-model/steinberger/radial-visc-continuous.txt) in purple.
```

To change from the interpolated profile (which is the one being used in the
input file of this cookbook) requires changing the `Radial viscosity file name`
parameter in the Steinberger material model to the desired value.

```{literalinclude} rheology2.part.prm
```

In order to improve solver convergence, the material model has additional
parameters that allow it to limit the viscosity variations. Because of the
resolution in this cookbook we limit the lateral viscosity variations to three
orders of magnitude in both directions (for a total of six orders of
magnitude), and we additionally limit the overall viscosity between
$10^{20}$&nbsp;Pa&nbsp;s and $5 \times 10^{23}$&nbsp;Pa&nbsp;s. This allows
the features of the flow field to be resolved.

```{literalinclude} rheology.part.prm
```

In the Earth, we would expect higher viscosities in the lithosphere and lower
viscosities in plumes and near the core-mantle boundaries. This type of
viscosity formulation is appropriate for global convection models. However, it
does not approximate lithospheric deformation well. The model only accounts
for diffusion creep, so the lithosphere has a high viscosity and forms a
stagnant lid on top of the sublithospheric mantle. In order to achieve more
realistic subduction in a model like this, one would have to either prescribe
plate velocities at the surface (forcing plates to subduct) or take into
account plastic yielding (so that the lithosphere can break).

:::{note}
If the model takes too long to run, increase the minimum viscosity.
:::

## The projected density approximation

Since our model is compressible, the most accurate way to solve the mass
conservation equation implemented in ASPECT is to use the 'projected density
approximation.' This way, ASPECT will compute the density gradients in the mass
conservation directly from the density field (interpolated onto the finite
element grid) rather than approximating it with a reference profile or
temperature/pressure derivatives of the density.

To use the projected density approximation, we need to specify it as the form
of the equations we want to use, and we need to provide a field that the
density values can be interpolated on. The first part is handled in the
'Formulation' section of the input file. This is where we can
select the projected density approximation as the formulation we want to use
for the mass conservation equation. The temperature equation uses the real
density (rather than a reference profile) as well.

To allow for the interpolation, we create a compositional field that we call
'density_field.' We assign the field the type
'density,' so that ASPECT knows that this is the field it should use to
compute the density gradient required to solve the equations. ASPECT does not need to
solve an equation for this field, it only needs to interpolate the density
values onto it. This is covered by the compositional field method
'prescribed field.' For fields of this type, the material model
provides the values that should be interpolated onto the field.

```{literalinclude} projected_density.part.prm
```

The complete input file can be found in
[cookbooks/steinberger/doc/steinberger.prm](https://www.github.com/geodynamics/aspect/blob/main/cookbooks/steinberger/steinberger.prm).

## Results

We run the model for 300 million years. Over the time of the model evolution,
some plumes rise and spread beneath the base of the lithosphere, and some cold
downwellings detach from the base of the lithosphere. The temperature at the
end of the model run and some of the material properties are shown in
{numref}`fig:steinberger-end-state`.


```{figure-md} fig:steinberger-end-state
<img src="endstate.png" style="width:96.0%" />

 End state of the model. From left to right and top to bottom: Temperature, viscosity, density, and specific heat capacity.
```

For comparison, we also show the end state of the model with the continuous
viscosity profile.
```{figure-md} fig:steinberger-end-state-2
<img src="endstate_continuous.png" style="width:96.0%" />

 End state of the model with continuous viscosity profile. Left: temperature. Right: viscosity.
```
