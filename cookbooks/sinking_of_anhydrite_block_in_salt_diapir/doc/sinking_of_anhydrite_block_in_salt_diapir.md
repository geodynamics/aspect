```{tags}
category:cookbook
```

(sec:cookbooks:sink-anydr)=
# Sinking of anhydrite blocks within a Newtonian salt diapir

*This section was contributed by Cedric Thieulot.*

The setup for this experiment originates in {cite:t}`burchardt:etal:2012a`.
Similar experiments are conducted by the same authors in
{cite:t}`burchardt:etal:2011` and {cite:t}`burchardt:etal:2012b`.

Based on the principal structural configurations of the Main Anhydrite segments
within the Gorleben diapir, as shown in {numref}`fig:sinking-anhydrite-diapir`,
the authors set up three series of models.

```{figure-md} fig:sinking-anhydrite-diapir
<img src="diapir.*" width="80%" />

NW–SE cross-section through the Gorleben diapir.
Taken from {cite:t}`burchardt:etal:2012a`.
```

Each model analyses the deformation associated with the gravity-driven sinking
of one anhydrite block. The block geometry was simplified to rectangular shapes
of different aspect ratios which vary from 10:1 to 1:10,
with a 1:1 ratio corresponding to a $100 \times 100~\text{ m}$ block as shown in
the fifth panel of {numref}`fig:sinking-anhydrite-setups`.

```{figure-md} fig:sinking-anhydrite-setups
<img src="setups.*" width="90%" />

Scaled sketches of model set-ups.
Taken from {cite:t}`burchardt:etal:2012a`.
```

In the article the authors consider three cases A,B,C with potentially different
domain sizes or viscosity stratification. We focus here on experiment C with a domain
of size $L_x=2500~\text{ m}$ and $L_y=5000~\text{ m}$. The gravity is not specified
in the paper but it is set to $g=9.81~\text{ m s}^{-2}$. Boundary conditions are
free slip on all four sides. Simulations are run up to about $500~\text{ ka}$.
There are three materials in the domain as shown in {numref}`fig:sinking-anhydrite-setup`.
The full input file that contains these modifications and that was used for the simulations
can be found at [cookbooks/sinking_anhydrite_block.prm](https://github.com/geodynamics/aspect/blob/main/cookbooks/sinking_of_anhydrite_block_in_salt_diapir/sinking_of_anhydrite_block_in_salt_diapir_particle_in_cell.prm).

```{figure-md} fig:sinking-anhydrite-setup
<img src="setup.*" width="80%" />

Geometry, density and viscosity values of the three compositions in the case
of a block of aspect ratio 1:1, i.e. $b_x=b_y=100~\text{ m}$.
Note that the top edge of the block is always at distance $h=100~\text{ m}$
below the surface at the beginning of the simulation.
```

In all three articles the authors use the Finite Differences code FDCON from {cite:t}`weinberg:schmeling:1992`.
In the current case the original resolution of the experiment is $200 \times 400$ cells,
i.e. cells of $12.5 \times 12.5~\text{ m}$ in size. The particles are located
every 10 to $12.5~\text{ m}$ in vertical and horizontal direction, i.e. about 100,000 particles in total.
We here make use of adaptive mesh refinement with a criterion based on the compositions interface and
seed the system with 2,000,000 randomly distributed particles. Elements of the finest are
then $5 \times 5~\text{ m}$ in size.
Unfortunately the authors fail to mention the type of averaging that they use for the density and
viscosity carried by the particles. We here choose the harmonic one, keeping in mind that
it may affect the dynamics of the system somewhat (e.g. {cite:t}`schmeling:etal:2008`).

Rather importantly, the authors add that their models are based on a number of
assumptions and simplifications:

-   All materials are homogeneous and isotropic, neglecting any stratigraphic
heterogeneities within the salt formations or the anhydrite.

-   The materials used are incompressible and entirely viscous, that is,
no elastic behavior of, for example, the anhydrite is enabled.

-   The salt rheology in the models is Newtonian. However, salt rheology is a
complex product of, for example, composition, grain size, fluid content, temperature and
strain rate (e.g. {cite:t}`urai:etal:1986`, {cite:t}`vankeken:etal:1993`, {cite:t}`jackson:etal:1994`).

-   The models are isothermal, neglecting temperature effects on the rheology.

-   Limitations regarding geometry include the simplified, rectangular shape of
the anhydrite blocks with thicknesses of 100 m (instead of 70 m in case of the Main
Anhydrite) and the perfectly straight interface between the two salt
types. Hence, pre-existing deformation caused by salt ascent and
emplacement along with the entrainment of the Main Anhydrite layer is neglected.

On {numref}`fig:sinking-anhydrite-results`, the panels of Fig. 6 of the paper
are placed next to ASPECT's results at three different times.
We find a reasonable agreement between both codes, again
acknowledging that the size of the sinking block probably requires an
even higher resolution that what is here used and that material
averaging is also likely to change results.

```{figure-md} fig:sinking-anhydrite-results
<img src="results.*" width="90%" />

Results obtained at three different times: a) $t=104~\text{ ka}$;
b) $t=228~\text{ ka}$; c) $t=352~\text{ ka}$. Blue-green panels are
from the original paper while black \& white plots are obtained with
aspect. The last panel on the right shows the mesh and highlights
Adaptive Mesh Refinement working optimally to capture the materials interface.
```
