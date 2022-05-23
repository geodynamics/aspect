This is a cookbook demonstrating an application of fixing the temperature
in user specified compositional fields. The model simulates 2D subduction of an arbitrary slab
generated using the Geodynamic WorldBuilder (Fraters et al. 2019) {numref}`fig:initial-composition`,
where the temperature within the slab is fixed to the temperature it was initialized with.
The model runs for 15 Myr, and during this time the slab maintains its thermal state
while the mantle wedge convects and reaches a more realistic thermal state.

```{figure-md} fig:initial-composition
<img src="initial_composition.*" alt="Compositions defined in the fixed slab temperature cookbook"/>
```

The plugin works by resetting the temperature degrees of before each temperature solve
to the temperature that the compositional field was initialized with. This ensures that
all points within the specified field do not evolve their temperature with time, and this can be easily
toggled in the input file by setting `Prescribe internal temperature = true` in the
`Prescribed internal temperature model` subsection. In this cookbook, the slab not only has its
temperature fixed, but its location is also fixed by making the field a `static` field. This allows for
the location of the mantle wedge to also be constant, and thus allows a steady state to be reached. This
is great for those interested in modeling the dynamics of the mantle wedge. A setback of this method,
however, is that while temperature inside the slab remains fixed, the temperature within the slab
can still be advected out of the slab and into the surrounding mantle. In the case of the
mantle wedge, this is fine because the velocities in the slab point largely straight down
and away from the mantle wedge due to the slabs large negative buoyancy.

The initial temperature and compositional fields are defined by the Geodynamic WorldBuilder, and a
visco plastic material model is used to simulate slab subduction dynamics. This creates a compositie rheology incorporating deformation from brittle failure, and the following creep mechanisms: diffusion, dislocation, and Peierls. Parameters for each of these mechanisms are choosen based on laboratory experiments deforming olivine grains {cite}`HK04,mei_et_al_2010_JGR`.
```{literalinclude} ../fixed_slab_temperature.prm
    :lines: 94-143
```

The starting temperature field and final temperature field can be seen in {numref}`fig:temperature-fields`
. Initially, the temperature in the mantle is not influenced by the presence of the slab, consequently if
a model were to start from this state the thermal structure of the mantle wedge would be incorrect which
could negatively impact model results. After 15 Myr, the temperature structure has evolved due to the
convection in the mantle wedge, and has taken on a more realistic shape. The advection of the temperature
from within the slab can be seen on the bottom side of the slab in the mantle, where the temperature
contours bulge out.

```{figure-md} fig:temperature-fields
<img src="temperature_fields.*" alt="Temperature fields of subducting slab."/>

Temperature fields at start of model time (left), and after 15 Myr (right). Vectors show velocity
direction and colour shows velocity magnitude.
```