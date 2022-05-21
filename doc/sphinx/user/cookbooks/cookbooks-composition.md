(sec:cookbooks-composition)=
# Using passive and active compositional fields

One frequently wants to track where material goes, either because one simply
wants to see where stuff ends up (e.g., to determine if a particular model
yields mixing between the lower and upper mantle) or because the material model
in fact depends not only pressure, temperature and location but also on the
mass fractions of certain chemical or other species. We will refer to the first
case as *passive* and the latter as *active* to indicate the role
of the additional quantities whose distribution we want to track. We refer to
the whole process as *compositional* since we consider quantities that
have the flavor of something that denotes the composition of the material at any
given point.
There are basically two ways to achieve this: one can advect a set of
particles ("tracers") along with the velocity field, or one can advect along a
field. In the first case, where the closest particle came from indicates the
value of the concentration at any given position. In the latter case, the
concentration(s) at any given position is simply given by the value of the
field(s) at this location.
ASPECT implements both strategies, at least to a certain degree. In this
cookbook, we will follow the route of advected fields.

:::{toctree}
cookbooks/composition_passive/doc/composition_passive.md
cookbooks/composition_active/doc/composition_active.md
cookbooks/composition-reaction/doc/composition-reaction.md
:::
