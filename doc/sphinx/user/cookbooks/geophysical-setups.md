(sec:cookbooks:geophysical-setups)=
# Geophysical setups

Having gone through the ways in which one can set up problems in rectangular
geometries, let us now move on to situations that are directed more towards
the kinds of things we want to use <span class="smallcaps">ASPECT</span> for:
the simulation of convection in the rocky mantles of planets or other
celestial bodies.

To this end, we need to go through the list of issues that have to be
described and that were outlined in {ref}`sec:cookbooks:overview`, and address them one
by one:

-   *What internal forces act on the medium (the equation)?* This may in fact
    be the most difficult to answer part of it all. The real material in
    Earth's mantle is certainly no Newtonian fluid where the stress is a
    linear function of the strain with a proportionality constant (the
    viscosity) $\eta$ that only depends on the temperature. Rather, the real
    viscosity almost surely also depends on the pressure and the strain rate.
    Because the issue is complicated and the exact material model not entirely
    clear, for the next few subsections we will therefore ignore the issue and
    start with just using the &ldquo;simple&rdquo; material model where the
    viscosity is constant and most other coefficients depend at most on the
    temperature.

-   *What external forces do we have (the right hand side)* There are of
    course other issues: for example, should the model include terms that
    describe shear heating? Should it be compressible? Adiabatic heating due
    to compression? Most of the terms that pertain to these questions appear
    on the right hand sides of the equations, though some (such as the
    compressibility) also affect the differential operators on the left.
    Either way, for the moment, let us just go with the simplest models and
    come back to the more advanced questions in later examples.

    One right hand side that will certainly be there is that due to
    gravitational acceleration. To first order, within the mantle gravity
    points radially inward and has a roughly constant magnitude. In reality,
    of course, the strength and direction of gravity depends on the
    distribution and density of materials in Earth &ndash; and, consequently,
    on the solution of the model at every time step. We will discuss some of
    the associated issues in the examples below.

-   *What is the domain (geometry)?* This question is easier to answer. To
    first order, the domains we want to simulate are spherical shells, and to
    second order ellipsoid shells that can be obtained by considering the
    isopotential surface of the gravity field of a homogeneous, rotating
    fluid. A more accurate description is of course the geoid for which
    several parameterizations are available. A complication arises if we ask
    whether we want to include the mostly rigid crust in the domain and simply
    assume that it is part of the convecting mantle, albeit a rather viscous
    part due to its low temperature and the low pressure there, or whether we
    want to truncate the computation at the asthenosphere.

-   *What happens at the boundary for each variable involved (boundary
    conditions)?* The mantle has two boundaries: at the bottom where it
    contacts the outer core and at the top where it either touches the air or,
    depending on the outcome of the discussion of the previous question, where
    it contacts the lithospheric crust. At the bottom, a very good
    approximation of what is happening is certainly to assume that the
    velocity field is tangential (i.e., horizontal) and without friction
    forces due to the very low viscosity of the liquid metal in the outer
    core. Similarly, we can assume that the outer core is well mixed and at a
    constant temperature. At the top boundary, the situation is slightly more
    complex because in reality the boundary is not fixed but also allows
    vertical movement. If we ignore this, we can assume free tangential flow
    at the surface or, if we want, prescribe the tangential velocity as
    inferred from plate motion models. <span class="smallcaps">ASPECT</span>
    has a plugin that allows to query this kind of information from the
    `GPlates` program.

-   *How did it look at the beginning (initial conditions)?* This is of course
    a trick question. Convection in the mantle of earth-like planets did not
    start with a concrete initial temperature distribution when the mantle was
    already fully formed. Rather, convection already happened when primordial
    material was still separating into mantle and core. As a consequence, for
    models that only simulate convection using mantle-like geometries and
    materials, no physically reasonable initial conditions are possible that
    date back to the beginning of Earth. On the other hand, recall that we
    only need initial conditions for the temperature (and, if necessary,
    compositional fields). Thus, if we have a temperature profile at a given
    time, for example one inferred from seismic data at the current time, then
    we can use these as the starting point of a simulation.

This discussion shows that there are in fact many pieces with which one can
play and for which the answers are in fact not always clear. We will address
some of them in the cookbooks below. Recall in the descriptions we use in the
input files that <span class="smallcaps">ASPECT</span> uses physical units,
rather than non-dimensionalizing everything. The advantage, of course, is that
we can immediately compare outputs with actual measurements. The disadvantage
is that we need to work a bit when asked for, say, the Rayleigh number of a
simulation.

:::{toctree}
cookbooks/shell_simple_2d/doc/shell_simple_2d.md
cookbooks/shell_simple_3d/doc/shell_simple_3d.md
cookbooks/shell_3d_postprocess/doc/shell_3d_postprocess.md
cookbooks/initial-condition-S20RTS/doc/initial-condition-S20RTS.md
cookbooks/gplates/doc/gplates.md
cookbooks/burnman/doc/burnman.md
cookbooks/steinberger/doc/steinberger.md
cookbooks/multicomponent_steinberger/doc/multicomponent_steinberger.md
cookbooks/morency_doin_2004/doc/morency_doin_2004.md
cookbooks/crustal_deformation/doc/crustal_deformation.md
cookbooks/continental_extension/doc/continental_extension.md
cookbooks/inner_core_convection/doc/inner_core_convection.md
cookbooks/lower_crustal_flow/doc/lower_crustal_flow.md
cookbooks/global_melt/doc/global_melt.md
cookbooks/mid_ocean_ridge/doc/mid_ocean_ridge.md
cookbooks/grain_size_ridge/doc/grain_size_ridge.md
cookbooks/transform_fault_behn_2007/doc/transform_fault_behn_2007.md
cookbooks/kinematically_driven_subduction_2d/doc/kinematically_driven_subduction_2d.md
cookbooks/allken_et_al_2012_rift_interaction/doc/allken.md
cookbooks/tian_parameterization_kinematic_slab/doc/tian_parameterization_kinematic_slab.md
cookbooks/mantle_convection_with_continents_in_annulus/doc/mantle_convection_in_annulus.md
cookbooks/inclusions/doc/inclusions.md
cookbooks/subduction_initiation/doc/subduction_initiation.md
cookbooks/vankeken_subduction/doc/vankeken_subduction.md
cookbooks/2d_annulus_visualization/doc/2d_annulus_visualization.md
cookbooks/tomography_based_plate_motions/doc/tomography_based_plate_motions.md
cookbooks/future/README.md
:::
