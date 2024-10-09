# Mesh refinement criteria

Despite research since the mid-1980s, it isn't completely clear how to
refine meshes for complex situations like the ones modeled by
ASPECT. The basic problem is that mesh refinement
criteria either can refine based on some variable such as the temperature, the
pressure, the velocity, or a compositional field, but that oftentimes this by
itself is not quite what one wants. For example, we know that Earth has
discontinuities, e.g., at 440km and 610km depth. In these places, densities
and other material properties suddenly change. Their resolution in computation
models is important as we know that they affect convection patterns. At the
same time, there is only a small effect on the primary variables in a
computation &ndash; maybe a jump in the second or third derivative, for
example, but not a discontinuity that would be clear to see. As a consequence,
automatic refinement criteria do not always refine these interfaces as well as
necessary.

To alleviate this, ASPECT has plugins for mesh
refinement. Through the parameters in {ref}`parameters:Mesh_20refinement/Strategy`, one can select
when to refine but also which refinement criteria should be used and how they
should be combined if multiple refinement criteria are selected. Furthermore,
through the usual plugin mechanism, one can extend the list of available mesh
refinement criteria (see the parameter "Strategy" in
{ref}`parameters:Mesh_20refinement/Strategy`). Each such plugin is responsible for producing a
vector of values (one per active cell on the current processor, though only
those values for cells that the current processor owns are used) with an
indicator of how badly this cell needs to be refined: Large values mean that
the cell should be refined, small values that the cell may be coarsened away.

To implement a new mesh refinement criterion, you need to overload the
`aspect::MeshRefinement::Interface` class and use the
`ASPECT_REGISTER_MESH_REFINEMENT_CRITERION` macro to register your new class.
The implementation of the new class should be in namespace
`aspect::MeshRefinement`.

The primary function of this interface, `execute()`, computes the set of refinement criteria (one per
cell) and returns through its function argument. Typical examples can be found in
the existing implementations in the `source/mesh_refinement` directory. As
usual, your termination criterion implementation will likely need to be
derived from the `SimulatorAccess` class to get access to the current state of the
simulation.
