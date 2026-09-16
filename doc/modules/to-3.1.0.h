/**
 * @page changes_between_3.0.0_and_3.1.0 Changes between version 3.0.0 and version 3.1.0
 *
 * <p> This is the list of changes made after the release of ASPECT version
 * 3.0.0 for version 3.1.0. All entries are signed with the names of the author.
 * </p>
 *
 * <ol>
 *
 * <li> New: The mesh deformation handler can now apply the initial mesh deformation
 * in a configurable number of substeps (new parameter ``Number of initial mesh
 * deformation substeps'', default 1) to avoid inverted cells for steep initial
 * topography.
 * <br>
 * (Timo Heister, 2026/09/10)
 *
 * <li> New: The Landlab mesh deformation plugin can now provide the model time as
 * the boundary value for a compositional field named ``sediment_age'', enabling
 * the age of deposited sediments to be tracked in ASPECT models.
 * <br>
 * (Timo Heister, 2026/09/02)
 *
 * <li>  Changed: The output directory of the gravity point values
 *  postprocessor has been renamed from 'output_gravity' to
 *  'gravity_point_values' to avoid name conflicts with
 *  other postprocessors and output directories.
 *  <br>
 *  (Rene Gassmoeller, 2026/09/01)
 *
 * <li> Added: a new particle property ``particle generation time`` that stores the model time at which each particle is generated.
 * <br>
 * (Ranpeng Li 2026/07/27)
 *
 * <li>  New: There is now a Boundary Composition plugin ``mesh deformation''
 *  that queries the active mesh deformation plugins for
 *  compositional field values to be prescribed on the boundaries.
 *  By default, the mesh deformation plugins return 0.
 *  <br>
 *  (Anne Glerum 2026/08/24)
 *
 * <li>  Changed: The parameter ``Model name'' has been removed
 *  from the Boundary Composition and Boundary
 *  Temperature subsections after being deprecated for
 *  a long time. In addition, using ``Model name'' will now
 *  throw an error when used in the Initial Composition and
 *  Initial Temperature subsections. In all cases, ``List of
 *  model names'' should be used instead. A script is available
 *  for updating existing prm files in the contrib/utilities/
 *  folder.
 *  <br>
 *  (Anne Glerum 2026/08/07)
 *
 * <li> Added: A mesh deformation plugin which couples to the landscape
 * evolution library Landlab to deform the surface of an ASPECT model.
 * <br>
 * (Timo Heister, Daniel Douglas, 2026/07/26)
 *
 * <li> New: Added a new particle property called 'general composition reaction' that
 * allows compositional fields to be tracked on particles and update them using
 * the 'reaction terms' returned by the material model. While this particle
 * property allows different fields to be tracked on different particle managers, it is
 * recommended to have all fields on a single manager when dealing with sharp interfaces.
 * <br>
 * (Arijit Chakraborty and Rene Gassmoeller, 2026/09/09)
 *
 * <li>  Fixed: The Initial Composition and Initial
 *  Temperature plugins ``Ascii data layered'' now
 *  use their own subsection header ``Ascii data
 *  layered'' instead of ``Ascii data model''.
 *  Therefore the plugins now show up in the manual
 *  with their own subsection.
 *  <br>
 *  (Anne Glerum 2026/08/01)
 *
 * <li> New: Added a ``particle information`` postprocessor that prints the property
 * names of every particle world at time zero.
 * <br>
 * (Timo Heister, 2026/07/30)
 *
 * <li> Changed: Models with melt transport can now use discontinuous elements for
 * compositional fields, with two exceptions: the porosity field and fields
 * that are advected with the melt velocity still require continuous elements.
 * Compositional fields that are advected as finite element fields with the
 * solid velocity now use the same discontinuous Galerkin (DG) face terms as
 * models without melt transport.
 * <br>
 * (Ryan Stoner, 2026/07/30)
 *
 * <li> New: Initial topography can now be prescribed through the
 * Geodynamic World Builder. The plugin reads the surface topography
 * and maximum expected elevation from the configured World Builder file.
 * <br>
 * (Michael Pons, 2026/07/29)
 *
 * <li> The heating model 'Tidal Heating' with the option 'latitudinal variation' now works for 2-dimensional spherical geometries.
 * The model assumes that the 2D geometry represents the simplified latitudinal variation from Nimmo et al. (2007).
 * <br>
 * (Hyunseong Kim 2026/07/30)
 *
 * <li> Added: Swapped std::min(std::max(val, low), high) for std::clamp(val, low, high)
 * <br>
 * (Buchanan Kerswell, 2026/07/30)
 *
 * <li> New: Added a ``finite element information`` postprocessor that prints the
 * finite element space of every solution variable at time zero. For
 * compositional fields, the output includes the field index, name, and type.
 * <br>
 * (Timo Heister, 2026/07/29)
 *
 * <li> Changed: The melt fraction visualization postprocessor has been restructured to allow
 * more general usage in the future. To compute the melt fraction of pyroxenite,
 * "pyroxenite" must now be listed not only in the compositional fields subsection but also
 * in the melt fraction visualization subsection.
 * <br>
 * (Qianyi Lu, 2026/07/29)
 *
 * <li> Fixed: Mesh deformation now uses geometry-specific periodicity constraints,
 * including the rotation of vector components in curved geometries. This fixes
 * free surfaces in phi-periodic spherical shell geometries, which previously
 * failed while constructing constraints.
 * <br>
 * (Pons Michaël, 2026/07/29)
 *
 * <li> The gravity model 'Radial with tidal potential' now works for 2-dimensional spherical geometries.
 * The model assumes that the 2D geometry represents the equatorial plane.
 * <br>
 * (Hyunseong Kim 2026/07/29)
 *
 * <li> New: The matrix-free Stokes solver with global coarsening geometric multigrid
 * (GMG-GC) now supports Q1-projected viscosity.
 * <br>
 * (Timo Heister, 2026/07/29)
 *
 * <li> Added: 'Coordinate system' to prescribed stokes solution.
 * <br>
 * (Buchanan Kerswell, 2026/07/29)
 *
 * <li> New: Added a new subsection for Initial composition
 * particle property. It can now be specified which
 * compositional fields are tracked on which particle manager.
 * <br>
 * (Arijit Chakraborty, 2026/07/29)
 *
 * <li> New: Prescribed dilation plugin system and a first implementation (function plugin).
 * The feature is switched on by selecting any model in subsection Prescribed dilation (set List of model names = ...).
 * <br>
 * (Alexandr Dizov, 2026/07/28)
 *
 * <li> New: The matrix-free Stokes solver with global coarsening geometric multigrid (GMG-GC) now supports computations without material averaging and with project to Q1 averaging.
 * <br>
 * (Timo Heister, 2026/07/28)
 *
 * <li> New: Add a treatment for a thin, sub-grid resolution layer by assuming it weakens an interface between two compositions by a constant factor. Although possible on other compositions, this weakening should primarily be used on fields or DG fields where interfaces are sharp.
 * <br>
 * (Ryan Stoner, 2026/07/28)
 *
 * <li> Added: There is now a new class which allows the creation of
 * ScratchSpaces, which are reusuable memory blocks. These can
 * be useful in function which are called in a loop and need
 * the same size temporary vectors.
 * <br>
 * (Menno Fraters and Wolfgang Bangerth, 2026/07/28)
 *
 * <li> Fixed: There was a factor 2 missing in the implementation of the viscous dissipation post processor.
 * <br>
 * (Cedric Thieulot, 2026/07/28)
 *
 * <li> Added: 'Coordinate system' to particles 'Function' class and mesh deformation
 * 'BoundaryFunction' class.
 * <br>
 * (Buchanan Kerswell, 2026/07/28)
 *
 * <li> New: Added an option to advect different particle worlds
 * with different chosen velocity. The current velocity
 * choices are solid and fluid velocity.
 * <br>
 * (Arijit Chakraborty, 2026/07/28)
 *
 * <li> New: Integration tests can now be organized in explicitly configured
 * subdirectories (categories). CTest test names and test dependencies use paths
 * relative to the ``tests/`` directory.
 * <br>
 * (Timo Heister, 2026/07/27)
 *
 * <li> Changed: In the crust and lithosphere formation reaction model, the upwelling angle has been changed from a hardcoded value (30 degrees from the horizontal direction) to a user-defined input parameter.
 * Added: An option to select whether the harzburgite profile in the lithosphere is constant or decreases linearly with depth.
 * <br>
 * (Ranpeng Li 2026/07/27)
 *
 * <li> Added: A benchmark comparing ASPECT to the analytical solution of a
 * line load on a semi-infinite elastic half-space.
 * <br>
 * (Daniel Douglas, Cedric Thieulot, Prajakta Mohite, John Naliboff, 2026/07/27)
 *
 * <li> Fixed: A bug in adiabatic/compute_profile.cc, adiabatic/compute_entropy_profile.cc,
 * and cookbooks/tomography_based_plate_motions/plugins/reference_profile.cc when using
 * the get_property() function. When asking for a value at a depth very close to a
 * calculated profile point, the code previously always returned the value from a point
 * at shallower depth. This is now fixed and returns the value at the correct point.
 *
 * Added: A new test that computes an adiabatic profile using a pressure-temperature lookup table.
 * <br>
 * (Qianyi Lu and Ranpeng Li, 2026/07/26)
 *
 * <li> Added: A parallel unstructured interface class for mesh
 * deformation plugins for coupling ASPECT to external tools
 * that utilize independent meshes.
 * <br>
 * (Timo Heister, Wolfgang Bangerth, Daniel Douglas, 2026/07/26)
 *
 * <li> New: The online documentation now includes the Doxygen API reference, with
 * links between the Sphinx user and developer documentation and the documented
 * source code, see https://aspect-documentation.readthedocs.io/en/latest/doxygen/index.html.
 * <br>
 * (Timo Heister, 2026/07/25)
 *
 * <li> Fixed: We no longer rely on std::isnan() for deciding when to run postprocessors. This logic would break with compilers and -O3 optimization levels.
 * <br>
 * (Timo Heister, 2026/07/24)
 *
 * <li>  Added: Checkpointing can now also be
 *  requested at specific model times, in
 *  combination with, or without, checkpointing
 *  after a specific amount of wall time or
 *  number of time steps.
 *  <br>
 *  (Anne Glerum 2026/07/24)
 *
 * <li> New: The entropy method is now stabilized by two additions:
 * the heat capacity lookup is now interpolated in log space,
 * and the heat capacity is also optionally limited,
 * asymptotically approaching a maximum value at high exact values.
 * This particularly stabilizes problems dominated by thermal
 * diffusion involving near-univariant reactions.
 * <br>
 * (Bob Myhill, Ranpeng Li, 2026/07/24)
 *
 * <li> New: In nearly a thousand places, the declaration of member variables of plugin classes are now annotated with the run-time parameter from the input file they are initialized with.
 * <br>
 * (Wolfgang Bangerth, 2026/07/24)
 *
 * <li> Added: Implemented kinetic models are interface-controlled growth and
 * eutectoid decomposition from Cahn (1956; https://doi.org/10.1016/0001-6160(56)90041-4).
 * <br>
 * (Buchanan Kerswell, 2026/07/23)
 *
 * <li> Fixed: Some of the particle generator plugins had no means of serializing themselves. As a consequence, they lose state between checkpointing and resuming. This is now fixed.
 * <br>
 * (Wolfgang Bangerth, 2026/07/23)
 *
 * <li> Fixed: Some of the termination criteria and particle properties plugins had no means of serializing themselves. As a consequence, they lose state between checkpointing and resuming. This is now fixed.
 * <br>
 * (Wolfgang Bangerth, 2026/07/23)
 *
 * <li> Fixed: We forgot to serialize the termination criteria and so they lose state between checkpointing and resuming. This is now fixed.
 * <br>
 * (Wolfgang Bangerth, 2026/07/23)
 *
 * <li> Added: The darcy field advection scheme can now be advected based on the pressure gradient,
 * instead of just the buoyancy force.
 * <br>
 * (Daniel Douglas, 2026/07/22)
 *
 * <li> Updated: The Geodynamic World Builder has now been updated
 * to version 1.1.
 * <br>
 * (Menno Fraters, 2026/07/21)
 *
 * <li> Added: An option to specify chemical compositions and their entropies
 * via functions when computing adiabatic profile of the entropy model.
 * <br>
 * (Ranpeng Li, 2026/07/17)
 *
 * <li> New: Elasticity can now be combined with two-phase melt transport. The melt
 * Stokes assembler adds the stored elastic (deviatoric) shear stress as a force
 * on the solid momentum equation (after Keller et al. 2013). Setting porosity to
 * zero recovers purely viscous behavior predicted by single-phase flow (i.e. in
 * cells without melt the melt momentum operator reduces exactly to the
 * incompressible Stokes equation). In addition, visco-plastic strain weakening and the elastic-stress
 * and strain-invariant particle properties work when the visco-plastic model is
 * wrapped by the reactive fluid transport model.
 * <br>
 * (Ryan Stoner, 2026/07/14)
 *
 * <li> Added: A new composition type of "reaction progress"
 * which is used to keep track of the state of reactions,
 * for example in models considering reaction kinetics.
 * <br>
 * (Haoyuan Li, 2026/07/13)
 *
 * <li> Changed: add support for updating the InitialComposition particle property using the prescribed solution plugin.
 * <br>
 * (Haoyuan Li, 2026/07/11)
 *
 * <li> New: The global coarsening GMG Stokes solver now supports periodic boundary
 * conditions, including the rotational periodicity of spherical shells and
 * hanging nodes on the periodic boundary. This solver must be selected
 * explicitly by setting the Stokes solver type to `block GMG` and the Stokes
 * GMG type to `global coarsening`.
 * <br>
 * (Francesco Radica, 2026/07/09)
 *
 * <li> Fixed: The periodicity constraints of phi-periodic spherical shells were
 * silently wrong for scalar fields such as the pressure of the matrix-free
 * Stokes solvers, because the rotation matrix was misinterpreted as a face
 * interpolation matrix. Forcing the GMG solver on periodic spherical shells
 * used to produce wrong velocities near the periodic boundary.
 * <br>
 * (Francesco Radica, 2026/07/09)
 *
 * <li> Fixed: When checkpointing a model with FastScape mesh deformation, the
 * basement and silt fraction fields were saved using the surface elevation
 * accessor (fastscape_copy_h_) rather than the basement and silt fraction
 * accessors. On restart this could corrupt the basement and silt fraction
 * state, making the computed sediment thickness (elevation - basement)
 * collapse to nearly zero. The basement and silt fraction accessors
 * (fastscape_copy_basement_ and fastscape_copy_f_) are now used.
 * <br>
 * (Z. Lei, 2026/07/08)
 *
 * <li> Changed: The viscoplastic strain invariant particle property now supports the compositing
 * material model when it uses the visco-plastic material model to compute viscosity.
 * The interface has also been extended to allow checking for the use of a specific
 * material model, even when that model is used as a subordinate component within
 * a compositing material model.
 * <br>
 * (Haoyuan Li, 2026/03/11)
 *
 * <li> Changed: The dynamic core statistics postprocessor no longer stores or
 * serializes a separate copy of the dynamic core state. It now reads the current
 * core data directly from the dynamic core boundary temperature plugin, making the
 * boundary temperature plugin the single owner of the restart state.
 * <br>
 * (Francesco Radica, 2026/06/04)
 *
 * <li> Changed: If specifying friction angles as a function, we now
 * limit the number of function expressions to the number of
 * chemical compositional fields plus one, instead of
 * specifying functions for all types of compositional fields.
 * <br>
 * (Arushi Saxena, 2026/06/03)
 *
 * <li>  Added: Different particle managers can
 *  now also have different postprocessing
 *  output intervals and formats and exclude
 *  different particle properties.
 *  <br>
 *  (Anne Glerum 2026/05/27)
 *
 * <li> New: ASPECT can now use the Tpetra packages provided by
 * Trilinos if configured with deal.II 9.8.0 or higher.
 * This means ASPECT can now use Trilinos 17 for its vector,
 * matrix, and solver types. While accuracy and compatibility
 * have been confirmed, performance may deteriorate
 * compared to the older Epetra packages, until solver
 * parameters have been investigated and tuned.
 * <br>
 * (Rene Gassmoeller, 2026/05/26)
 *
 * <li> Changed: Added a test documenting that the dynamic core boundary temperature
 * model can be combined with the standard function boundary temperature model to
 * prescribe a spatially variable outer boundary temperature. Set
 * `Dynamic core/Outer temperature = 0` when the function model should provide the
 * full outer boundary temperature, and make the function return zero on the CMB so
 * that the dynamic core model remains responsible for the evolving CMB
 * temperature.
 * <br>
 * (Francesco Radica, 2026/05/23)
 *
 * <li> Fixed: Particle managers are now initialized and serialized when particles are
 * used to advect compositional fields, even if the particles postprocessor is not
 * active.
 * <br>
 * (Francesco Radica, 2026/05/21)
 *
 * <li> New: The 2d spherical shell geometry model now supports
 * 180 degree phi-periodic shells, including particle advection across
 * the periodic boundary.
 * <br>
 * (Francesco Radica, 2026/05/21)
 *
 * <li> Fixed: grain size exponent was not divided by stress exponent.
 * Exponent on grain size is now grain_size_exponent / stress_exponent.
 * Also default values based on Goldsby & Kohlstedt (2001) was changed to match the unit in ASPECT code.
 * The unit of prefactor in Goldsby & Kohlstedt (2001) is MPa**-n * meter**m /s. The changed default is Pa**-n * meter**m / s.
 * The unit of activation energy in Goldsby & Kohlstedt (2001) is kJ/mol. The changed default is J/mol.
 * <br>
 * (Hyunseong Kim, 2026/05/20)
 *
 * <li> Added: a configurable linear solver failure strategy to allow iterative Stokes solver schemes to continue nonlinear iterations after linear solver failures, including checkpoint/restart support, and tests for both AMG and GMG Stokes solvers.
 * <br>
 * (Haoyuan Li, 2026/05/13)
 *
 * <li> Added: additional named outputs for diffusion and dislocation viscosities in the visco-plastic material model, including rheology-specific handling of inactive flow laws and visualization support through the additional material output system.
 * <br>
 * (Haoyuan Li, 2026/05/12)
 *
 * <li> Added: A new visualization postprocessor 'prescribed solution', which outputs
 * whether the solution components are prescribed by
 * the prescribed solution plugin system and if so to which value.
 * <br>
 * (Haoyuan Li, 2026/05/11)
 *
 * <li> Fixed: The deal.II memory-leak workaround in the StrainDependent and
 * Elasticity rheologies (used with deal.II older than 9.8.0-pre, see
 * dealii/dealii#19328) now also resets the cached velocity-gradient
 * FEPointEvaluation, which was previously missed and kept leaking.
 * <br>
 * (Ninghui Tian, 2026/05/06)
 *
 * <li>  Fixed: The order of silt and sand parameters has
 *  been switched in the FastScape mesh deformation
 *  plugin call that sets the marine parameters in
 *  FastScape. Now silt parameters are used for the
 *  finegrained fraction and sand parameters for the
 *  coarsegrained sediments in the marine domain.
 *  <br>
 *  (Anne Glerum 2026/05/05)
 *
 * <li> Fixed: Fields advected with the `darcy field` advection method
 * now correctly use the Darcy velocity instead of the solid velocity
 * to compute the necessary entropy viscosity stabilization.
 * <br>
 * (Daniel Douglas, 2026/04/21)
 *
 * <li> Changed: Added checkpointing parameters to configure how many checkpoint
 * slots to keep, restart from a specific checkpoint slot, or restart from
 * the checkpoint whose saved time is closest to a requested model time.
 * <br>
 * (Ninghui Tian, 2026/03/27)
 *
 * <li> Changed: Clarify the default filename for user termination in the manual
 * <br>
 * (Max Rudolph, 2026/03/18)
 *
 * <li> Added: a new instance of 'initial composition' in the prescribed solution plugin system.
 * This plugin prescribes compositional fields from the initial composition model within a
 * region defined by an indicator function.
 * <br>
 * (Haoyuan Li, 2026/03/11)
 *
 * <li> Added: a new instance of 'initial temperature' in the prescribed solution plugin system.
 * <br>
 * (Haoyuan Li, 2025/12/30)
 *
 * <li> Fixed: A memory leak in deal.II 9.7 and older for certain material models
 * has been fixed with a work-around.
 * <br>
 * (Timo Heister, 2026/03/06)
 *
 * <li> New: Initial topography for mesh deformation plugins
 * can now alternatively be provided by modifying constraints
 * directly.
 * <br>
 * (Timo Heister, 2026/02/04)
 *
 * <li> Improved: The mesh deformation solver now supports higher-order
 * polynomial degrees in GMG-based setups. A new parameter,
 * `Mesh deformation mapping order`, controls the polynomial degree
 * used by the mesh deformation mapping; the default value `auto`
 * uses the larger of 4 and the Stokes velocity polynomial degree
 * for curved geometries, and 1 otherwise.
 * <br>
 * (Ninghui Tian, Rene Gassmoeller, Timo Heister, 2026/03/30)
 *
 * <li> Fixed: Tangential velocity boundary conditions along
 * boundaries with mesh deformation did not
 * correctly include the mesh deformation when computing
 * the tangential direction, leading to flow into and out of
 * the domain. This is fixed now.
 *
 * <br>
 * (Bob Myhill, Rene Gassmoeller and Timo Heister, 2026/01/12)
 *
 * <li> Fixed: The prescribed "velocity function" feature to set internal velocities now works correctly with the GMG solver.
 * <br>
 * (Timo Heister, 2026/01/01)
 *
 * <li> Added: a new instance of 'temperature function' in the prescribed solution plugin system.
 * <br>
 * (Haoyuan Li, 2025/12/30)
 *
 * <li> Changed: The particle interpolator 'bilinear least squares'
 * was renamed to 'linear least squares' and all references in the
 * repository have been renamed as well. This is to indicate that
 * while the interpolator was implemented as a bi/trilinear function
 * in the past, it had been changed since ASPECT 2.4.0.
 * <br>
 * (Rene Gassmoeller, 2025/12/12)
 *
 * <li> ASPECT now requires deal.II 9.6 or newer and a C++ compiler with C++17 support.
 * <br>
 * (Timo Heister, 2025/12/10)
 *
 * <li> Changed: The location of files output by the topography,
 * sea_level, heat_flux_map, dynamic_topography, and geoid
 * postprocessors are now output to postprocessor specific
 * subdirectories.
 * <br>
 * (Daniel Douglas, 2025/11/04)
 *
 * <li> New: ASPECT's cookbooks and benchmarks are now marked
 * with a number of tags that simplify finding models which
 * make use of particular features. There is a new page
 * in the documentation called 'Page index' that lists
 * all existing tags and their corresponding documentation
 * pages.
 * <br>
 * (Rene Gassmoeller, 2025/09/16)
 *
 * <li> Added: Options to add particles to cells during model runs using
 * a point density function or histogram based method to add particles
 * to regions of cells which are lacking in particles.
 * The point of the new particle addition method is to reduce imbalances in
 * the particle distribution which can result from randomly selecting locations
 * in which to add particles to cells.
 * <br>
 * (Jarett Baker-Dunn, 2025/09/26)
 *
 * <li> Changed: Improved the maximum horizontal compressive stress
 * postprocessor by adding stress averaging to ensure accurate
 * angle calculations and support for elastic rheology.
 * <br>
 * (Ninghui Tian, Rene Gassmoeller, 2025/09/27)
 *
 * <li> Added: An option to delete excess particles from cells using a point
 * density function which selects the most clustered particles to delete first.
 * The point of the new particle removal method is to avoid some issues which
 * can result from randomly selecting excess particles to delete.
 * <br>
 * (Jarett Baker-Dunn, 2025/09/26)
 *
 * <li> Changed: ASPECT's radial gravity plugins used to live in the same
 * header and source files contrary to most of our other plugins, which
 * each have one file per plugin. These plugins have been split into
 * separate files. The old header file 'gravity_model/radial.h' is
 * now deprecated and will be removed in the future. In addition,
 * the 'radial earth-like' plugin, which was disabled a long time ago,
 * was finally removed from the code.
 * <br>
 * (Rene Gassmoeller, 2025/09/22)
 *
 * <li> Changed: Added a scaling factor to allow time steps to be specified
 * in years instead of seconds.
 * <br>
 * (Ninghui Tian, 2025/09/20)
 *
 * <li> Added: a new prescribed solution plugin system that allows adding regional
 * solutions as constraints to the model. Introduced a first example plugin,
 * 'velocity function', which prescribes the velocity based on user-defined functions.
 * <br>
 * (Haoyuan Li, 2025/09/20)
 *
 * <li> Fixed: ASPECT used several features that could fail in MPI
 * communication in case an assert was triggered. The model
 * would end up in a communication deadlock (a model that hangs
 * without output) instead of correctly crashing and producing
 * an error message. This was fixed.
 * <br>
 * (Rene Gassmoeller, 2025/09/16)
 *
 * <li> New: Added a new nonlinear solver scheme 'iterated Advection, no Stokes'
 * that iterates the advection equation and uses a prescribed Stokes
 * solution.
 * <br>
 * (Rene Gassmoeller, 2025/09/02)
 *
 * <li> Added: A cookbook for CPO-induced anisotropic viscosity.
 * <br>
 * (Yijun Wang, 2025/08/27)
 *
 * <li> Fixed: On Mac systems, ASPECT sometimes crashed at the very end of the
 * program run with an exception triggered by the CGAL library that was
 * caused by the settings for floating point exceptions. This is now
 * fixed.
 * <br>
 * (Timo Heister, Wolfgang Bangerth, 2025/08/04)
 *
 * <li> Changed: Deprecate the input parameter 'Use years in output instead of seconds'
 * and replace with a new parameter called 'Use years instead of seconds' to indicate
 * it is also used in input parameters.
 * <br>
 * (Rene Gassmoeller, 2025/08/04)
 *
 * <li> Added: A new postprocessor 'timing statistics' that writes the
 * information about wall time spent in different sections of the
 * computation into the statistics file.
 * <br>
 * (Rene Gassmoeller, 2025/07/29)
 *
 * <li> Added: ASPECT now has a new gravity model plugin called 'Radial with tidal potential' for tidal forces.
 * This plugin is useful for modeling long-term interior and surface evolution in moons orbiting a large planet.
 * <br>
 * (Hyunseong Kim, Antoniette Greta Grima, Wolfgang Bangerth 2025/07/16)
 *
 * <li> Added: A boolean parameter in the bingham
 * average plugin to choose rotation matrix or
 * Euler angle representation.
 * <br>
 * (Yijun Wang, 2025/07/09)
 *
 * <li> Changed: The entropy method can now be used without loading additional shared
 * libraries. If there is at least one compositional field of type 'entropy'
 * the entropy equation will govern the thermodynamic evolution of the model
 * and the temperature equation will only be used to compute heat conduction.
 * <br>
 * (Rene Gassmoeller, 2025/07/08)
 *
 * <li> Changed: The entropy method can now also be used with the 'compositing'
 * material model. Additionally, the adiabatic conditions plugin
 * 'compute entropy profile' now considers initial composition to compute
 * the initial reference profile, identical to the 'compute profile' plugin.
 * <br>
 * (Rene Gassmoeller, 2025/07/04)
 *
 * <li> Fixed: The default value for the parameter 'Lateral viscosity file name'
 * in the material model 'entropy model' referred to a non-existent file.
 * This parameter now refers to an example data file.
 * <br>
 * (Rene Gassmoeller, 2025/07/04)
 *
 * <li> Added: The option for the user to use
 * the adiabatic pressure instead of the full
 * pressure in the Tian 2019 reaction model.
 * <br>
 * (Daniel Douglas, 2026/06/30)
 *
 * <li> New: The geometric multigrid (GMG) solver now also supports running
 * without viscosity averaging. This has the same memory requirement and
 * computational cost per iteration as Q1 averaging. Note that we still
 * need to average the viscosity for the GMG hierarchy, which likely
 * results in higher iteration counts compared to computations with
 * viscosity averaging enabled.
 * <br>
 * (Timo Heister, 2025/06/29)
 *
 * <li> Added: Initial topography can now be imposed on the
 * 'box with lithosphere boundary indicators' geometry
 * model.
 * <br>
 * (Arushi Saxena, 2025/06/24)
 *
 * <li> Fixed: A thin 3D spherical shell with an opening angle of 90 degrees
 * (octant of a sphere) would crash, because of a bug of how boundary
 * indicators were assigned. This bug is fixed in deal.II 9.7 and the
 * fix was backported into ASPECT.
 * <br>
 * (Rene Gassmoeller, 2025/06/20)
 *
 * <li> Added: An entropy statistics postprocessor to the multicomponent entropy averaging material model.
 * This postprocessor calculates for every time step the average number of iterations to equilibrate the temperature between the multiple components
 * <br>
 * (Ranpeng Li 2025/06/19)
 *
 * <li> Added: ParticleDistributionStatistics postprocessor which calculates
 * some statistics about the clustering and distribution of particles within
 * cells. Statistics are computed from a point-density function of the particles
 * and include standard deviation and maximum and minimum of the point-density
 * function.
 * <br>
 * (Jarett Baker-Dunn, 2025/06/19)
 *
 * <li> Changed: The tomography_based_plate_motions cookbook now has
 * added functionality to import slabs from the Slab2 model and
 * an initial imposed topography.
 * <br>
 * (Arushi Saxena, Juliane Dannberg, Rene Gassmoeller, 2025/06/19)
 *
 * <li> Added: Implemented Olivine D-type fabric with
 * the major slip system {0kl}[001] using an updated
 * generalized CPO (Crystallographic Preferred Orientation) algorithm.
 * This enhancement now reproduced lab results.
 * <br>
 * (Xiaochuan Tian, 2025/06/18)
 *
 * <li> Added: there is now a new particle property that saves the local velocity gradient into a particle property.
 * <br>
 * (Agi Kiraly and Yijun Wang, 2025/06/18)
 *
 * <li> Added: Grain boundary sliding flow law for ice-1 from Goldsby & Kohlstedt, 2001
 * <br>
 * (Antoniette Greta Grima, 2025/06/18)
 *
 * <li> Changed: The deviator and second invariant of symmetric tensors are modified to
 * be consistent with the plane strain assumption in 2D. This fixes the
 * inconsistency between some material models and the Stokes assemblers when
 * computing the deviatoric strain rate. It also results in better nonlinear solver
 * convergence due to some fixes to the Newton solver.
 * <br>
 * (Yimin Jin, 2025/06/17)
 *
 * <li> New: Add support for specifying
 * bedrock_river_incision_rate (Kf) and bedrock_transport_coefficient (Kd)
 * in Fastscape Fortran using user-defined functions.
 * These parameters represent climate and rock erodibility conditions.
 * They can now vary in both space and time through 2D and time-dependent functions.
 *
 * (Liang Xue, Derek Neuharth, 2025/06/17)
 *
 * <li> Added: The Steinberger material model now includes Drucker-Prager plasticity,
 * using the same yield strength for all chemical compositions.
 * <br>
 * (Qianyi Lu, 2025/06/17)
 *
 * <li> Added: Generalized CPO algorithm for slip systems
 * that supports additional crystal structure types
 * with existing Olivine fabrics benchmarked.
 * The monoclinic mineral clinopyroxene is implemented
 * as an example.
 * <br>
 * (Xiaochuan Tian, 2025/06/17)
 *
 * <li> Changed: The prescribed dilation is now split into two parts: a
 * pressure-dependent part that is moved to the left-hand side of the mass
 * conservation equation, and a pressure-independent part that stays on the
 * right-hand side. In addition, the effect of dilation on the momentum
 * conservation equation for incompressible models is no longer represented by
 * the dilation term on the right-hand side, but by the volumetric strain rate
 * term on the left-hand side. These changes enable the Drucker-Prager model
 * to produce associated plastic flows, in which the shear bands develop along
 * the Coulomb angle.
 * <br>
 * (Yimin Jin, 2025/06/16)
 *
 * <li> Added: The ability for specifying a sea level curve in the
 * FastScape Fortran plugin, allowing it to be defined either
 * as a constant or as a user-defined time-dependent function.
 * <br>
 * (Liang Xue, John Naliboff, Rene  Gassmoeller, Wolfgang Bangerth, 2025/06/16)
 *
 * <li> Added: The entropy model now supports multiple components.
 * The method assumes there is a background field, whose percentage doesn't need to be specified.
 * The adiabatic reference profile is calculated based on the background entropy.
 * <br>
 * (Ranpeng Li, Juliane Dannberg, Rene Gassmoeller, Bob Myhill, 2025/06/16)
 *
 * <li> Fixed: The Crystal Preferred Orientation plugin was not serializing
 * the internal state
 * <br>
 * (Menno Fraters, Yijun Wang, Agi Kiraly, 2025/06/16)
 *
 * <li> Added: A cookbook of crystal preferred orientation calculation
 * for olivine based on Fraters and Billen 2021
 * <br>
 * (Xiaochuan Tian, 2025/06/16)
 *
 * <li> New: There is now a new reaction model and particle property
 * that allows it to model the formation of basaltic oceanic crust
 * and harzburgitic lithosphere as mantle material reaches the
 * surface.
 * <br>
 * (Juliane Dannberg, 2025/06/16)
 *
 * <li> Added: A cookbook which replicates some results
 * of the publication by Burchardt et al (2012)
 * about the sinking of anhydrite blocks within
 * a Newtonian salt diapir.
 * <br>
 * (Cedric Thieulot, 2025/06/15)
 *
 * <li> New: A new prefactor has been added to the Drucker-Prager rheology model
 * that makes it possible to reduce (or increase) the yield strength of specific
 * compositions given in the input file
 * <br>
 * (Antoniette G. Grima, 2025/06/14)
 *
 * <li> Fixed: The conservation of water mass when being partitioned
 * within the solid phase and the free fluid phase in the
 * reactive fluid transport model.
 * <br>
 * (Daniel Douglas, 2025/06/14)
 *
 * <li> Added: A cookbook that demonstrates how to govern
 * phase transformations using a non-equilibrium thermodynamic
 * formulation for calculating reaction rates (using operator
 * splitting). The simple subduction model shows a dynamic
 * metastable olivine layer forming near the olivine -->
 * wadsleyite phase transition.
 * <br>
 * (Buchanan Kerswell, 2025/06/14)
 *
 * <li> Added: A density anomaly visualization postprocessor that can output density
 * anomaly from the lateral average density or from the adiabatic density reference profile.
 * <br>
 * (Qianyi Lu, 2025/06/13)
 *
 * <li> New: A "current surface" postprocessor was added that stores the
 * current surface including mesh deformation. This quantity can be visualized using
 * the new 'depth with mesh deformation' postprocessor. This function is
 * currently only available for 2D box geometries.
 * (Derek Neuharth, 2025/06/13)
 *
 * <li> New: The multicomponent compressible equation of state
 * and material model now support the use of phase transitions.
 * <br>
 * (John Naliboff, 2025/06/13)
 *
 * <li> Added: ASPECT now has a new heating model plugin called 'tidal heating' for diurnal tides.
 * This plugin is useful for modeling long-term interior evolution in moons orbiting a large planet.
 * Distribution of tidal strain rate can be selected between 'constant' and 'latitudinal variation'.
 * <br>
 * (Hyunseong Kim, 2025/06/13)
 *
 * <li> Added: ParticleDistributionScore postprocessor which scores how clustered
 * particles are within cells on a scale from 0-1, with 1 being most clustered.
 * <br>
 * (Jarett Baker-Dunn, 2025/06/13)
 *
 * <li> New: There is now a new type of boundary condition for the
 * temperature equation: Robin boundary conditions, which allow
 * for prescribing a linear combination of a prescribed
 * temperature (Dirichlet) and a prescribed heat flux (Neumann).
 * In other words, we can now prescribe a heat flux that depends
 * on how much the current temperature differs from a prescribed
 * temperature.
 * <br>
 * (Juliane Dannberg, 2025/06/12)
 *
 * <li> Added: A post-processor which outputs the max
 * fluid velocity.
 * <br>
 * (Daniel Douglas, 2025/06/12)
 *
 * <li> Fixed: The velocity statistics postprocessor to use
 * the correct quadrature formula when determining the
 * maximum velocity within the model.
 * <br>
 * (Daniel Douglas, 2025/06/12)
 *
 * <li> Changed: The file layout of checkpointing files has been changed. Now,
 * files are stored in folders rotating between output/restart/01/,
 * output/restart/02/, and output/restart/03/. To reuse checkpoint files
 * written earlier, you will need to: 1. create output/restart/01/,
 * 2. move restart files output/restart.mesh to output/restart/01/mesh
 * (similarly for other files), 3. create a file
 * output/restart/last_good_checkpoint.txt with the content "1".
 * <br>
 * (Timo Heister, 2025/06/11)
 *
 * <li> Changed: The MaterialModel::MaterialModelInputs and
 * MaterialModel::MaterialModelOutputs classes had member functions
 * `get_additional_input()` and `get_additional_output()` functions that
 * simply returned C-style pointers. The use of such pointers leaves it
 * entirely unclear who now owns the object pointed to. As a consequence,
 * these functions have now been deprecated and replaced by functions
 * `get_additional_input_object()` and `get_additional_output_object()`
 * that return their results in the form of `std::shared_ptr` values that
 * make clear that the calling place receiving the pointer now shares
 * ownership of the input or output object with the place that the object
 * is requested from.
 *
 * The old functions have been retained for backward compatibility
 * purposes, but they are now deprecated.
 * <br>
 * (Wolfgang Bangerth, 2025/06/11)
 *
 * <li> Added: The ability to define the maximum yield stress as a
 * compositionally dependent variable.
 * <br>
 * (Daniel Douglas, 2025/06/10)
 *
 * <li> Fixed: The random number generator used to create and delete particles was
 * being redeclared for every cell, causing particle indexes to be deleted in
 * the same order in every cell. Declaring the random number generator as a
 * member variable of the particle manager fixed this problem.
 * <br>
 * (Jarett Baker-Dunn, Rene Gassmoeller, 2025/06/10)
 *
 * <li> Fixed: The evaluation of the minimum/maximum value in
 * melt_statistics and in pressure_statistics was inaccurately being
 * done at the Gauss quadrature points. Use a Gauss-Lobatto quadrature
 * instead to improve the accuracy when these values lie at the
 * boundaries of a cell (e.g. the surface).
 * <br>
 * (Daniel Douglas, 2025/06/12)
 *
 * <li> Changed: Added support for specifying initial topography in
 * spherical coordinates in the 'function' initial topography plugin.
 * This enables users to prescribe initial surface perturbations
 * naturally on spherical geometries. The new feature extends the
 * plugin's flexibility for global-scale mantle convection models.
 * <br>
 * (Ninghui Tian, Rene Gassmoeller, 2025/06/09)
 *
 * <li> Changed: Added support for prescribing either pressure-only or
 * full traction vector boundary conditions in AsciiDataBoundary.
 * Introduced the input parameter "Prescribe pressure instead of
 * full traction" (default: true) to maintain backward compatibility.
 * Traction components can now also be specified in spherical coordinates.
 * This improves flexibility for boundary traction specifications.
 * <br>
 * (Ninghui Tian, Rene Gassmoeller, 2025/05/01)
 *
 * <li> Fixed: The grain size material model would locally compute a
 * wrong grain size growth term if the grain size was reset by
 * a phase transition, and was reduced below the minimum grain
 * size (by the phase reset or otherwise) in the same time step.
 * This is fixed now.
 * <br>
 * (Rene Gassmoeller, 2025/04/17)
 *
 * <li> Fixed: The entropy material model did under some conditions not fill the additional
 * material model outputs correctly. This would not affect the model results, but could
 * lead to floating point exceptions in the postprocessor 'named additional outputs'.
 * This is fixed now.
 * <br>
 * (Rene Gassmoeller, 2025/03/04)
 *
 * <li> Fixed: The grain size material model expected the wrong
 * number of entries for some input parameters if multiple
 * compositional fields were active. This is fixed now.
 * Also when new particles were generated in a boundary cell
 * with Dirichlet boundary conditions for the compositional
 * fields, the grain size particle property
 * would be interpolated incorrectly. This is fixed as well.
 * <br>
 * (Rene Gassmoeller, Menno Fraters, 2025/03/27)
 *
 * <li> Fixed: The 'surface stress' postprocessor used the wrong sign for stresses,
 * leading to incorrect stress directions in models with elastic deformation.
 * Additionally, the 'stress second invariant' postprocessor had incorrect
 * magnitude values due to a sign convention issue.
 * <br>
 * (Ninghui Tian, Rene Gassmoeller, 2025/03/21)
 *
 * <li> New: There is now a new particle property called 'composition
 * reaction' that tracks the initial composition but also allows
 * for reactions between the different compositions (or reactions
 * just including one composition) at specific points in time
 * given in the input file.
 * <br>
 * (Juliane Dannberg, 2025/03/14)
 *
 * <li> Fixed: The statistics file contained an invalid number of iterations
 * for variables that were not solved (e.g. because compositions were
 * prescribed, or Stokes was only solved using cheap or only expensive
 * iterations). This is confusing, because it looks like an error and for the
 * Stokes solver it caused an incorrect iteration number to be reported. It is
 * more accurate to report 0 iterations if none were made, which is what the
 * postprocessor does now.
 * <br>
 * (Rene Gassmoeller, 2025/03/04)
 *
 * <li> Fixed: Adds current cell to material model inputs
 * in the cpo particle property to fix the Olivine: Karato
 * 2008 mineral type. Also added a test to test this feature.
 * <br>
 * (Menno Fraters and Daniel Douglas, 2025/02/25)
 *
 * <li> Fixed: If the nonlinear solver fails in a timestep with mesh
 * refinement and "cut timestep size" is selected as the
 * Nonlinear solver failure strategy, ASPECT now correctly
 * repeats the timestep first before refining the mesh. Before,
 * ASPECT would execute refinement/coarsening of the mesh first,
 * which could lead to changes of the mesh in each failure
 * cycle.
 * <br>
 * (Juliane Dannberg, 2025/02/14)
 *
 * <li> Added: A cookbook which demonstrates how to use the
 * extract_local_velocity.py script to take the output of
 * a global convection model and apply the velocity as
 * boundary conditions within a 3D regional spherical chunk.
 * <br>
 * (Daniel Douglas, 2025/02/11)
 *
 * <li> Fixed: The boundary temperature plugin 'dynamic core' used to
 * crash when the inner core was completely molten or completely
 * solid. This is fixed now.
 * <br>
 * (Francesco Radica, Rene Gassmoeller, 2025/02/06)
 *
 * <li> Changed: Separated out the tian2019 solubility
 * reaction model into it's own module independent
 * of the reactive fluid transport material model.
 * <br>
 * (Daniel Douglas, 2025/02/06)
 *
 * <li> Changed: ASPECT's boundary traction and boundary velocity manager classes are
 * now also derived from the common classes Plugins::ManagerBase. In order to
 * standardize the interface, the functions get_active_boundary_traction_conditions()
 * and get_active_boundary_traction_names() (and their velocity counterparts)
 * have been deprecated. They have been replaced by the new functions get_active_plugins(),
 * get_active_plugin_boundary_indicators(), get_prescribed_boundary_traction_indicators(),
 * and get_component_mask().
 * <br>
 * (Rene Gassmoeller, 2025/02/02)
 *
 * <li> Fixed: Constant modes are now used in the melt solver.
 * <br>
 * (Quang Hoang, Timo Heister, 2025/01/31)
 *
 * <li> Changed: ASPECT now re-generates particles in each initial
 * adaptive refinement cycle instead of only once after global
 * refinement. This means that particle locations during initial
 * adaptive refinement are chosen according to the generator
 * instead of randomly.
 * <br>
 * (Juliane Dannberg, 2025/01/29)
 *
 * <li> New: ASPECT now has an interface for plugins that describe thermal
 * conductivity. This is useful to share functionality to compute
 * thermal conductivity across material models. Material models can,
 * but do not have to make use of these plugins.
 * <br>
 * (Rene Gassmoeller, 2025/01/27)
 *
 * <li> Changed: Prescribed compositional fields are now copied from the material
 * model output in one combined operation (instead of being copied field by field).
 * This change improves the efficiency of the copy process.
 * <br>
 * (Rene Gassmoeller, 2025/02/28)
 *
 * <li> Changed: The entropy reader throws when the range of the provided look-up table
 * does not fully cover the entropy-pressure range in the model.
 * <br>
 * (Ranpeng Li, 2025/01/10)
 *
 * <li> Added: ASPECT can now use the 'include' keyword in parameter
 * files to include other parameter files, even if those other
 * parameter files include parameter files themselves. ASPECT now
 * also respects the parameters 'Dimension' and 'Additional shared
 * libraries' in the included parameter files, which was not the
 * case before.
 * <br>
 * (Rene Gassmoeller, 2024/12/04)
 *
 * <li> Changed: Renamed the variable 'PhaseFunctionInputs::phase_index' to
 * PhaseFunctionInputs::phase_transition_index'. The new variable name
 * is more precise since it is used to index phase transitions rather than phases.
 * Material models that make use of the old name will have to be adjusted.
 * <br>
 * (Haoyuan Li, 2024/12/01)
 *
 * <li> Added: there is now a new class of phase function that handles discrete phase transitions
 * by looking up the most dominant phases in a lookup table. This function can be used
 * to make the rheology of the visco-plastic material model dependent on the dominant
 * mineral phase.
 * <br>
 * (Haoyuan Li, 2024/11/07)
 *
 * <li> Changed: The function simulator::replace_outflow_boundary_ids to allow outflow boundary conditions
 * to be determined based on the velocity that each compositional field is advected with.
 * <br>
 * (Daniel Douglas, 2024/10/24)
 *
 * <li> New: We implement an option to use the Weighted BFBT preconditioner
 * introduced by Rudi et al (2017). This preconditioner was designed
 * for problems involving highly heterogeneous viscosities.
 * <br>
 * (Quang Hoang, Timo Heister, 2024/06/06)
 *
 * <li> Fixed: The latent heat material model now correctly reads in
 * the Viscosity prefactors that change the viscosities of
 * individual phases for each compositional field. To make this
 * work in a consistent way, the format of this input parameter
 * is now the same as for other phase transition inputs (it is
 * parsed as a map with keywords rather than a comma-separated
 * list), which is an incompatible change in the input file.
 * <br>
 * (Juliane Dannberg, 2024/01/24)
 *
 * <li>  Changed: The implementation of visco-elasticity and
 *  visco-elasto-plasticity has been updated to properly
 *  track stresses over time. This means that iterative
 *  advection schemes need to be used, as well as the DG
 *  method for compositions. For fields, operator splitting
 *  has to be switched on, while in the case of particles,
 *  the particle property 'elastic stress' takes care of
 *  the stress update. If a fixed elastic time step is
 *  chosen, two sets of stress fields have to be tracked
 *  one for the current time step, one for the previous time
 *  step.
 *  <br>
 *  (Anne Glerum, Robert Myhill, Rene Gassmoeller, Juliane Dannberg, John Naliboff, Gerry Puckett, Esther Heckenbach 2023/12/04)
 *
 * </ol>
 */
