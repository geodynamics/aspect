/**
 * @page changes_between_2.1.0_and_2.2.0 Changes between version 2.1.0 and version 2.2.0
 *
 * <p> This is the list of changes made after the release of ASPECT version
 * 2.1.0 for version 2.2.0. All entries are signed with the names of the author.
 * </p>
 *
 * <ol>
 *
 * <li> Changed: The Geodynamic World Builder has been updated to version 0.3.0.
 * <br>
 * (Menno Fraters, 2020/06/18)
 *
 * <li> New: Project to Q1 viscosity averaging implemented for the matrix-free
 * GMG Stokes solver. Using "project to Q1 only viscosity" averaging scheme
 * with the GMG solver recovers optimal convergence rates for velocity in
 * all benchmarking problems.
 * <br>
 * (Thomas C. Clevenger, 2020/06/16)
 *
 * <li> New: The IDR(s) Krylov method can now be used for solving the Stokes
 * system, only when using the matrix-free GMG solver. The parameter
 * s can be chosen by the user. This method has a short term recurrence
 * and can greatly reduce the memory requirement as compared to GMRES.
 * <br>
 * (Thomas C. Clevenger, 2020/05/25)
 *
 * <li> New: Modified benchmark of elastic stress build-up in a viscoelastic cantilever flexing under its own weight.
 *      The boundary conditions are modified to have an open top and bottom, which allows the ambient fluid to be more
 *      readily displaced while the cantilever flexes. Four python scripts are included which create various plots that are
 *      useful for comparison to an analytic solution from Geodynamics, Turcotte and Schubert. These 4 scripts plot:
 *      1) The stress profile near the attachment point of the beam
 *      2) The maximum/minimum stress with time
 *      3) The maximum deflection with time
 *      4) The shape of the beam when it is at its maximum flexed state.
 * <br>
 * (D. Douglas, G. Ito, J. Naliboff, 2020/05/13)
 *
 * <li> New: Added a benchmark that illustrates how to solve an entropy equation
 * instead of a temperature equation.
 * <br>
 * (Rene Gassmoeller, Juliane Dannberg, 2020/05/07)
 *
 * <li> New: Added flexibility to the parse_map_to_double_array input method. The new
 * input method will allow a single value for a given key (e.g. density of a
 * composition), even if multiple values (e.g. for different phases) are expected.
 * The single value will be expanded to be the same for all values associated with
 * that key. For examples of applying this in a input file, see the test
 * visco_plastic_phases.
 * <br>
 * (Haoyuan Li, 2020/05/06)
 *
 * <li> New: I added a Fibonacci spiral sampling scheme to the 'gravity point values'
 * postprocessor, which produces a uniformly distributed map. It can only be used
 * for global sampling. This sampling scheme is based on the Fibonacci suit and
 * generates a user defined number of sampling locations per radius rather than
 * the longitude-latitude grid used for the default non-uniform map sampling
 * scheme.
 * <br>
 * (Ludovic Jeanniot, 2020/04/07)
 *
 * <li> New: There is a new module in the MaterialModel::Rheology namespace
 * (ConstantViscosityPrefactors) for multiplying the viscosity of each
 * compositional field by a specific factor.
 * <br>
 * (John Naliboff, 2020/04/04)
 *
 *
 * <li> New: Benchmark of buoyancy-driven viscoelastic plate
 * flexure with a free surface from Choi et al., 2013,
 * JGR, v.118, p.2429-2444, doi:10.1002/jgrb.50148.
 * <br>
 * (J. Naliboff, C. Thieulot, 2020/04/03)
 *
 * <li> Remove: The Viscoelastic Plastic material model and test suite
 * has been removed. It has been superseded by the Visco Plastic
 * material model, which now has the option to simulate
 * viscoelastic-plastic deformation.
 * <br>
 * (John Naliboff, 2020/03/31)
 *
 *
 * <li> Added: the precision of gravity acceleration, potential and gradients written in the output
 * and statistics files of the "gravity point values" postprocessor may be changed by the user
 * via the parameter file under the name "set Precision in gravity output".
 * The default value is 12.
 * <br>
 * (Ludovic Jeanniot, 2020/03/30)
 *
 * <li> Changed: The named additional outputs of the visco plastic material model
 * now output the current internal angle of friction in degrees instead
 * of radians.
 * <br>
 * (Anne Glerum, 2020/03/19)
 *
 * <li> New: Average- minimum- and maximum gravity acceleration and potential are
 * written into the statistic file by the 'gravity point values' postprocessor.
 * <br>
 * (Ludovic Jeanniot, 2020/03/18)
 *
 * <li> New: The visco plastic material model now has an option to simulate
 * viscoelastic and viscoelastic-plastic deformation.
 * <br>
 * (John Naliboff, Dan Sandiford, 2020/03/12)
 *
 * <li> Changed: The sampling scheme names 'map' and 'list' of the 'gravity point values'
 * postprocessor were renamed more appropriately 'uniform distribution'
 * and 'list of points'. The previous names still work for backward compatibility.
 * Additionally, this postprocessor now computes the theoretical gravity potential
 * in addition to the existing outputs.
 * <br>
 * (Ludovic Jeanniot, 2020/03/10)
 *
 * <li> New: The geometric multigrid preconditioner for the Stokes system now uses
 * a non-flexible version of the  preconditioner for the cheap iterations,
 * allowing for the use of GMRES, which requires roughly half of the memory of
 * FGMRES.
 * <br>
 * (Thomas C. Clevenger, 2020/02/28)
 *
 * <li> Fixed: We used to output the entire statistics file at the end of each
 * time step and overwrite its previous content. But these files can be
 * large (in the tens of megabytes) and this could put a significant
 * strain on the file system for simulations that did a large number of
 * time steps with relatively small meshes (and consequently a short time
 * between time steps). The new scheme now implemented only writes
 * incremental updates to the statistics file.
 * <br>
 * (Wolfgang Bangerth, 2020/02/26)
 *
 * <li> Fixed: The hdf5 visualization output would write incorrect data in parallel if
 * output was written after resuming from a checkpoint, but before the next
 * adaptive mesh refinement was performed. This is fixed now.
 * <br>
 * (Max Rudolph, Timo Heister, 2020/02/25)
 *
 * <li> New: Any material model using the strain_dependent rheology
 * module now has the option to track plastic strain in a
 * compositional field where the initial plastic strain is removed.
 * To use this option, users must include a compositional field
 * with the name 'noninitial_plastic_strain'. This feature will be
 * useful in models where the initial plastic strain field is used
 * to localize deformation, but may then mask the initial stages of
 * fault growth.
 * <br>
 * (John Naliboff, 2020/02/25)
 *
 * <li> Fixed: The 'chunk' geometry model would produce very inaccurate normal vectors
 * at horizontal boundaries when a free surface boundary was prescribed. This lead
 * to pressure and strain-rate oscillations at tangential velocity boundaries and
 * was fixed in the underlying chunk geometry description.
 * <br>
 * (Rene Gassmoeller, Fiona Clerc, 2020/02/24)
 *
 * <li> Fixed: The 'dynamic topography' postprocessor did not respect the manifold when writing the dynamic topography into separate data files. The topography was computed correctly, but for models with curved manifolds the output location would be the arithmetic average of all surface vertices instead of the actual face midpoint on the spherical surface. This is fixed now.
 * <br>
 * (Rene Gassmoeller, Marie Kajan, 2020/02/14)
 *
 * <li> New: The geometric multigrid preconditioner for the Stokes system now supports
 * Picard iterations for the defect correction form of the Stokes equations.
 * <br>
 * (Thomas C. Clevenger, 2020/02/14)
 *
 * <li> New: Benchmark of deformation of free surface by traction boundary
 * conditions, over viscous/viscoelastic halfspace.
 * <br>
 * (F. Clerc, J. Naliboff, C. Thieulot, D. Douglas, G. Ito, 2020/02/05)
 *
 * <li> Fixed: There was a bug in the 'box with lithosphere' geometry model plugin
 * that lead to a crash when using periodic boundary conditions on any of the
 * additional lithospheric boundaries. This is fixed now.
 * <br>
 * (Anne Glerum, 2020/02/04)
 *
 * <li> Changed: The reading of the ascii data file in the "adiabatic boundary" plugin (initial temperature model) was initially using the mid-point approach. Now it has been changed to use the existing ascii-data functionality.
 * Existing parameter files and data files for the "adiabatic boundary" plugin will have to be checked and converted to the new format.
 * <br>
 * (Emmanuel Njinju and  Tahiry Rajaonarison,  2020/05/05)
 *
 *
 *
 * <li> Incompatibility: The option to use PETSc for linear algebra has been removed
 * until further notice.
 * <br>
 * (Timo Heister, 2020/01/24)
 *
 * <li> Fix: Fixed a bug in the particle generator probability density function when
 * using fewer particles than mpi processes.
 * <br>
 * (Menno Fraters and Rene Gassmoeller, 2020/01/24)
 *
 * <li> New: The 'read_and_distribute_file_content' function now checks if the requested file lookup is a URL.
 * If the user has the libdap libraries installed and is requesting data from a libdap server then the
 * data is pulled from the server, converted to the proper format, and passed in.
 * <br>
 * (Kodi Neumiller, 2020/01/16)
 *
 * <li> Added: Link to libdap if ASPECT_WITH_LIBDAP is specified at configure time.
 * <br>
 * (Kodi Neumiller and Timo Heister, 2020/01/16)
 *
 * <li> New: Implement the "no Advection, single Stokes" solver scheme.
 * <br>
 * (Timo Heister, 2020/01/07)
 *
 * <li> New: The chunk geometry model can now incorporate initial
 * topography from an ascii data file.
 * <br>
 * (Anne Glerum, 2019/12/06)
 *
 * <li> Fixed: The 'AsciiDataLookup' class had a bug that caused models to crash
 * if they loaded new non-equidistant data files after a previous file was
 * already loaded. This could only happen when using the 'ascii data' plugins for
 * boundary conditions with time-dependent non-equidistant data files. This is
 * fixed now.
 * <br>
 * (Rene Gassmoeller, Sibiao Liu, 2019/11/27)
 *
 * <li> New: The 'depth average' postprocessor now additionally computes the laterally
 * averaged density of vertical mass flux for each depth slice in the model.
 * <br>
 * (Rene Gassmoeller, 2019/11/08)
 *
 * <li> Added: The 'depth average' postprocessor now supports writing more than one
 * file in different output formats. The default has changed to write a
 * human-readable ascii text file in addition to the previous gnuplot file.
 * <br>
 * (Rene Gassmoeller, 2019/11/08)
 *
 * <li> Fixed: Resuming from a checkpoint of a box model with initial topography
 * was broken, because the grid transformation was not applied
 * again after restart. Now the box geometry model connects to the signal after
 * resuming the computation such that it applies the initial topography again.
 * <br>
 * (Anne Glerum, 2019/11/06)
 *
 * <li> Fixed: The ellipsoidal natural coordinate conversion functions now use/output
 * radius, lon, lat coordinates instead of radius, lat, lon.
 * <br>
 * (Anne Glerum, 2019/10/30)
 *
 * <li> New: The 'multicomponent incompressible' equation of state now allows multiple
 * different phases per compositional field. So far the 'visco plastic' material
 * model is the only material model that makes use of this new functionality.
 * <br>
 * (Rene Gassmoeller, 2019/09/29)
 *
 * <li> Changed: There is now a general class
 * `MaterialModel::Utilities::PhaseFunction` that can be used to model
 * phase transitions using a smooth phase function. The class allows one or
 * more transitions for each composition and the background material. As a
 * consequence the syntax for the 'latent heat' material model has changed if
 * more than one composition or more than one phase transition is used.
 * Please refer to the documentation of the `Phase transition depths` parameter
 * for the new syntax.
 * <br>
 * (Rene Gassmoeller, 2019/10/29)
 *
 * <li> Fixed: The volume of fluid composition advection algorithm the boundary fluxes
 * on inflow boundaries with fixed composition incorrectly. This has been fixed by
 * correcting the calculation for the amount of composition which is advected into
 * the modeled region during the assembly step.
 * <br>
 * (Jonathan Robey, 2019/10/19)
 *
 * <li> Changed: The viscoelastic plastic material has been reordered so
 * that the viscoelastic viscosity is calculated prior to yielding
 * and viscoelastic stresses are taken into account when comparing
 * the current stress against the plastic yield stress. This approach
 * is consistent with the approach outlined in Moresi et al., 2003,
 * Journal of Computational Physics, v. 184, p. 476-697.
 * <br>
 * (John Naliboff, Daniel Sandiford, 2019/10/13)
 *
 * <li> New: There is a new model in the MaterialModel::Rheology namespace for plastic
 * material behavior based on the Drucker Prager yield criterion. This functionality
 * replaces existing functionality in MaterialUtilities and is currently used in the
 * visco plastic, drucker prager, and dynamic friction material models.
 * <br>
 * (John Naliboff, 2019/10/10)
 *
 * <li> New: Geometric Multigrid now supports no normal flux boundaries for 2d and 3d spherical and shell domains.
 * <br>
 * (Thomas C. Clevenger and Timo Heister, 2019/10/05)
 *
 * <li> New: Added a new parameter 'Temperature method' that controls how the
 * temperature equation is solved. The parameter works analogous to the
 * 'Compositional field method' parameter, and the current options are 'field'
 * (the default temperature equation), and 'prescribed field', which allows the
 * material model to fill a PrescribedTemperatureOutputs object that will set the
 * temperature field to defined values. This parameter is mostly useful if you
 * have other ways to determine temperature than by solving the temperature
 * equation (e.g. deriving from a compositional field, or from particles).
 * <br>
 * (Rene Gassmoeller, 2019/09/25)
 *
 * <li> New: Added a class MaterialModel::Utilities::PhaseFunction that can handle
 * generic phase functions for phase transitions that are approximated by an
 * analytic function.
 * <br>
 * (John Naliboff, Rene Gassmoeller, 2019/08/22)
 *
 * <li> New: There are two new modules in the MaterialModel::Rheology namespace for
 * viscous creep material behavior, diffusion creep and dislocation creep.
 * The functionality in these modules was transferred from the visco plastic
 * material model. In the future other material models can use this functionality
 * without duplicating code.
 *
 * <li> Removed: Input parameters and functionality related to the old syntax
 * for strain softening in the ViscoPlastic material model and the
 * Rheology::StrainDependent namespace.
 * <br>
 * (John Naliboff, 2019/09/17)
 *
 *
 * <li> New: There is a new model in the MaterialModel::Rheology namespace for material behavior
 * that depends on finite strain. At present, the functionality in this
 * model deals strictly with strain softening behavior, which was
 * transferred from the visco plastic material model. However, in the future
 * other material models can use this functionality without duplicating code.
 *
 * <li> New: Added a statistics postprocessor 'load balance statistics' that computes
 * statistics about the number of owned cells, and particles if active, across all
 * MPI processes.
 * <br>
 * (Rene Gassmoeller, 2019/08/16)
 *
 * <li> New: The matrix-free, block GMG Stokes solver can now be used for compressible
 * flow computations (all except those involving implicit reference density).
 * <br>
 * (Conrad Clevenger, 2019/08/06)
 *
 * <li> New: There is now a solver scheme 'single Advection, iterated Newton Stokes'
 * that solves the temperature and composition equations once per timestep
 * and uses defect corrector Picard and Newton iterations to iterate out
 * the Stokes solution.
 * <br>
 * (Anne Glerum, 2019/07/10)
 *
 * <li> New: There is now a new approximation for the compressible convection
 * models that is called 'projected density field'. This approximation is
 * more accurate (but also more non-linear) than the existing approximation
 * ALA, and isentropic compression. The details of the approximation are
 * described in an article currently in review in Geophysical Journal
 * International (Gassmoeller, et al., 2019).
 * <br>
 * (Rene Gassmoeller, 2019/07/09)
 *
 * <li> New: There is a new material model plug-in that sets viscosity to
 * a specified constant above the lithosphere-asthenosphere boundary,
 * which is specified by an ascii file or a maximum lithosphere depth
 * value. Below this the viscosity is taken from any other available
 * material model. All other material properties are taken from the
 * base model.
 * <br>
 * (Sophie Coulson, 2019/06/28)
 *
 * <li> New: ASPECT now includes initial temperature and initial
 * composition plugins that use ASCII data files to define
 * the initial temperature or composition at a series of
 * layer boundaries. The first one or two columns define a
 * structured grid of horizontal points (x, y or lat, lon),
 * and the third gives the vertical position (z or r) of
 * the layer at that point. The remaining column(s) provide
 * the temperature/composition of the layer boundary at that
 * point. The plugin first linearly interpolates values
 * on the surface of the layer boundary, and then either
 * linearly interpolates between layer boundaries or uses the
 * value obtained on the overlying layer boundary.
 * <br>
 * (Bob Myhill, 2019/06/18) *
 * <li> New: There is now a new module structure for mesh deformation
 * plugins. The original free surface has become one of these plugins;
 * a function plugin is also available. The new structure enables the
 * user to prescribe movement of the mesh independent of the
 * Stokes velocity by setting mesh velocity boundary conditions
 * on one or more mesh deformation boundaries. These boundary
 * conditions are used to solve the Laplace problem for the velocity
 * of the internal mesh nodes, as was originally done for the free
 * surface alone.
 * <br>
 * (Rene Gassmoeller, Anne Glerum, Derek Neuharth, Marine Lasbleis 2019/06/18)
 *
 * Incompatibility: Due to the new mesh deformation module, there
 * is no backwards compatibility for input files that include a
 * free surface. The "Free surface" section has been moved into
 * the new "Mesh deformation" section. A free surface boundary can be
 * selected through "set Mesh deformation boundary indicators = top: free surface"
 * in the "Mesh deformation" section.
 * <br>
 * (Rene Gassmoeller, Anne Glerum, Derek Neuharth, Marine Lasbleis 2019/06/18)
 *
 * <li> Changed: The 'isothermal compression' approximation for compressible mantle
 * convection has been renamed to 'isentropic compression', which is a more
 * appropriate name. The approximation is only correct if one uses the isentropic
 * rather than the isothermal compressibility in material models.
 * <br>
 * (Rene Gassmoeller, 2019/06/14)
 *
 * <li> Fixed: There was a bug in the heating model manager that occurred when a
 * simulation used several heating models and some of them filled the rates
 * of temperature change heating model outputs, but others did not. Because
 * many heating models do not touch the rates of temperature change heating
 * model outputs at all, the values need to be reset to zero before those
 * heating models are evaluated, otherwise they will contain rates of
 * temperature change outputs from a heating model that was evaluated
 * previously (which would lead to unintended temperature changes). This is
 * fixed now by resetting the appropriate outputs to zero between the
 * evaluation of the different heating models.
 * <br>
 * (Juliane Dannberg, 2019/06/12)
 *
 * <li> Changed: The Geodynamic World Builder is not longer a git submodule, but the files are
 * now distributed with ASPECT.
 * <br>
 * (Menno Fraters, 2019/05/31)
 *
 * <li> New: ASPECT now includes a multicomponent compressible model,
 * where the parameters of each component and the bulk are computed
 * self-consistently. The pressure and temperature dependencies of the
 * material parameters are based on the Murnaghan equation of state,
 * with a thermal pressure term that is linear with temperature and an
 * isochoric heat capacity that is globally constant.
 * <br>
 * (Bob Myhill, 2019/05/30)
 *
 * <li> New: Extended spherical shell geometry model to
 * include custom mesh schemes.
 * <br>
 * (Ludovic Jeanniot, Marie Kajan and Wolfgang Bangerth, 2019/05/30)
 *
 * <li> New: Viscosity grooves benchmark (multiple low viscosity shear bands).
 * <br>
 * (Cedric Thieulot, 2019/05/29)
 *
 * <li> New: Advection in annulus benchmark modified from annulus benchmark.
 * Used in conjunction with SUPG to show that coarse meshes using EV creates
 * much higher diffusion.
 * <br>
 * (G. Euen, R. Gassmoeller, T. Heister, 2019/05/29)
 *
 * <li> New: Slab detachment benchmark first presented in Schmalholtz (2011) and
 * run with ASPECT by A. Glerum (2018).
 * <br>
 * (C. Thieulot and A. Glerum, 2019/05/28)
 *
 * <li> New: Fix a rare deadlock caused by adaptive refinement in parallel runs.
 * <br>
 * (Timo Heister, 2019/05/27)
 *
 * <li> New: ASPECT now includes a thermodynamically self-consistent
 * compressible material model that implements the
 * Modified Tait equation of state that is described in
 * Holland and Powell, 2011 (doi:10.1111/j.1525-1314.2010.00923.x).
 * <br>
 * (Bob Myhill, 2019/05/27)
 *
 * <li> New: There is a new matrix-free Stokes solver which uses geometric multigrid for
 * preconditioning in the A-block. Currently the new method can handle no slip BCs for velocity
 * on any mesh, and free slip boundary conditions for a box. Free surface, melt transport, and
 * compressible flow are not yet implemented.
 * <br>
 * (Conrad Clevenger, Timo Heister, 2019/05/27)
 *
 * <li> New: Polydiapirs benchmark from Weinberg and Schmeling 1992.
 * <br>
 * (Cedric Thieulot, 2019/05/26)
 *
 * <li> New: There is a new initial topography model that sets the
 * initial topography based on a function specified in the input
 * file.
 * <br>
 * (Anne Glerum, 2019/05/26)
 *
 * <li> New: Several benchmarks for advection problems got added to benchmarks/advection/<br>
 * (Timo Heister 2019/05/25)
 *
 * <li> New: Viscoplastic_strain_invariants particle plugin added to track
 * the plastic and/or viscous, or total strain with particles.
 * <br>
 * (Derek Neuharth and Anne Glerum, 2019/05/25)
 *
 * <li> New: The material models can now outsource the computation of the viscosity
 * into a separate rheology model. This can be used to remove duplicated code
 * between material models. As a consequence the default parameter for viscosity
 * in the 'Simpler' material model has been changed to 1e21 Pa s.
 * <br>
 * (Rene Gassmoeller, 2019/05/25)
 *
 * <li> New: There is a new initial temperature plug-in that sets a specified
 * lithosphere temperature above the lithsphere-asthenosphere boundary,
 * which is specified by an ascii file or a maximum lithosphere depth
 * value. Below this the initial temperature is set as NaN. This
 * plug-in can be combined with another using the 'replace if valid'
 * operator.
 * <br>
 * (Sophie Coulson, 2019/05/25)
 *
 * <li> New: There is a new material model that combines viscous, elastic
 * and plastic deformation for one or more materials. At present,
 * the model does not take into account diffusion or dislocation creep
 * flow laws and strain softening.
 * <br>
 * (John Naliboff, 2019/05/24)
 *
 * <li> New: There is a new initial temperature plugin that
 * sets up a steady-state continental geotherm for a 3-layer
 * lithosphere.
 * <br>
 * (Anne Glerum, 2019/05/24)
 *
 * <li> New: Added a "rigid shear" benchmark that is based on an analytical
 * solution of the Stokes equations, and can be used to measure the influence
 * of different advection methods (particles, fields) on the Stokes
 * solution. The benchmark is described in the submitted manuscript: Gassmoeller,
 * Lokavarapu, Bangerth, Puckett. "Evaluating the Accuracy of Hybrid Finite
 * Element/Particle-In-Cell Methods for Modeling Incompressible Stokes Flow".
 * <br>
 * (Rene Gassmoeller, 2019/05/24)
 *
 * <li> New: This adds a new, alternative stabilization method for the advection called SUPG. This is preliminary work.<br>
 * (Thomas C. Clevenger, Rene Gassmoeller, Timo Heister, Ryan Grove, 2019/05/24)
 *
 * <li> New: Relocate Lookup namespace and the associated material property readers from grain_size to material model utlilities.
 * <br>
 * (Paul Bremner, 2019/05/24)
 *
 * <li> Remove: The
 * SimulatorAccessor::simulator_is_initialized() function has been
 * removed. It was not used anywhere before. It has been superseded by
 * the SimulatorAccess::simulator_is_past_initialization() function.
 * <br>
 * (Wolfgang Bangerth, 2019/05/24)
 *
 * <li> New: There is now a function
 * SimulatorAccess::simulator_is_past_initialization() that can be used
 * to test whether the simulator has been completely initialized and has
 * started the time loop.
 * <br>
 * (Wolfgang Bangerth, 2019/05/22)
 *
 * <li> Changed: The
 * BoundaryVelocity::Manager::get_prescribed_boundary_velocity() function
 * has been removed. It had previously already been removed.
 * <br>
 * (Wolfgang Bangerth, 2019/05/22)
 *
 * <li> Changed: The
 * BoundaryVelocity::Manager::get_active_boundary_velocity_conditions()
 * function used to return a `std::vector` of `std::shared_ptr` objects
 * to heating models enabled for the current run. This has been changed:
 * It now returns a `std::vector` of `std::unique_ptr` objects.
 * <br>
 * (Wolfgang Bangerth, 2019/05/22)
 *
 * <li> Changed: The
 * InitialComposition::Manager::get_active_initial_composition_conditions()
 * function used to return a `std::vector` of `std::shared_ptr` objects
 * to heating models enabled for the current run. This has been changed:
 * It now returns a `std::vector` of `std::unique_ptr` objects.
 * <br>
 * (Wolfgang Bangerth, 2019/05/22)
 *
 * <li> Changed: The entropy viscosity method for stabilizing the advection equations
 * for temperature and composition was substantially improved. A bugfix allowed to
 * reduce the stabilization parameters to the old defaults that were published in
 * the Kronbichler et al 2012 article. Additionally, boundary layers that are
 * conduction dominated are not longer unnecessarily stabilized, which results in
 * a more accurate heat flux, a better overall heat-flux balance, and closer
 * agreement with existing benchmark cases. The old implementation was not wrong,
 * but it required an unnecessarily high resolution to correctly resolve boundary
 * layers. The 'artificial viscosity' postprocessor can always be used to check
 * the amount of applied artificial diffusion.
 * <br>
 * (Rene Gassmoeller, 2019/05/13)
 *
 * <li> Fixed: There was an interpolation bug in the 'gplates' boundary velocity plugin
 * that lead to constant spherical velocities in a small spherical wedge right
 * below the 0 meridian. The missing interpolation lead to a small velocity jump
 * right at the 0 meridian. This is fixed now.
 * <br>
 * (Rene Gassmoeller, 2019/05/08)
 *
 * <li> Changed: The lateral averaging functions, e.g. used in the 'depth average'
 * postprocessor were optimized for some geometry models (spherical shell, box,
 * two merged boxes). Their lateral averaging is now slightly less accurate (the
 * difference is typically on the order of 1e-3 at the coarsest resolutions and
 * disappears with higher refinement), but they will be around 5 (2D) to 25 (3D)
 * times faster.
 * <br>
 * (Rene Gassmoeller, 2019/05/05)
 *
 * <li> Changed: The
 * BoundaryTemperature::Manager::get_active_boundary_temperature_conditions()
 * function used to return a `std::vector` of `std::shared_ptr`
 * objects. This has been changed: The function now returns a
 * `std::vector` of `std::unique_ptr` objects.
 * <br>
 * (Wolfgang Bangerth, 2019/05/02)
 *
 * <li> Changed: The HeatingModels::Manager::get_heating_models() function
 * used to return a `std::vector` of `std::shared_ptr` objects to heating
 * models enabled for the current run. This has been changed: It now
 * returns a `std::vector` of `std::unique_ptr` objects.
 * <br>
 * (Wolfgang Bangerth, 2019/05/02)
 *
 * <li> Changed: The test project will no longer automatically re-configure,
 * when cmake is running. To force re-detection of new tests, you need
 * to run "make setup_tests" again.
 * <br>
 * (Timo Heister, 2019/05/26)
 *
 * <li> New: ASPECT now requires deal.II version 9.0.0 or newer.
 * <br>
 * (Timo Heister and Rene Gassmoeller, 2019/04/25)
 *
 * <li> New: There is a new termination criterion that cancels the model run
 * when a steady state average temperature is reached.
 * <br>
 * (Rene Gassmoeller, Juliane Dannberg, Eva Bredow, 2019/01/09)
 *
 * <li> Incompatibility: When selecting the 'dynamic core' plugin for the
 * boundary temperature, the plugin previously used the "Dynamic core"
 * spelling with a leading upper-case letter. However, all other plugins
 * use leading lower-case letters. Consequently, the plugin was changed
 * to follow the usual convention using all lower-case letters, even
 * though this is a change that requires changes to users' input files.
 * <br>
 * (Wolfgang Bangerth, 2018/06/21)
 *
 * </ol>
 */
