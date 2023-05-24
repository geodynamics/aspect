/**
 * @page changes_between_1.5.0_and_2.0.0 Changes between version 1.5.0 and version 2.0.0
 *
 * <p> This is the list of changes made after the release of ASPECT version
 * 1.5.0 for version 2.0.0. All entries are signed with the names of the author.
 * </p>
 *
 * <ol>
 *
 * <li> Fixed: The mass reaction term on the right-hand side of the
 * Stokes system in models with melt transport is now computed
 * correctly (instead of being zero) in models with operator
 * splitting.
 * <br>
 * (Juliane Dannberg, 2018/05/08)
 *
 * <li> Fixed: The nullspace removal had a minor bug that caused a small
 * remaining rotation in three dimensional models even if the net rotation
 * or angular momentum should have been removed. The bug was fixed.
 * <br>
 * (Rene Gassmoeller, Shangxin Liu, 2018/05/03)
 *
 * <li> Fixed: The material models <code>visco plastic</code> and <code>diffusion
 * dislocation</code> now use the reference density to calculate thermal
 * conductivity if the <code>Temperature equation</code> formulation is set
 * to <code>reference density profile</code>.
 * <br>
 * (John Naliboff, 2018/04/30)
 *
 * <li> Changed: The default value for AdditionalSeismicOutputs was changed from -1 to
 * a signaling Nan. This ensures that every plugin that creates seismic outputs
 * needs to use a material model that actually fills those outputs. Note that the
 * 'depth average' postprocessor now only outputs seismic velocity columns if they
 * are available. If they are not available, but explicitly requested, the
 * postprocessor will throw an error, and if 'all' quantities are requested the
 * seismic velocities will only be written if provided by the material model.
 * Previously, the columns were always written and filled with -1 if no velocities
 * were available.
 * <br>
 * (Rene Gassmoeller, 2018/04/27)
 *
 * <li> Removed: The 'solidus' initial temperature plugin was made somewhat redundant
 * by the 'adiabatic conditions/ascii_data' plugin, and since the inclusion of
 * melting the name 'solidus' is misleading. It was moved into the tests directory
 * until it is made more general, or used in a cookbook.
 * (Rene Gassmoeller, 2018/04/27)
 *
 * <li> New: Two new nonlinear solvers are now available: a defect correction Picard
 * iteration and a Newton iteration nonlinear solver, both available trough the
 * nonlinear solver called <code>Newton Stokes</code>.
 * <br>
 * The defect correction Picard (also called approximate Newton) iteration usually
 * follows the same convergence pattern as the normal Picard iteration (linear
 * convergence), but due to it's defect correction form can be more accurate and
 * can make use of the Eisenstat Walker (1994) method of determining the optimal
 * linear tolerance. It can be used by setting the <code>Max pre-Newton nonlinear
 * iterations</code> larger or equal to the <code>Max nonlinear iterations</code>
 * and the <code>Nonlinear Newton solver switch tolerance</code> to an unreachable
 * tolerance.
 * <br>
 * The Newton iteration is a defect correction Picard which also uses the
 * derivatives of the viscosity to the strain-rate and viscosity to the pressure
 * to converge faster to the correct solution. In contrast to the (defect
 * correction) Picard iteration, this iteration is not globally convergent. This
 * means that only when the current solution is close enough to the real solution,
 * the iteration converges quadratically. In many cases it is therefore advisable
 * to first perform a few defect correction Picard iteration, before turning on
 * the Newton solver.
 * <br>
 * Because the linear system resulting from the Newton iteration is not always
 * stable (Symmetric Positive Definite), we have options available to stabilize
 * it (see Fraters et al., in prep). This stabilization prevents the linear
 * iterative solver from crashing by forcing the system to be SPD, but may slow
 * down the convergence. We therefore implemented an option for a fail safe,
 * which automatically turns on the stabilization when the linear solver crashes.
 * <br>
 * (Menno Fraters and Wolfgang Bangerth, 2018/04/26)
 *
 * <li> Fixed: The melt simple model now correctly always uses the Thermal
 * bulk viscosity exponent parameter to compute the temperature
 * dependence of the compaction/bulk viscosity (instead of the Thermal
 * viscosity exponent that describes the temperature dependence of
 * the shear viscosity, which was used before in certain cases and
 * was not the correct parameter for the computation).
 * <br>
 * (Juliane Dannberg, 2018/04/25)
 *
 * <li> Changed: A number of Stokes solver parameters have been moved from the global
 * namespace into the subsection 'Solver parameters/Stokes solver parameters'.
 * The advection and composition solver tolerances have also been moved to the
 * 'Solver parameters' subsection. The update script was improved to handle all
 * changes automatically.
 * <br>
 * (John Naliboff, Rene Gassmoeller, 2018/04/23)
 *
 * <li> Changed: The update scripts for old parameter and source files have been moved
 * and renamed. The official way to update old files is now to call
 * doc/update_prm_files.sh or doc/update_source_files.sh with all file names that
 * should be updated.
 * <br>
 * (Rene Gassmoeller, 2018/04/20)
 *
 * <li> Changed: All of the nonlinear solver schemes have been renamed:
 * IMPES -> single Advection, single Stokes
 * Stokes only -> no Advection, iterated Stokes
 * iterated Stokes -> single Advection, iterated Stokes
 * Newton Stokes -> iterated Advection and Newton Stokes
 * iterated IMPES -> iterated Advection and Stokes
 * Advection only -> single Advection, no Stokes
 * The old names still work, but should be considered as deprecated.
 * <br>
 * (Juliane Dannberg, 2018/04/17)
 *
 * <li> Removed: The deprecated parameters 'Heating model/Model name', 'Model
 * settings/Include shear heating', 'Model settings/Include adiabatic heating',
 * and 'Model settings/Include latent heat' have been removed. Use 'Heating
 * model/List of model names' instead.
 * <br>
 * (Rene Gassmoeller, 2018/04/11)
 *
 * <li> Removed: The 'Model settings' subsection in the parameter file has been
 * removed. It contained a selection of parameters that conceptually belonged to
 * other subsections and were never moved to the correct places. The update script
 * in doc/update_prm_files_to_2.0.0.sed was updated to correctly identify and move
 * these parameters to their new subsections, and all files in the repository have
 * been updated using this script.
 * <br>
 * (Rene Gassmoeller, 2018/04/10)
 *
 * <li> Overhaul of the melt solver.
 * Changed: Models with melt transport now have a new preconditioner,
 * and the linear system is solved in a different way: Based on the
 * Darcy coefficient, each cell is classified as a melt cell
 * (where the melt transport equations are solved) or not a melt
 * cell (in this case, the compaction pressure dofs are constrained
 * to zero, and the equations that are solved are the Stokes system
 * without any melt-related terms). To achieve better convergence
 * for low melt fractions, the compaction pressure p_c is replaced
 * by a scaled compaction pressure p_c_bar, so that the linear
 * system we are solving is different than before, and the solution
 * vector now contains p_c_bar instead of p_c.
 * The advantages of these changes are much lower iteration counts in
 * models with low or zero porosity, and that we are now solving
 * the Stokes system if no melt is present.
 * <br>
 * (Juliane Dannberg, 2018/04/06)
 *
 * <li> Removed: The 'radial earth-like' gravity profile was not really earth-like, and
 * was in fact worse than the default values for most of the other available
 * gravity plugins. This plugin was removed and all occurrences replaced by the
 * 'ascii data' plugin, which by default uses a data file with values from the
 * PREM reference model.
 * <br>
 * (Rene Gassmoeller, 2018/04/04)
 *
 * <li> Changed: The 'radial linear' gravity model no longer goes towards a zero
 * gravity at the bottom of the domain, but allows setting the gravity at
 * the bottom to a fixed value.
 * The new default value was chosen to represent mantle gravity profiles.
 * If you used the old model, please restore the old behavior by setting the
 * parameter to 0.0.
 * <br>
 * (Rene Gassmoeller, 2018/04/04)
 *
 * <li> Removed: Due to low usage, and the lack of realistic applications, the code
 * paths for using anisotropic viscosity have been removed from the main assembly.
 * They are still available and tested as plugins for the anisotropic_viscosity
 * test and can be utilized if needed. This optimization speeds up the Stokes
 * assembly by up to 25%.
 * <br>
 * (Rene Gassmoeller, 2018/03/22)
 *
 * <li> Improved: The 'depth average' postprocessor was significantly optimized, and
 * now uses approximately ten times less computing time.
 * <br>
 * (Rene Gassmoeller, 2018/03/21)
 *
 * <li> Overhaul of the melt simple material model.
 * Fixed: The freezing of melt in the melt simple material model is
 * now more consistent: The freezing rate parameter determines the
 * freezing of all melt, and melt freezes if both the porosity and
 * the depletion are larger than the equilibrium melt fraction.
 * Changed: The melt simple material model can now only be used
 * with operator splitting.
 * <br>
 * (Juliane Dannberg, 2018/03/16)
 *
 * <li> New: We now also write a machine readable .json file that contains all
 * parameters of each model additionally to the .prm and .tex files. This
 * allows automated evaluation of model output (e.g. plotting statistics
 * over input parameters).
 * <br>
 * (Rene Gassmoeller, 2018/03/11)
 *
 * <li> Fixed: The melt material models now work correctly with operator
 * splitting and latent heat of melting. This is achieved by setting
 * a separate time scale for these reactions in the material model.
 * <br>
 * (Juliane Dannberg, 2018/02/21)
 *
 * <li> Changed: If the velocity is zero, the time step size is now set to
 * the maximum time step input parameter instead of arbitrarily being
 * set to 1 second.
 * <br>
 * (Juliane Dannberg, 2018/02/20)
 *
 * <li> New: There is now a new mesh refinement criterion that adds the
 * option to refine cells based on the compaction length, a
 * typical length scale of features in models with melt migration.
 * <br>
 * (Juliane Dannberg, 2018/02/12)
 *
 * <li> New: ASPECT now supports using the free surface in models with
 * melt migration, and with all nonlinear solver schemes.
 * <br>
 * (Juliane Dannberg, 2018/01/30)
 *
 * <li> New: The 'initial profile' adiabatic conditions plugin was
 * renamed to 'compute profile' as it can
 * now use a function to change the adiabatic surface pressure
 * and temperature over time. If this option is used the
 * adiabatic reference profile will be recomputed for every
 * timestep based on the new function values. If the new
 * option is not set the behavior is unchanged.
 * <br>
 * (Rene Gassmoeller, 2018/01/25)
 *
 * <li> New: There is now a gravity profile that follows the values of the Preliminary
 * Reference Earth Model (PREM), which is also the new default value for the
 * 'ascii data' gravity model.
 * <br>
 * (Rene Gassmoeller, 2018/01/19)
 *
 * <li> Changed: Removed min/max bounds on the specific heat capacity and
 * thermal expansivity inside the 'Steinberger' material model.
 * These bounds were only used if the 'Use latent heat' parameter
 * was set to true. To prevent silently changing the solution these
 * bounds were removed.
 * <br>
 * (Rene Gassmoeller, 2018/01/16)
 *
 * <li> Fixed: The 'heat flux map' visualization postprocessor used to compute
 * the heat flux through every boundary face and average the computed fluxes per
 * cell. At side boundaries of boxes the computed flux usually makes no sense,
 * and one is usually only interested in the surface or bottom heat flux. Thus
 * the plugin was changed to only compute fluxes at top and bottom boundaries.
 * <br>
 * (Rene Gassmoeller, 2018/01/15)
 *
 * <li> New: We now reuse the memory for matrices for compositional fields instead of
 * allocating one for each field. This will result in memory savings if a large
 * number of them are used.
 * <br>
 * (Timo Heister, 2018/01/12)
 *
 * <li> Changed: The boundary velocity plugins are now controlled by
 * a manager class that allows assigning several plugins at the same boundary.
 * This behavior is identical to the boundary temperature and boundary
 * composition plugins with the additional possibiltiy of assigning different
 * plugins to different boundaries, and only prescribing certain
 * components of the velocity.
 * <br>
 * (Rene Gassmoeller, 2018/01/05)
 *
 * <li> Changed: The boundary composition plugins are now controlled by
 * a manager class that allows assigning several plugins at the same
 * time. This behavior is identical to the boundary temperature plugins.
 * <br>
 * (Rene Gassmoeller, 2018/01/04)
 *
 * <li> New: A new material model 'damage rheology' is added.
 * This model uses a compositional field to store and evolve a
 * quantity that represents the grain size, which influences
 * the viscosity. All other material properties can be read from
 * data tables generated by Perplex or Hefesto, or set in a similar
 * way to the 'simple' material model. Details about the model
 * are available in Dannberg et al., 2017 in G-Cubed.
 * <br>
 * (Rene Gassmoeller, Juliane Dannberg, Bob Myhill, 2018/01/02)
 *
 * <li> Fixed: Memory consumption for matrix storage got reduced significantly by not allocating unused entries in the preconditioner matrix.
 * <br>
 * (Timo Heister, 2017/12/22)
 *
 * <li> Changed: The ASPECT header (the first lines of each output and the
 * beginning of log.txt files) now contains more information about the
 * version of the software and the used dependencies. This will
 * improve the reproducibility of old model results.
 * <br>
 * (Rene Gassmoeller, 2017/12/12)
 *
 * <li> New: Added a new "static" field type for compositional fields that are not
 * subject to advection/diffusion/reactions.
 * <br>
 * (Timo Heister, 2017/11/23)
 *
 * <li> New: Implement new mass conservation formulation "hydrostatic compression".
 * <br>
 * (Juliane Dannberg, Rene Gassmoeller, Timo Heister, 2017/11/22)
 *
 * <li> New: The option to prescribe a Stokes solution when using the
 * advection only solver scheme now also works in models where melt
 * transport is used.
 * <br>
 * (Juliane Dannberg, 2017/11/21)
 *
 * <li> Changed: The 'Morency and Doin' material model was removed from the available
 * material models, and moved to a cookbook folder. This material model is
 * really only useful in that particular setup.
 * <br>
 * (Rene Gassmoeller, 2017/11/15)
 *
 * <li> Changed: The screen output of the different solver schemes now has
 * the same style. All lines that output the total nonlinear residual
 * of the system iterated over begin with 'Nonlinear residual' and
 * then specify the system being solved.
 * <br>
 * (Juliane Dannberg, 2017/11/01)
 *
 * <li> New: There is a new cmake option ASPECT_COMPARE_TEST_RESULTS,
 * which allows to run the test suite without comparing its results.
 * This is useful for other machines than the reference tester, and to
 * check if one of the tests crashes. The option defaults to ON, which
 * preserves the old behavior of comparing the results.
 * <br>
 * (Rene Gassmoeller, 2017/10/31)
 *
 * <li> Fixed: A forgotten factor of 2 in assembling the preconditioner for
 * the Stokes equation, and a forgotten term for compressible models,
 * led to a suboptimal performance of the solver. The
 * answers were always correct, but it took more iterations than
 * necessary. This is now fixed, and results in significantly fewer
 * GMRES iterations when solving the Stokes equations.
 * <br>
 * (Menno Fraters, Wolfgang Bangerth, Rene Gassmoeller, Timo Heister, 2017/10/27)
 *
 * <li> Changed: The 'depth' function of geometry models now uses any existing initial topography instead of the undeformed geometry as reference level to compute the depth.
 * <br>
 * (Menno Fraters 2017/05/16)
 *
 * <li> Added: A new particle interpolation scheme called
 * harmonic average. If a cell contains zero particles,
 * then the harmonic average of the particle properties
 * in the neighboring cells is returned.
 * <br>
 * (Harsha Lokavarapu, 2017/10/21)
 *
 * <li> Changed: The hdf5 visualization output no longer automatically merges
 * vertices of adjacent cells to save disk space. This has lead to confusion
 * when investigating discontinuous output properties, and is now an input
 * parameter that needs to be explicitly enabled.
 * <br>
 * (Rene Gassmoeller, 2017/10/19)
 *
 * <li> New: The new 'function' boundary composition plugin allows
 * to prescribe compositional field boundary conditions using
 * an analytical function. It supports the
 * coordinate systems cartesian, spherical, and depth.
 * <br>
 * (Matt Weller, Rene Gassmoeller, 2017/10/01)
 *
 *
 * <li> New: The core of the Newton solver has been added to ASPECT and is functional.
 * This includes several tests to test that the Newton solver is functioning properly.
 * <br>
 * (Menno Fraters 2017/05/16)
 *
 * <li> New: Much of the non ASPECT specific particle functionality has been moved into
 * the new ParticleHandler class, that will eventually be transferred into
 * deal.II. Some particle related interfaces and plugins need to be adjusted, in
 * particular the Particle::Integrator plugins.
 * <br>
 * (Rene Gassmoeller, 2017/09/20)
 *
 * <li> New: Traction boundary conditions now also work in models
 * with melt transport.
 * <br>
 * (Juliane Dannberg, 2017/09/13)
 *
 * <li> New: a new boundary temperature plugin that the core mantle boundary temperature evolves through time
 * by calculating the heat flux through core mantle boundary and solving the core energy balance.
 * A postprocessor of the core evolution is also added. It is used in our paper Zhang & O'Neill [2016]
 * and it could be useful for people interested in long term planetary evolution.
 * <br>
 * (Siqi Zhang & Craig O'Neill, 2017/08/22)
 *
 * <li> Fixed: Visualization postprocessors now get the current cell, and output
 * material properties that the material model computes using information
 * related to the current cell correctly.
 * <br>
 * (Juliane Dannberg, 2017/08/21)
 *
 * <li> New: There is now a visualization postprocessor that outputs
 * the volumetric strain rate (the divergence of the velocity).
 * <br>
 * (Juliane Dannberg, 2017/08/17)
 *
 * <li> New: There is now a new cookbook that describes how to use melt
 * transport in a convection model. This change also includes a new
 * postprocessor that visualizes the melt fraction material property
 * that is computed by some material models.
 * <br>
 * (Juliane Dannberg, 2017/08/14)
 *
 * <li> New: Sparsity patterns of the system matrix now use the current_constraints
 * instead of the static constraints. This allows new kind of constraints in the
 * system and reduces the number of nonzero entries in the system matrix.
 * <br>
 * (Timo Heister, 2017/08/01)
 *
 * <li> New: With a recent developer version of deal.II (9.0.0.pre) ASPECT can now
 * output particle data using deal.II functionality. This allows compressed vtu
 * output, and the same MPI I/O grouping options that are already available for
 * the visualization postprocessor. The ascii particle output format no longer
 * exists and instead generates (nearly identically formatted) gnuplot output.
 * <br>
 * (Rene Gassmoeller, 2017/07/26)
 *
 * <li> Fixed: The diffusion dislocation and drucker prager material models
 * can now better handle cases where parameters or material model inputs
 * are slightly outside the range of what is expected for the mantle,
 * and in these cases they now compute reasonable viscosities instead
 * of triggering floating point exceptions.
 * <br>
 * (Juliane Dannberg, 2017/07/18)
 *
 * <li> Changed: Rename the name of dynamic topography subsection from "Dynamic Topography" to "Dynamic topography".
 * <br>
 * (Shangxin Liu, 2017/07/14)
 *
 * <li> Changed: In models with melt transport, the stabilization now
 * uses the maximum of solid and melt velocity, instead of using
 * the melt velocity alone.
 * <br>
 * (Juliane Dannberg, 2017/07/12)
 *
 * <li> New: ASPECT now supports operator splitting as a new solver scheme.
 * This allows it to decouple advection and reactions between
 * compositional fields, and can also be used for temperature changes
 * related to these reactions. As part of this scheme, material models
 * now have an optional additional output that computes reaction rates,
 * and heating models have a new output that compute the corresponding
 * rate of change in temperature. An example is given in a new cookbook.
 * <br>
 * (Juliane Dannberg, 2017/07/04)
 *
 * <li> New: There is now an 'install' target, so ASPECT can be installed using
 * "make install". The installation path is given by CMAKE_INSTALL_PREFIX.
 * <br>
 * (Timo Heister, 2017/06/18)
 *
 * <li> New: Added a cookbook that uses the dynamic topography and geoid postprocessor in a simple
 * spherical shell harmonic perturbation model run. Added the geoid output and lower boundary
 * dynamic topography output to the existing 'S20RTS initial condition' cookbook.
 * <br>
 * (Jacky Austermann, 2017/06/12)
 *
 * <li> New: Adds a Compositing material model to allow selecting material properties
 * as generated from a number of other material models.
 * <br>
 * (Jonathan Robey, 2017/05/26)
 *
 * <li> New: Boundary temperatures can now be prescribed as a combination of several
 * plugins in the same way it is possible for initial temperatures. Operators
 * to combine the models are also available.
 * <br>
 * (Rene Gassmoeller, 2017/05/24)
 *
 * <li> New: The viscoplastic material model now outputs the cohesion and friction angles
 * which can depend on the strain.
 * <br>
 * (Anne Glerum, 2017/05/23)
 *
 * <li> New: There is a new initial composition plugin called 'porosity' that initializes
 * the porosity to the equilibrium melt fraction if the material model provides one.
 * Other fields are left unchanged, and the plugin does not work if there is no
 * equilibrium melt fraction provided by the material model.
 * <br>
 * (Rene Gassmoeller, 2017/05/19)
 *
 * <li> New: The adiabatic initial conditions plugin 'initial profile' can now use a
 * reference function to initialize the composition, instead of using the initial
 * composition plugin. This is useful, if the initial composition plugin depends
 * on the adiabatic conditions and can therefore not be used for computing
 * the initial profile.
 * <br>
 * (Rene Gassmoeller, 2017/05/18)
 *
 * <li> New: Boundary names for all geometries are now consistently named "top" and "bottom". Old names are accepted for now.
 * <br>
 * (Timo Heister, 2017/05/17)
 *
 * <li> New: ASPECT can now initialize the compositional field by supplying
 * a list of plugins and operators determining how each plugin modifies the
 * field. The operators are specified as a new input file parameter
 * 'List of model operators', which takes a comma-separated list
 * taken from the selection add, subtract, minimum and maximum. This list
 * defaults to add if it is not given in the parameter file. If it is given,
 * it may either be of length 1 (in which case all plugins modify the
 * compositional field in the same way), or of the same length as
 * 'List of model names'.
 * <br>
 * (Bob Myhill, 2017/05/16)
 *
 * <li> New: ASPECT can now create initial temperature conditions by supplying
 * a list of plugins and operators determining how each plugin modifies the
 * temperature field. The operators are specified as a new input file
 * parameter 'List of model operators', which takes a comma-separated list
 * taken from the selection add, subtract, minimum and maximum. This list
 * defaults to add if it is not given in the parameter file. If it is given,
 * it may either be of length 1 (in which case all plugins modify the
 * temperature field in the same way), or of the same length as
 * 'List of model names'.
 * <br>
 * (Bob Myhill, 2017/05/16)
 *
 * <li> New: The sphere geometry model can now be used with the
 * spherical constant boundary composition and ellipsoidal
 * chunk can now be used with the spherical shell initial
 * temperature models.
 * <br>
 * (Bob Myhill, 2017/05/16)
 *
 * <li> New: There is now a visualization postprocessor for geoid data.
 * <br>
 * (Ian Rose, 2017/05/16)
 *
 * <li> New: Added a postprocessor which visualizes the spd factor for the newton solver. This factor is used in
 * the Newton solver to make sure the Jacobian stays Positive definite. For more info see the Utilities
 * function compute_spd_factor.
 * <br>
 * (Menno Fraters 2017/05/16)
 *
 * <li> New: Rewrote the Drucker Prager material model into the new evaluate style and added the derivatives of the
 * viscosity wrt the strain-rate and pressure. Also created a test to test the derivatives against a finite
 * difference derivative.
 * <br>
 * (Menno Fraters 2017/05/16)
 *
 * <li> New: Added a simple nonlinear test, which contains a simple nonlinear material model
 * which has the derivatives of the viscosity wrt the strain-rate and pressure and
 * tests the derivatives against a finite difference derivative. This will be useful later
 * on for benchmarks.
 * <br>
 * (Menno Fraters 2017/05/16)
 *
 * <li> New: A new particle property is added which indicates the presence of melt
 * greater than the melt transport threshold at the particles position. If melt
 * is not present a 0 is recorded. If melt is present a 1 is recorded. Only works
 * if there is a compositonal field named "porosity". A second new particle property
 * is also added which indicates the value for each compositional field. This can be
 * used as a preliminary way to track the petrological evolution of material.
 * <br>
 * (Joe Schools, Rene Gassmoeller, 2017/05/15)
 *
 * <li> New: The dynamic topography postprocessor now uses the consistent
 * boundary flux method for computing surface stresses, which is
 * significantly more accurate. The postprocessor also exposes a method
 * for getting the dynamic topography vector, so that the visualization
 * and geoid postprocessors do not need to duplicate effort.
 * <br>
 * (Ian Rose, 2017/05/15)
 *
 * <li> Changed: The default number of grouped files has been changed from 0 (i.e. the nunmber of
 * processors) to 16 in order to keep the number of output files small if not explicitly
 * stated otherwise.
 * <br>
 * (Jacky Austermann, 2017/05/15)
 *
 * <li> Changed: The boussinesq approximation formulation is now
 * renamed to Boussinesq approximation.
 * <br>
 * (Bob Myhill, 2017/05/14)
 *
 * <li> New particle interpolation scheme bilinear least squares
 * using singular value decomposition algorithm. Currently
 * only 2D models are supported. We chose a simple overshoot
 * and undershoot correction scheme, mainly to truncate based
 * on local cell maximum and minimum property value.
 * <br>
 * (Elbridge Gerry Puckett, Ying He, Harsha Lokavarapu 2017/05/13) *
 * <li> Fixed: Non-zero prescribed velocity boundary conditions used a linear instead
 * of the correct mapping. This caused incorrect values on curved boundaries.
 * <br>
 * (Timo Heister, 2017/05/12)
 *
 * <li> Changed: ASPECT's seismic anomalies post-processor produced spurious
 * anomalies as a result of laterally averaging velocities in depth slices,
 * and then smoothing slices. The new post-processor provides two options,
 * given by the Average velocity scheme variable. The option
 * reference profile provides a percentage anomaly relative to the
 * Vp and Vs calculated by evaluating the material model at the P-T
 * conditions given by adiabatic conditions. The lateral average provides
 * an anomaly relative to the laterally averaged velocity in the model at
 * that depth. This velocity is calculated by linear interpolation between
 * volumetric averages calculated within depth slices. The number of depth
 * slices in the domain is user-defined. This second option may give
 * odd results if there is a jump in seismic velocity not aligned with
 * the depth slices, so in general, the first option is preferred.
 * <br>
 * (Bob Myhill, 2017/05/12)
 *
 * <li> New: Functions [boundary - temperature; velocity; traction] now allow user defined Coordinate Systems [cartesian; spherical; depth].
 * <br>
 * (Matt Weller, 2017/05/12)
 *
 * <li> Changed: The names of ASPECT's data directories now follow the same naming
 * scheme used for the source and include directories. The conversion scripts
 * in the doc/ directory have been adjusted to correctly modify parameter and
 * source files to the new naming scheme.
 * <br>
 * (Rene Gassmoeller, 2017/05/12)
 *
 * <li> New: ASPECT can now compute geoid in 3D spherical shell geometry. The geoid is
 * calculated in spherical harmonic domain and the final results are transferred
 * into spatial domain. The users can also choose to output the spherical harmonic
 * coefficients of geoid, density contribution part, surface dynamic topography
 * contribution part, and CMB dynamic topography contribution part respectively.
 * <br>
 * (Shangxin Liu, Ian Rose, 2017/05/11)
 *
 * <li> New: Added functions to compute weighted p norm averages (Utilities::weighted_p_norm_average)
 * and the derivatives of these weighted p norm averages (Utilities::derivative_of_weighted_p_norm_average)
 * to namespace Utilities. The derivatives function is templated to be able to compute derivatives which
 * are doubles and tensors.
 * <br>
 * (Menno Fraters 2017/05/10)
 *
 * <li> Refactored the SolKz benchmark and added the SolKz benchmark with
 * compositional fields and the SolKz benchmark using active
 * particles. For all three cases, we used Q2_Q1 elements with
 * Q2 temperature and Q2 compositional fields for SolKz and
 * discontinuous piecewise constant composition field for the
 * active particle case.
 * <br>
 * (Elbridge Gerry Puckett, Ying He, Harsha Lokavarapu 2017/05/10) *
 * <li> Changed: The chunk geometry now removes the manifold id only from the boundaries and sets the correct
 * boundary objects for deal.II versions prior to version 9. For deal.II 9, only the manifold is used,
 * for which the push_forward_gradient function was implemented to compute the normal vectors. This way,
 * velocity anomalies along mesh refinement levels are removed.
 * <br>
 * (Anne Glerum, 2017/05/10)
 *
 * <li> New: Allow additional RHS force terms in the Stokes system by enabling the
 * parameter "Enable additional Stokes RHS" and filling
 * AdditionalMaterialOutputsStokesRHS in the material model.
 * <br>
 * (Timo Heister, 2017/05/09)
 *
 * <li> New: 'nearest neighbor' particle interpolator: gets properties from
 * the nearest particle in the same cell. If the cell is empty it gets
 * properties from the nearest particle in the nearest cell.
 * <br>
 * (Jonathan Perry-Houts, 2017/05/09)
 *
 * <li> Changed: The material model interface now contains a base class for
 * named additional outputs and a derived class that that can be filled
 * with output for the seismic velocities in the evaluate() function.
 * This replaces the member functions seismic_Vs and seismic_Vp in the
 * material model interface, and can now be used to generate graphical
 * output for other quantities as well. This change also removes the
 * postprocessors 'seismic_vp', and 'seismic_vs'. Instead there is a
 * postprocessor 'named additional outputs' that outputs all available
 * named output quantities. The file doc/update_prm_files_to_2.0.0.sed
 * was updated to modify parameter files to use the new postprocessor.
 * <br>
 * (Juliane Dannberg, Timo Heister, 2017/05/09)
 *
 * <li> Added: A new heating model (compositional_heating) that
 * allows users to specify a different internal heating rate
 * for each compositional field. Unlike the radioactive decay
 * heating model, the internal heating rate for each
 * compositional field is constant through time.
 * <br>
 * (John Naliboff, 2017/05/08)
 *
 * <li> New: The dynamic topography post processor now also calculates the
 * topography at the bottom of the domain (and not only the upper surface).
 * <br>
 * (Jacky Austermann, 2017/05/08)
 *
 * <li> Changed: The unit tests are now a separate cmake project that is
 * configured automatically. To speed up scanning of dependencies,
 * running ctest will no longer rebuild aspect if the binary is out of
 * date. Issues with using "ninja" as the generator have been fixed.
 * <br>
 * (Timo Heister, 2017/04/29)
 *
 * <li> Changed: ASPECT's 'particles' were initially introduced as 'tracers',
 * and the naming scheme was never unified. In order to
 * use a more consistent and simpler to remember structure all remaining
 * references to 'tracers' have been replaced by 'particles'.
 * This change is unfortunately incompatible with existing user plugins,
 * and existing parameter files. The provided scripts
 * doc/update_source_files_to_2.0.0.sed and doc/update_prm_files_to_2.0.0.sed
 * allow an easy conversion of existing user plugins and parameter files to the
 * new naming scheme.
 * <br>
 * (Rene Gassmoeller, 2017/04/24)
 *
 * <li> Changed: The names for all initial condition
 * plugins are now unified to the universal scheme 'initial_property' that
 * is planned for ASPECT 2.0.
 * This change is unfortunately incompatible with existing user plugins,
 * and existing parameter files. The provided scripts
 * doc/update_source_files_to_2.0.0.sed and doc/update_prm_files_to_2.0.0.sed
 * allow an easy conversion of existing plugins and parameter files to the
 * new naming scheme.
 * <br>
 * (Rene Gassmoeller, 2017/04/24)
 *
 * <li> Changed: The names for all boundary condition
 * plugins are now unified to the universal scheme 'boundary_property' that
 * is planned for ASPECT 2.0.
 * This change is unfortunately incompatible with existing user plugins,
 * but the provided script doc/update_source_files_to_2.0.0.sed allows
 * an easy conversion of existing plugins to the new naming scheme.
 * Parameter files are not affected, because the naming of
 * subsections already follows the new structure.
 * <br>
 * (Rene Gassmoeller, 2017/04/24)
 *
 * <li> New: Initial composition conditions can now be specified as a list
 * of plugins. Their initial compositions are summed to generate the
 * final initial condition. For this purpose there is a new input file
 * parameter 'List of model names', and a new manager class that owns all
 * of the initial composition plugins. Old input file parameters and the
 * Simulator access function 'get_compositional_initial_condition' have been
 * deprecated. Note that there is no longer a default initial composition
 * (previously it was 'function'), but a plugin needs to be specified.
 * <br>
 * (Rene Gassmoeller, 2017/04/23)
 *
 * <li> New: Initial temperature conditions can now be specified as a list
 * of plugins. Their initial temperatures are summed to generate the
 * final initial condition. For this purpose there is a new input file
 * parameter 'List of model names', and a new manager class that owns all
 * of the initial temperature plugins. Old input file parameters and the
 * Simulator access function 'get_initial_condition' have been deprecated.
 * <br>
 * (Rene Gassmoeller, 2017/04/23)
 *
 * <li> Changed: With a sufficiently new deal.II (9.0.0.pre),
 * the spherical shell geometry no longer uses boundary objects,
 * but lets the SphericalManifold class handle all geometry tasks.
 * Using older version of deal.II the manifold is used for mesh refinement,
 * and boundary objects are used otherwise.
 * This results in improved meshes for the side boundaries of octants
 * of 3D spheres in both cases.
 * <br>
 * (Anne Glerum, Rene Gassmoeller, 2017/04/19)
 *
 * <li> Changed: Entries about new contributions are no longer stored
 * in doc/modules/changes.h, but instead each contribution is stored as
 * a single file of the form YYYYMMDD_name in doc/modules/changes/.
 * This creates less conflicts when multiple parallel branches are
 * merged.
 * <br>
 * (Rene Gassmoeller, 2017/04/17)
 *
 * <li> Changed: The visco_plastic material model now uses a single
 * compositional field for strain weakening, which can be calculated
 * with the integrated strain invariant particle property or an
 * equivalent compositional field plugin located in the benchmark
 * buiter_et_al_2008_jgr.
 * <br>
 * (John Naliboff, 2017/04/14)
 *
 * <li> Changed: The folder benchmark/ was renamed into benchmarks/
 * to be consistent with our directory tests/.
 * <br>
 * (Rene Gassmoeller, 2017/04/14)
 *
 * <li> Changed: The global statistics information (time, timestep, degrees of
 * freedom, solver iterations), is now written at one point in a new
 * postprocessor 'global_statistics' instead of wherever they happened
 * to be available. This allows to control the order and format of the
 * columns in the statistics file, and fixes a bug in the statistics file,
 * in which data from later timesteps was written into lines
 * of previous timesteps in case of nonlinear solver schemes. Additionally,
 * the postprocessor by default now outputs only a single line per
 * timestep also for nonlinear solver schemes. To restore the previous
 * behavior (one line per nonlinear iteration) there is a new input
 * parameter "Postprocess/Global statistics/Write statistics for each nonlinear
 * iteration". A consequence of this change is that statistics about
 * the initial refinement steps is only written if the "Mesh refinement/
 * Run postprocessors on initial refinement" parameter is set to true.
 * <br>
 * (Juliane Dannberg, Rene Gassmoeller, 2017/04/14)
 *
 * <li> Changed: aspect now outputs only relative nonlinear residuals of
 * nonlinear iterations in the screen output. In addition, the number of
 * the nonlinear iteration is included in the screen output.
 * <br>
 * (Juliane Dannberg, 2017/04/13)
 *
 * <li> New: ASPECT now supports postprocessing of nonlinear iterations,
 * in particular generating graphical output for each iteration.
 * <br>
 * (Juliane Dannberg, 2017/04/12)
 *
 * <li> New: aspect now provides as script that -- when used together with
 * the deal.II parameter GUI program -- allows for a graphical creation and
 * modification of input parameter files. All available parameters are listed,
 * including their documentation, type, allowed range, and default value. The
 * parameter file written by the GUI will only contain values that are
 * different from the default values, to keep the file easily readable.
 * <br>
 * (Timo Heister, Rene Gassmoeller, Juliane Dannberg, 2017/03/28)
 *
 * <li> New: aspect now supports the --output-xml flag to generate .xml parameter
 * files that can be edited using the deal.II parameter GUI.
 * <br>
 * (Timo Heister, 2017/03/24)
 *
 * <li> New: Compressible computations now work with discontinuous pressure
 * elements.
 * <br>
 * (Timo Heister, 2017/03/23)
 *
 * <li> New: The box geometry can now include initial topography.
 * <br>
 * (Anne Glerum, 2017/03/16)
 *
 * </ol>
 */
