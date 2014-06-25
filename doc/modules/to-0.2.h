/**
 * @page changes_between_0.1_and_0.2 Changes between version 0.1 and version 0.2
 *
 * <p> This is the list of changes made after the release of Aspect version
 * 0.1 and before 0.2. All entries are signed with the names of the author.
 * </p>
 *
 * <ol>
 * <li> Changed: Some models declared parameters that indicate the names of
 * files and that contained underscores. These have been renamed to contain
 * dashes instead, to avoid perpetual difficulties in creating the manual.
 * <br>
 * (Timo Heister, 2013/04/28)
 *
 * <li> Fixed: cleanup the various depth_average functions. Fix that all
 * quantities got averaged the wrong way if the mesh was adaptively refined.
 * Additionally, the velocity_magnitude and sinking_velocity was computed
 * incorrectly.
 * <br>
 * (Timo Heister, 2013/04/28)
 *
 * <li> Fixed: Writing in any output format other than VTU and HDF5 did not
 * yield any output files. This is now fixed.
 * <br>
 * (Wolfgang Bangerth, 2013/04/28)
 *
 * <li> Fixed: There were many places where we indiscriminately used
 * MPI_COMM_WORLD, rather than the communicator used for the actual
 * simulation. This is not a problem almost all the time, except for cases
 * where Aspect is run as part of a bigger simulation.
 * <br>
 * (Wolfgang Bangerth, 2013/04/28)
 *
 * <li> New: All parts of the code now use the new MaterialModel::evaluate().
 * The individual functions like viscosity() got removed. To get an old
 * material model to work, you can either move to the new setup or derive from
 * InterfaceCompatibility instead of Interface.
 * <br>
 * (Timo Heister, 2013/04/28)
 *
 * <li> New: Checking whether the end time has been reached is now a regular
 * termination criterion plugin. Unlike other plugins, however, it is always
 * active.
 * <br>
 * (Wolfgang Bangerth, 2013/04/28)
 *
 * <li> New: The following plugin subsystems are now documented: mesh
 * refinement criteria, termination criteria, and velocity boundary
 * conditions.
 * <br>
 * (Wolfgang Bangerth, 2013/04/28)
 *
 * <li> New: There is now a new top-level parameter, "Number of cheap Stokes
 * solver steps" that allows choosing how many iterations the Stokes solver
 * should attempt using the "cheap" preconditioner described in the ASPECT
 * paper. In particular, a value of zero indicates that the first phase should
 * be skipped altogether. This is a useful strategy for problems with strongly
 * varying viscosity for which the cheap preconditioner does not help very
 * much and the time used for this first phase is often wasted.
 * <br>
 * (Wolfgang Bangerth, 2013/04/22)
 *
 * <li> Fixed: Entropy exponent stabilization_alpha ("alpha") should be an
 * integer, not a double. The recommended setting is the value 2. Note that
 * the value 1 is not explained in the ASPECT paper, only a value of 2. The
 * various other stabilization parameters have also been better explained and
 * cross-referenced to the ASPECT paper.
 * <br>
 * (Timo Heister, 2013/04/17)
 *
 * <li> Fixed: When running an input file that simply forgot to provide names
 * to certain plugins (for example, that does not list anything at all for the
 * geometry model), the resulting error message was pretty obscure and non-
 * informative. This has now been fixed: the error message is much more
 * descriptive and at least says where the error may have happened and what
 * the missing parameter is.
 * <br>
 * (Wolfgang Bangerth, 2013/04/05)
 *
 * <li> New: A postprocessor to output the sparsity pattern of the block
 * matrix. Intended mainly for debugging purposes.
 * <br>
 * (Ted Studley, 2012/12/21)
 *
 * <li> New: a module for specifying termination criteria has been added. The
 * first termination criteria is detecting steady state RMS velocity. This
 * scheme also allows checkpointing the simulation right before termination.
 * <br>
 * (Eric Heien, 2012/12/11)
 *
 * <li> New: The way we do mesh refinement has been completely rewritten. This
 * is now handled via plugins that provide different refinement strategies and
 * that can be combined if so desired. The names of refinement criteria have
 * also been changed: "Temperature", "Velocity", and "Density c_p temperature"
 * have been renamed to "temperature", "velocity", and "thermal energy
 * density". Furthermore, the effect of the old name "Normalized density and
 * temperature" can be obtained by passing the list "density, temperature".
 * For more detail, see the manual's section on parameters that describe the
 * mesh refinement process.
 * <br>
 * (Wolfgang Bangerth, 2012/12/08)
 *
 * <li> Fixed: In nonlinear models, we did not recompute the matrix and
 * preconditioner for the Stokes system in every time step. This is now fixed.
 * <br>
 * (Wolfgang Bangerth, 2012/12/04)
 *
 * <li> Changed: Write the parameter file only on the root node to avoid file
 * system locking issues.
 * <br>
 * (Eric Heien, 2012/11/28)
 *
 * <li> Changed: Adjusted the stabilization parameters to stabilize advection
 * benchmarks without natural diffusion.
 * <br>
 * (Rene Gassmoeller, 2012/11/16)
 *
 * <li> Changed: Unified the routines for advecting temperature and
 * compositional fields. This includes a different scaling for the temperature
 * stabilization that may change results.
 * <br>
 * (Rene Gassmoeller, 2012/11/15)
 *
 * <li> New: If the "Use conduction timestep" parameter is true, the timestep
 * is calculated as the minimum of the convection *and* heat conduction
 * timesteps. This is to allow models where flow is very slow relative to
 * conduction.
 * <br>
 * (Eric Heien, 2012/11/06)
 *
 * <li> New: Aspect now catches if the output directory specified in the
 * parameter file does not exist, rather than providing a cryptic error
 * message upon first attempt to write to this directory.
 * <br>
 * (Wolfgang Bangerth, 2012/10/21)
 *
 * <li> New: Material models can now depend on compositional fields.
 * <br>
 * (Juliane Dannberg, 2012/10/19)
 *
 * <li> New: Introduce compositional fields, that are advected by the
 * velocity.
 * <br>
 * (Juliane Dannberg, 2012/10/19)
 *
 * <li> New: There is now a plugin that allows using GPlates velocity models
 * as input for the velocity boundary conditions.
 * <br>
 * (Rene Gassmoeller, 2012/10/18)
 *
 * <li> Changed: SimulatorAccess moved from aspect::Postproccessor to aspect::
 * and can be used in material/gravity/... models too. If needed, derive from
 * SimulatorAccess.
 * <br>
 * (Timo Heister, 2012/09/24)
 *
 * <li> Fixed: Resume does no longer trigger additional refinement times in
 * the past.
 * <br>
 * (Timo Heister, 2012/09/20)
 *
 * <li> New: Add structures for iterated IMPES and nonlinear solvers.
 * <br>
 * (Markus Buerg, 2012/09/20)
 *
 * <li> New: InitialCondition::Function for the temperature given in the
 * parameter file.
 * <br>
 * (Timo Heister, 2012/09/17)
 *
 * <li> New: BoundaryTemperature::Box is more flexible because it has
 * parameters for the values for all sides of the box.
 * <br>
 * (Timo Heister, 2012/09/17)
 *
 * <li> New: Added parameter 'Temperature solver tolerance' to control the
 * accuracy of the temperature solver.
 * <br>
 * (Timo Heister, 2012/09/05)
 *
 * <li> New: Added MaterialModel::update() and GravityModel::update() that are
 * called before every time step.
 * <br>
 * (Timo Heister, 2012/09/05)
 *
 * <li> New: Support for HDF5/XDMF visualization output.
 * <br>
 * (Eric Heien, 2012/08/28)
 *
 * <li> New: It is now possible to select the tolerance for linear solvers
 * from the parameter file via the global parameter "Linear solver tolerance".
 * <br>
 * (Wolfgang Bangerth, 2012/06/30)
 *
 * <li> New: Tracer particle postprocesser can write in parallel to HDF5.
 * <br>
 * (Eric Heien, 2012/06/08)
 *
 * <li> New: Tracer particle postprocessor also outputs
 * <code>particle.pvd</code> for ParaView. Tracer particle data may also be
 * output to ASCII files for easy parsing and analysis.
 * <br>
 * (Eric Heien, 2012/06/01)
 *
 * <li> New: Aspect now writes a <code>solution.pvd</code> for Paraview that
 * contains a list of the files that jointly make up the entire simulation
 * (and not just a single time step) together with the simulation time each of
 * the files that describe a time step correspond to.
 * <br>
 * (Wolfgang Bangerth, 2012/05/30)
 *
 * <li> New: It is now selectable how the pressure should be normalized (and
 * if it should be normalized at all). If it should be normalized at the end
 * of each time step, one can select if the average pressure at the top
 * surface, or the average pressure throughout the entire domain should be set
 * to a given value (for example to zero).
 * <br>
 * (Wolfgang Bangerth, 2012/05/15)
 *
 * <li> New: The variables we output in graphical format files are now user
 * selectable from the input parameter file. Functions that compute something
 * from velocity, pressure and temperature (e.g., the viscosity, strain rate,
 * etc) are now implemented as plugins like many of the other parts of Aspect.
 * now implemented
 * <br>
 * (Wolfgang Bangerth, 2012/05/08)
 *
 * <li> Modified: Checkpoint frequency is now a user specified parameter.
 * Frequency may be specified as either wall time or number of time steps.
 * <br>
 * (Eric Heien, 2012/05/04)
 *
 * <li> Modified: Background writing of visualization data now uses the
 * mkstemp function rather than forking mktemp, and will continue as normal
 * (at lower efficiency) if temporary files cannot be created.
 * <br>
 * (Eric Heien, 2012/05/03)
 *
 * <li> New: The compressibility that functions implementing the interface
 * aspect::MaterialModel::Interface::compressibility() need to return was
 * previously defined incorrectly as $\frac{\partial\rho}{\partial p}$. This
 * is not what is commonly referred to as compressibility, and the function is
 * now supposed to return $\frac 1\rho \frac{\partial\rho}{\partial p}$
 * instead, following the common definition of the word.
 * <br>
 * Note that all currently implemented compressible models already did that.
 * <br>
 * (Wolfgang Bangerth, Timo Heister, 2012/05/03)
 *
 * <li> New: The number of space dimensions in which a simulation happens is
 * now a parameter that is set in the input parameter filer, rather than
 * statically at compile time as before.
 * <br>
 * (Wolfgang Bangerth, 2012/04/03)
 *
 * <li> New: A postprocessor module for particles was added. This module
 * creates a random uniform distribution of particles in the mesh and moves
 * them based on the velocity field using an RK2 scheme.
 * <br>
 * (Eric Heien, 2012/04/02)
 * </ol>
 *
 *
 */
