/**
 * @page changes_between_2.3.0_and_2.4.0 Changes between version 2.3.0 and version 2.4.0
 *
 * <p> This is the list of changes made after the release of ASPECT version
 * 2.3.0 for version 2.4.0. All entries are signed with the names of the author.
 * </p>
 *
 * <ol>
 *
 * <li> Fixed: The prescribed temperature is now being used in the
 * computation of material properties when assembling the
 * Stokes equations. Before, the Stokes equations used the
 * initial temperature rather than the prescribed temperature
 * because it was only copied into the solution and not the
 * current linearization point.
 * <br>
 * (Juliane Dannberg, 2022/07/06)
 *
 * <li> Fixed: The global volume of the mesh that was stored in the Simulator class was
 * not updated after mesh deformation. This is fixed now.
 * <br>
 * (Rene Gassmoeller, 2022/07/06)
 *
 * <li> New: The defect correction and Newton solvers now also support
 * the projected density method in the mass conservation equation.
 * <br>
 * (Menno Fraters, 2022/06/27)
 *
 * <li> New: There is now a new material model for melting in the
 * lowermost mantle. It can be used to reproduce the results
 * of Dannberg, J., Myhill, R., Gassmöller, R., Cottaar, S.
 * (2021). The morphology, evolution and seismic visibility
 * of partial melt at the core–mantle boundary: implications
 * for ULVZs. Geophysical Journal International, 227(2),
 * 1028-1059.
 * <br>
 * (Juliane Dannberg, 2022/05/24)
 *
 * <li> New: Added an output variable for depth average postprocessor called 'log viscosity' that computes the lateral average of log10(viscosity) in depth layers specified in the postprocessor.
 * <br>
 * (Arushi Saxena, 2022/05/24)
 *
 * <li> New: The DislocationViscosityOutputs now also contain the
 * diffusion viscosity in addition to the dislocation viscosity
 * and the boundary area change work fraction. This is useful
 * to compute the ratio of diffusion and dislocation creep in
 * cases where the viscosity is averaged cell-wise, because in
 * these cases this ratio can not be recovered accurately from
 * just the dislocation and effective viscosity.
 * <br>
 * (Juliane Dannberg, 2022/05/23)
 *
 * <li> New: The block GMG Stokes solver now works for problems with
 * free-surface boundaries and elasticity.
 * <br>
 * (Anne Glerum, Timo Heister, John Naliboff, Jiaqi Zhang, 2022/05/22)
 *
 * <li> New: Add manual documentation for the darcy
 * advection field method
 * <br>
 * (Daniel Douglas, 2023/05/22)
 *
 * <li> New: Add a benchmark for load induced flexure with options
 * for specifying sediment and rock material infilling the
 * flexural moat.
 *
 * <br>
 * (Daniel Douglas, 2022/05/22)
 *
 * <li> New: Add a test case for fixing the temperature in a user
 * specified compositional field to the temperature that the
 * field was initialized with.
 *
 * <br>
 * (Daniel Douglas, 2022/05/21)
 *
 * <li> Fixed: The 'gravity calculation' postprocessor has been made faster
 * by distributing the work to compute theoretical values among the
 * available processors, and due to other changes.
 * <br>
 * (Wolfgang Bangerth, 2022/05/20)
 *
 * <li> Changed: The reference viscosity function is no longer part
 * of the material model interface (because we do not require
 * it anymore for the pressure scaling or anywhere else). The
 * one material model affected by this is the depth-dependent
 * model, which now has its own reference viscosity parameter
 * rather than use the reference viscosity of the base model.
 * <br>
 * (Juliane Dannberg, 2022/05/20)
 *
 * <li> New: The `read_and_distribute_file_content` utility function can now also read
 * binary files compressed with gzip, which allows for example the StructuredData
 * class (and therefore all ascii data plugins) to read gzip compressed data.
 * Compressed files are detected based on their file ending '.gz'. I.e. if you want
 * to read a compressed file, make sure its name ends in '.gz', and if you want to
 * read a plain text file, make sure its name does not end in '.gz'. Compressed
 * files can save around 75% to 80% of disk space for typical ascii data files.
 * <br>
 * (Rene Gassmoeller, 2022/05/20)
 *
 * <li> Changed: Removed the reference viscosity functions from all
 * plugins. In plugins that did not use the "Reference
 * viscosity" input parameter, this parameter was removed as
 * well and can no longer be used in input files.
 * <br>
 * (Juliane Dannberg, 2022/05/20)
 *
 * <li> New: Visualization postprocessors now record the physical units of the
 * quantity they compute, and this information is also output into visualization
 * files with a sufficiently new version of deal.II.
 * <br>
 * (Wolfgang Bangerth, 2022/05/19)
 *
 * <li> New: Add an advection field method that
 * advects a compositional field according
 * to Darcy's Law.
 * <br>
 * (Daniel Douglas, 2022/05/17)
 *
 * <li> Changed: The global limiter for the quadratic least squares particle interpolation
 * plugin was replaced by a bounds preserving slope limiter that respects
 * local bounds on each cell. The global limiter options have been removed
 * from the parameter file.
 * <br>
 * (Mack Gregory, 2022/05/16)
 *
 * <li> New: Added a visualization postprocessor called 'boundary strain rate residual' that computes the residual of the strain rate invariant at the top boundary of the model domain. The residual is computed as the difference between the modeled strain rate and a reference strain rate that is read from an ascii data file, e.g. obtained by field measurements.
 * <br>
 * (Arushi Saxena, 2022/05/16)
 *
 * <li> New: The viscoplastic strain invariants particle
 * property can now also track noninitial plastic strain.
 * <br>
 * (John Naliboff, 2022/05/05)
 *
 * <li> New: There is now a cookbook of kinematically driven
 * oceanic subduction in 2D with isoviscous materials and
 * without temperature effects. The cookbook model setup
 * is based on Quinquis (2014).
 * <br>
 * (Anne Glerum, 2022/04/21)
 *
 * <li> New: There is now a new advection method for the temperature
 * called 'prescribed field with diffusion', which works in the
 * same way as the corresponding advection field method for
 * compositional fields.
 * <br>
 * (Juliane Dannberg, 2022/03/28)
 *
 * <li> Fixed: The entropy adiabat benchmark did not include
 * thermal conduction along the adiabat. This has been fixed
 * by splitting the energy equation into two parts: Conduction
 * is solved using the temperature field and the 'prescribed
 * field with diffusion' method, and all remaining terms are
 * included when solving the equation for the entropy (as a
 * compositional field). To show that this works correctly,
 * there are now three more variations of the benchmark,
 * modeling conduction starting from an adiabatic profile,
 * conduction of a half-space without adiabatic heating, and
 * a version of the entropy adiabat benchmark with increased
 * conduction (all using this new method).
 * <br>
 * (Juliane Dannberg, 2022/03/28)
 *
 * <li> Changed: The minimum and maximum refinement functions
 * now evaluate the center of each mesh cell to determine
 * whether to coarsen or refine a cell instead of evaluating
 * all vertices.
 * <br>
 * (Anne Glerum, 2022/03/17)
 *
 * <li> New: The mesh displacements are computed by the
 * GMG solver when the GMG Stokes solver is used.
 * <br>
 * (Timo Heister, Jiaqi Zhang, 2022/03/01)
 *
 * <li> Fixed: The scaling factor in the gplates boundary velocity model is now
 * applied correctly. Before, it was always taken to be 1, so it did not
 * scale the velocities at all.
 * <br>
 * (Juliane Dannberg, 2022/02/22)
 *
 * <li> Removed: The input parameter 'First data file model time' of the 'ascii data'
 * boundary condition plugin has been removed. It was not well defined what
 * the plugin returned before this time was reached.
 * <br>
 * (Rene Gassmoeller, 2022/02/17)
 *
 * <li> Changed: The global limiter for the least squares particle interpolation
 * plugin was replaced by a bounds preserving slope limiter that respects
 * local bounds on each cell. The global limiter options have been removed
 * from the parameter file options.
 * <br>
 * (Mack Gregory, Gerry Puckett, 2022/02/08)
 *
 * <li> Changed: The initial composition model called 'ascii data' can now
 * read in 3d ascii datasets into a 2d model and slice the dataset
 * in a user controlled plane. This
 * allows it to make high-resolution 2d models of problems that use
 * observational data (such as seismic tomography models).
 * <br>
 * (Juliane Dannberg, Rene Gassmoeller, 2022/02/07)
 *
 * <li> Fixed: The GMG solver could fail in problems with large viscosity variations
 * and Q1 viscosity averaging, because viscosities would become negative
 * on coarser GMG levels. This is now fixed.
 * <br>
 * (Timo Heister, Jiaqi Zhang, 2022/02/05)
 *
 * <li> Fixed: When using the full A block as Stokes preconditioner in a compressible
 * model using a Newton solver scheme we added an unnecessary assembler, which
 * would throw an assertion. This is fixed now by not adding the assembler.
 * <br>
 * (Rene Gassmoeller, 2022/02/02)
 *
 * <li> Changed: The Steinberger model now uses the hydrostatic pressure
 * rather than the full pressure to compute the material properties from
 * the thermodynamic lookup equation of state (density, etc.). This is
 * consistent with the projected density approximation and avoids
 * convergence issues of the nonlinear solver related to jumps of material
 * properties at phase transitions, which are common in thermodynamic
 * lookup tables. At the same time, using the full pressure does not really
 * make the results more accurate since the dynamic pressure is so small
 * that it changes material properties very little otherwise (not more than
 * 0.1% even in the most extreme cases).
 * <br>
 * (Juliane Dannberg, 2022/01/28)
 *
 * <li> New: There is now an averaging operation "geometric average only viscosity"
 * that can be used to average the viscosity values of the different
 * points in a cell. This can be used together with the GMG solver.
 * <br>
 * (Juliane Dannberg, 2022/01/28)
 *
 * <li> Fixed: The compressible terms for the DC Stokes and Newton solvers
 * were not correctly assembled. This is fixed now, and tests comparing
 * the Stokes and DC Stokes solver have been added.
 * <br>
 * (Menno Fraters, 2022/01/26)
 *
 * <li> Fixed: Compressible models with periodic boundaries that did not also
 * have one or more open boundaries did not correctly apply the necessary
 * pressure right-hand side compatibility modification when solving the
 * Stokes equations. This generally caused the linear solver to fail or to
 * take many more iterations than necessary. The right-hand side pressure
 * compatibility modification is now applied correctly, fixing this problem
 * and allowing it to use periodic boundaries in compressible models.
 * <br>
 * (Juliane Dannberg, 2022/01/26)
 *
 * <li> Fixed: The simple version of the steinberger viscosity profile data
 * file now correctly excludes both the top and bottom boundary layer.
 * These need to be removed from the original profile in case the
 * model uses a temperature-dependent viscosity based on the adiabatic
 * profile (rather than the laterally averaged temperature). Otherwise
 * the temperature effect will be included twice. The old version of
 * the file had excluded the top, but not the bottom boundary layer.
 * <br>
 * (Juliane Dannberg, 2022/01/24)
 *
 * <li> New: The steinberger material model now includes the option to use
 * a pressure- and temperature-dependent thermal conductivity.
 * <br>
 * (Juliane Dannberg, 2022/01/21)
 *
 * <li> New: The matrix-free GMG Stokes preconditioner is now
 * implemented for the free surface stabilization.
 * <br>
 * (Timo Heister, John Naliboff, Jiaqi Zhang, 2022/01/11)
 *
 * <li> Fixed: Boundary conditions are no longer applied to the projected
 * density field.
 * <br>
 * (Juliane Dannberg, 2022/01/07)
 *
 * <li> Changed: The Steinberger material model can now be used with the
 * projected density approximation (and with an arbitrary number of
 * compositional fields).
 * New: There is now a new type of compositional field called
 * 'density' that is intended to be used with the projected density
 * approximation.
 * <br>
 * (Juliane Dannberg, 2022/01/07)
 *
 * <li> Fixed: Boundaries with tangential velocity are no longer taken into
 * account when checking if a boundary is an outflow boundary, which is
 * important in the case where temperature or composition are only fixed
 * on inflow (but not outflow) parts of the boundary. The old behavior
 * led to temperature/composition sometimes not being prescribed
 * correctly, especially in the case of spherical geometries.
 * <br>
 * (Juliane Dannberg, 2021/07/07)
 *
 * <li> Fixed: Change assert throw for radiogenic heating such that it should only be evaluated if crust is defined by compostion.
 * <br>
 * (Elodie Kendall, 2022/01/03)
 *
 * <li> New: The variables "minimum viscosity" and "maximum viscosity" in the visco plastic material model can now have different values for different phases. This is useful in assigning shear zone viscosity in a subduction model.
 * <br>
 * (Haoyuan Li, 2021/12/21)
 *
 * <li> Fixed: The algorithm in the pressure normalization scheme 'surface' ignored
 * some surface faces of cells with much larger horizontal than vertical extent in
 * curved geometries. This could lead to non-zero averaged surface pressure for
 * these models. This is fixed now.
 * <br>
 * (Rene Gassmoeller, 2021/11/18)
 *
 * <li> New: Added a new postprocessor which computes the parameter "Mobility" following Lourenco et al., 2020 G3
 * <br>
 * (Elodie Kendall, Rene Gassmoeller, Anne Glerum and Bob Myhill, 2021/11/11)
 *
 * <li> Fixed: The 'random uniform' particle generator plugin would not consider the
 * volume of each cell correctly when generating a uniform particle density. This
 * was irrelevant for box models, but would lead to too many particles in cells
 * close to the bottom of curved models, and too few particles in cells close to
 * the top of curved models. Box models with initial adaptive refinement were not
 * affected by this bug. This is fixed now.
 * <br>
 * (Rene Gassmoeller, 2021/11/10)
 *
 * <li> Fixed: The spherical velocity postprocessor now correctly applies the parameter 'Use years in output instead of seconds' to output velocities either per year or per second.
 * <br>
 * (Elodie Kendall and Rene Gassmoeller, 2021/10/29)
 *
 * <li> New: There is now a postprocessor that computes the
 * second invariant of the deviatoric stress tensor.
 * <br>
 * (Anne Glerum, 2021/10/22)
 *
 * <li> Changed: The right-hand-side force term in the momentum
 * equation is now also computed in the first timestep (t0),
 * such that nonzero initial viscoelastic stresses can be set.
 * <br>
 * (Anne Glerum, 2021/09/27)
 *
 * <li> Changed/New: The material model dynamic_friction has been integrated
 * into a new rheology model friction_models that can be used together
 * with the visco_plastic material model.
 * <br>
 * (Esther Heckenbach, 2021/09/23)
 *
 * <li> New: Added a print statement into log file which tells the user the total wallclock time including restarts.
 * <br>
 * (Elodie Kendall, Timo Heister and Rene Gassmoeller, 2021/09/20)
 *
 * <li> New: ASPECT now has a cookbook which uses the gravity postprocessor to
 * compute gravity generated by S40RTS-based mantle density variations.
 * <br>
 * (Cedric Thieulot, 2021/08/31)
 *
 * <li> New: The "Compositional fields" subsection in ASPECT has a new parameter called
 * "Types of fields" that is used to specify the types of all the
 * "compositional" fields. Supported types are currently
 * "chemical composition", "stress", "grain size", "porosity", "generic"
 * and "unspecified". The parameter is used to fill a vector of
 * CompositionalFieldDescription objects that store metadata about each field.
 * This vector can be called from within material models
 * via this->introspection().get_field_descriptions().
 * <br>
 * (Bob Myhill, 2021/08/31)
 *
 * <li> Changed: We now require CMake version 3.1.0 or newer.
 * <br>
 * (Timo Heister, 2021/08/19)
 *
 * <li> New: The GMG solver now supports boundaries with mesh deformation,
 * although free surface boundaries are still not supported.
 * <br>
 * (Timo Heister, 2021/08/05)
 *
 * <li> New: Where possible, when using large data tables as input (e.g., for
 * initial conditions specified as tables), these data are now stored
 * only once on each node in memory areas that is accessible by all MPI
 * processes on that node.
 * <br>
 * (Wolfgang Bangerth, 2021/07/28)
 *
 * <li> New: The GMG solver now supports masked velocity boundary conditions specified as [xyz]:.
 * <br>
 * (Timo Heister, 2021/07/16)
 *
 * <li> Improved: The continental extension cookbook has been updated to
 * use a number of new features, including a faster solver scheme,
 * a more realistic method for initiating deformation, and
 * initial adaptive mesh refinement.
 * <br>
 * (John Naliboff, 2021/07/16)
 *
 * <li> New: There is now a cookbook that visualizes the phase diagram from results of a model run.
 * This includes examples from the Visco-Plastic and Steinberger material model.
 * <br>
 * (Haoyuan Li and Magali Billen, 2021/07/16)
 *
 * <li> New: The matrix-free GMG Stokes preconditioner is now
 * implemented for the Newton solver.
 * <br>
 * (Timo Heister, Menno Fraters, Jiaqi Zhang, 2021/07/15)
 *
 * <li> New: There is now a cookbook that reproduces convection models
 * with a phase function from Christensen and Yuen, 1985.
 * <br>
 * (Juliane Dannberg, 2021/07/14)
 *
 * <li> Changed: The latent heat material model is now consistent with the
 * density in the phase function formulation from Christensen and Yuen, 1985.
 * <br>
 * (Juliane Dannberg, 2021/07/14)
 *
 * <li> New: There is now a rising velocity output variable in the depth
 * average postprocessor.
 * <br>
 * (Juliane Dannberg, 2021/07/13)
 *
 * <li> New: CMake now detects and reports an error if the user tries to do an in-source ASPECT build. Please create a separate build directory before running cmake.
 * <br>
 * (Rene Gassmoeller, Timo Heister, 2021/07/13)
 *
 * <li> New: ASPECT now has a cookbook which shows how velocities can be
 * prescribed at positions specified by an ASCII input file.
 * <br>
 * (Bob Myhill, 2021/07/13)
 *
 * <li> Fixed: The sinking velocity depth average postprocessor now computes
 * the actual sinking velocity (and not the rising velocity).
 * <br>
 * (Juliane Dannberg, 2021/07/13)
 *
 * <li> New: The geoid postprocessor can now handle a deforming mesh (free surface), next to the already existing option from the dynamic topography postprocessor output.
 * <br>
 * (Maaike Weerdesteijn, Rene Gassmoeller, Jacky Austermann, 2021/07/12)
 *
 * <li> Changed: The "surface stress" visualization postprocessor used to just
 * output a number of independent components. Like the "stress" and
 * "strain rate tensor" postprocessors, it now outputs these as a tensor field.
 * <br>
 * (Wolfgang Bangerth, 2021/07/11)
 *
 * <li> New: The thermodynamic lookup equation of state can now read in
 * a column that contains the name of the dominant phase from the
 * data tables.
 * <br>
 * (Juliane Dannberg, 2021/07/10)
 *
 * <li> Fixed: The "prescribed_dilation" and "additional Stokes RHS" options now work correctly with the Newton solver.
 * <br>
 * (Sibiao Liu, 2021/07/09)
 *
 * <li> New: There is now a postprocessor that computes the maximum depth
 * of each compositional field at each timestep.
 * <br>
 * (Anne Glerum, 2021/07/09)
 *
 * <li> New: There is now a geometry model plugin "chunk with lithosphere boundary indicators"
 * that allows for two different types of boundary conditions on the side boundaries
 * of a chunk domain. On each side boundary, two boundary indicators are available
 * to set these different boundary conditions. A use-case can be prescribed plate motions
 * on the upper part of the domain and an open boundary underneath the plates.
 * <br>
 * (Anne Glerum, 2021/07/09)
 *
 * <li> Changed: The Boundary temperature model plugin now requires that the
 * "Fixed temperature boundary indicators" parameter is non-empty whenever
 * the "Model names" parameter is not empty.
 * <br>
 * (Bob Myhill, 2021/07/09)
 *
 * <li> Improved: Particle operations have been optimized for speed. Particles now
 * always carry a hidden internal property that stores properties necessary for
 * their advection scheme. These properties are not written into output but are
 * included if the particles are asked for their properties internally.
 * <br>
 * (Rene Gassmoeller, 2021/07/08)
 *
 * <li> New: There is now a ‘static’ option for the temperature field that is set-up similarly to the ‘static’ option for compositional fields. This allows the temperature field to be static while advection is on so you can still advect and build up elastic stresses.
 * <br>
 * (Rebecca Fildes, Magali Billen, 2021/07/08)
 *
 *
 * <li> Changed: When using the constant temperature boundary plugin, ASPECT
 * now checks that the Fixed temperature boundary indicators match the
 * indicators in the model subsection.
 * <br>
 * (Bob Myhill, 2021/07/09)
 *
 * <li> New: ASPECT now has a ThermodynamicTableLookup
 * equation of state plugin. This plugin allows material models
 * to read in one or more Perple_X or HeFESTo table files,
 * interpolate material properties at desired pressures and
 * temperatures, and use the interpolated properties as
 * material model outputs. The equation of state plugin is
 * currently used in the Steinberger material model.
 * <br>
 * (Bob Myhill, 2021/07/08)
 *
 * <li> New: It is now possible to use default or single values for Peierls creep
 * parameters in composite (viscoplastic) rheologies with phase transitions.
 * <br>
 * (Bob Myhill, 2021/07/08)
 *
 * <li> New: There is now a termination criterion based on the steady state
 * heat flux.
 * <br>
 * (Juliane Dannberg, 2021/07/07)
 *
 * <li> New: Tests can now depend on another test if a special tag 'DEPENDS-ON:'
 * is specified in the test parameter file.
 * <br>
 * (Rene Gassmoeller, 2021/07/07)
 *
 * <li> New: Added a new input parameter Elastic damper viscosity to Rheology::Elasticity.
 * This parameter corresponds to the viscosity of a viscous damper which deforms
 * at the same strain rate as the elastic element.
 * The default value of 0 Pas corresponds to purely elastic deformation.
 * <br>
 * (Bob Myhill, 2021/07/07)
 *
 * <li> New: It is now possible to only call the world builder to determine initial
 * compositions for selected compositional fields by specifying the parameter
 * 'List of relevant compositions' in the world builder initial composition
 * plugin.
 * <br>
 * (Rene Gassmoeller, Juliane Dannberg, 2021/06/25)
 *
 * <li> New: Added a cutoff stress to Rheology::PeierlsCreep.
 *
 * In parameterizations of the Peierls creep flow law where the
 * power law stress exponent is equal to zero, the strain rate
 * does not approach zero as the stress drops to zero. This is
 * a problem, because the iterative solve for the equilibrium
 * stress may yield a negative stress at low strain rates
 * (which results in a negative viscosity).
 *
 * The Peierls creep rheology module now includes a parameter
 * "Cutoff stresses for Peierls creep". At stresses below the
 * cutoff, the strain rate is modelled as a quadratic function
 * of the stress (edot_ii = astress^2 + bstress). This effectively
 * means that Peierls creep transitions into power law creep and
 * then a linear rheology as stress decreases below the cutoff.
 *
 * <br>
 * (Daniel Douglas, 2021/03/23)
 *
 * </ol>
 */
