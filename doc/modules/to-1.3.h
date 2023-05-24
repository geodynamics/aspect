/**
 * @page changes_between_1.2_and_1.3 Changes between version 1.2 and version 1.3
 *
 * <p> This is the list of changes made after the release of ASPECT version
 * 1.2 for version 1.3. All entries are signed with the names of the author.
 * </p>
 *
 * <ol>
 *
 * <li> Fixed: We now use the correct residual to calculate the stabilization
 * for the advection of temperature and compositional fields. Both the latent
 * heat term for the temperature and the reaction terms for the compositional
 * fields were missing.
 * <br>
 * (Juliane Dannberg, 2015/05/13)
 *
 * <li> New: There is now a mesh refinement criterion called Boundary based on
 * refining near a user-specified list of boundary indicators.
 * <br>
 * (Ian Rose, 2015/05/12)
 *
 * <li> Fixed: The advection solver had problems converging for constant
 * compositional fields in timestep 1. This is fixed now.
 * <br>
 * (Rene Gassmoeller, 2015/05/11)
 *
 * <li> New: There are now visualization postprocessors that can output the
 * shear stress and full stress tensors.
 * <br>
 * (Wolfgang Bangerth, 2015/05/02)
 *
 * <li> Changed: ASPECT accidentally used the type <code>unsigned int</code>
 * in denoting boundary indicators in
 * BoundaryComposition::Interface::composition() and
 * BoundaryTemperature::Interface::temperature(), as well as derived classes'
 * corresponding functions. It should have been the deal.II type
 * types::boundary_id. This is now fixed. It may require fixing user plugins
 * that overload these functions.
 * <br>
 * (Wolfgang Bangerth, 2015/04/30)
 *
 * <li> New: ASPECT can now average material properties between the quadrature
 * points of a cell. This greatly increases the stability of solutions in
 * simulations with spatially varying coefficients, and also greatly
 * accelerates the solution, at times up to a factor of ten.
 * <br>
 * (Wolfgang Bangerth, 2015/04/25)
 *
 * <li> Fixed: Removal rigid body translations and rotations when the
 * simulation has a nullspace was bugged. This feature has been fixed and
 * extended, and a new section has been added to the manual describing its
 * use.
 * <br>
 * (Ian Rose, 2015/04/23)
 *
 * <li> New: vtu visualization output can now be grouped into an arbitrary
 * number of files per time step by setting "Number of grouped files" to a
 * number larger than 1.
 * <br>
 * (Timo Heister, 2015/04/15)
 *
 * <li> Fixed: To make the right-hand side of the Stokes equation compatible
 * to the matrix we need to apply a correction for imbalanced in-/outflow
 * across the model boundaries. This correction was accidentally applied twice
 * in the first iteration of the iterated IMPES solver. This does not change
 * the results, because subsequent iterations will do it correctly, but it
 * might prevent the model to converge or slow down convergence. This is fixed
 * now by applying the correction directly after assembling the matrix and
 * right-hand side.
 * <br>
 * (Rene Gassmoeller, Timo Heister, Juliane Dannberg, 2015/04/14)
 *
 * <li> New: ASPECT can now read the list of input arguments from the default
 * input stdin, instead of a file, by specifying "--" for the name of the
 * input file.
 * <br>
 * (Wolfgang Bangerth, 2015/03/20)
 *
 * <li> Changed: We use the current linearization point instead of the old
 * solution in the assembly of the composition advection system now, and
 * update the current linearization point with the solution after we have
 * solved all of the fields. This allows for model parameters used in the
 * composition advection (like the reaction term) to depend on solution
 * variables and to be updated during nonlinear iterations.
 * <br>
 * (Juliane Dannberg, 2015/02/20)
 *
 * <li> Changed: The default values for the latent heat release of melting in
 * the 'latent heat melt' material model had the wrong sign. This likely
 * resulted from a change in the latent heat terms in the assembly long ago.
 * All tests were correctly set up, but the default values were forgotten.
 * This is fixed now.
 * <br>
 * (Rene Gassmoeller, 2015/02/18)
 *
 *
 * <li> Changed: The unused parameter 'Activation enthalpies' was removed from
 * the 'latent heat' material model. The active parameter for the same purpose
 * is 'Thermal viscosity exponent'. If there are parameter files specifying
 * the removed parameter it is safe to just remove the line.
 * <br>
 * (Rene Gassmoeller, 2015/02/16)
 *
 * <li> Fixed: The compressibility in the 'latent heat' and 'latent heat melt'
 * material models was used incorrectly to calculate the density. This did
 * only affect compressible models with these material models and is fixed
 * now.
 * <br>
 * (Rene Gassmoeller, 2015/02/12)
 *
 * <li> Changed: We use a new method to calculate the initial residuals used
 * for the iterated IMPES solver. Instead of using the first solution, we use
 * the norm of the right-hand side, so that the residual is small for cases
 * where the initial guess of a time step is already very good, and the number
 * of iterations is decreased for these cases (e.g. when using a smaller time
 * step size). Moreover, the screen output of the residuals is now in the same
 * order as the output before (about which system is solved and how many
 * iterations were needed).
 * <br>
 * (Juliane Dannberg, 2015/02/09)
 *
 * <li> New: A new plugin called 'ascii data' is available for boundary and
 * initial conditions that reads in text data from files provided by the user.
 * The plugin works for box and shell geometries and allows for time-dependent
 * or constant boundary conditions. The data files must provide data for the
 * whole model domain and the plugin interpolates the data from a structured
 * grid to ASPECT's mesh linearly.
 * <br>
 * (Eva Bredow, Rene Gassmoeller, 2015/02/03)
 *
 *
 * <li> New: There is now a mass flux statistics postprocessor that calculates
 * the mass flux through every boundary interface. This is helpful for
 * regional models with prescribed boundary velocities. In general one would
 * want to set up models with as much in- as out- flux, but in some cases it
 * is hard to know, if the boundary plugins actually ensure this. With this
 * plugins the net flux can be monitored.
 * <br>
 * (Rene Gassmoeller, 2015/02/02)
 *
 * <li> Fixed: There was a small bug in ASPECT 1.2 that only occurred in the
 * adiabatic initial temperature plugin, when no temperature was prescribed at
 * any boundary, the model was compressible and a bottom thermal boundary
 * layer was included. This is fixed now.
 * <br>
 * (Rene Gassmoeller, 2015/02/02)
 *
 * <li> Fixed: The residual of the iterated IMPES solver used to be nan when
 * the initial residual of one of the solution variable was zero. This is
 * fixed now, so that this variable is not considered for the overall residual
 * anymore.
 * <br>
 * (Juliane Dannberg, 2015/01/30)
 *
 * <li> Fixed: The Steinberger material model had a bug in case a material
 * table was used that had not the same number of points in temperature and
 * pressure direction. The table was truncated in pressure dimension to the
 * number of lines of the temperature dimension. This is fixed now.
 * <br>
 * (Rene Gassmoeller, 2015/01/29)
 *
 * </ol>
 */
