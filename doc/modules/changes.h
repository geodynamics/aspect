/**
 * @page changes_current Changes after the latest release (v1.2)
 *
 * <p> This is the list of changes made after the release of Aspect version
 * 1.2. All entries are signed with the names of the author. </p>
 *
 *
 * <ol>
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
 * material models was used incorrectly to calculate the density. This did only
 * affect compressible models with these material models and is fixed now.
 * <br>
 * (Rene Gassmoeller, 2015/02/12)
 *
 * <li> Changed: We use a new method to calculate the initial residuals 
 * used for the iterated IMPES solver. Instead of using the first solution,
 * we use the norm of the right-hand side, so that the residual is small
 * for cases where the initial guess of a time step is already very good, 
 * and the number of iterations is decreased for these cases (e.g. when
 * using a smaller time step size).
 * Moreover, the screen output of the residuals is now in the same order as 
 * the output before (about which system is solved and how many iterations were 
 * needed).  
 * <br>
 * (Juliane Dannberg, 2015/02/09)
 *
 * <li> New: A new plugin called 'ascii data' is available for boundary and
 * initial conditions that reads in text data from files provided by the user.
 * The plugin works for box and shell geometries and allows for time-dependent
 * or constant boundary conditions. The data files must provide data for the
 * whole model domain and the plugin interpolates the data from a structured
 * grid to Aspect's mesh linearly.
 * <br>
 * (Eva Bredow, Rene Gassmoeller, 2015/02/03)
 *
 *
 * <li> New: There is now a mass flux statistics postprocessor that
 * calculates the mass flux through every boundary interface. This 
 * is helpful for regional models with prescribed boundary velocities. 
 * In general one would want to set up models with as much in- as out-
 * flux, but in some cases it is hard to know, if the boundary
 * plugins actually ensure this. With this plugins the net flux can
 * be monitored.
 * <br>
 * (Rene Gassmoeller, 2015/02/02)
 *
 * <li> Fixed: There was a small bug in ASPECT 1.2 that only occured in
 * the adiabatic initial temperature plugin, when no temperature
 * was prescribed at any boundary, the model was compressible and a
 * bottom thermal boundary layer was included. This is fixed now.
 * <br>
 * (Rene Gassmoeller, 2015/02/02)
 *
 * <li> Fixed: The residual of the iterated IMPES solver used to be 
 * nan when the initial residual of one of the solution variable was
 * zero. This is fixed now, so that this variable is not considered
 * for the overall residual anymore.
 * <br>
 * (Juliane Dannberg, 2015/01/30)
 *
 * <li> Fixed: The Steinberger material model had a bug in case a
 * material table was used that had not the same number of points in
 * temperature and pressure direction. The table was truncated in pressure
 * dimension to the number of lines of the temperature dimension. This is
 * fixed now.
 * <br>
 * (Rene Gassmoeller, 2015/01/29)
 *
 * </ol>
 */
