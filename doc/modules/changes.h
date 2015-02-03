/**
 * @page changes_current Changes after the latest release (v1.2)
 *
 * <p> This is the list of changes made after the release of Aspect version
 * 1.2. All entries are signed with the names of the author. </p>
 *
 *
 * <ol>
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
