/**
 * @page changes_current Changes after the latest release (v1.2)
 *
 * <p> This is the list of changes made after the release of Aspect version
 * 1.2. All entries are signed with the names of the author. </p>
 *
 *
 * <ol>
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
