/**
 * @page changes_current Changes after the latest release (v1.1)
 *
 * <p> This is the list of changes made after the release of Aspect version
 * 1.1. All entries are signed with the names of the author. </p>
 *
 *
 * <ol>
 *
 * <li> New: Initial temperature condition "S40RTS perturbation" has been
 * added, that provides a perturbation for a shell geometry based on the
 * S20RTS or S40RTS global shear wave velocity model data.
 * <br>
 * (Jacky Austermann, 2014/06/06)
 *
 * <li> Fixed: When a linear solver fails to converge in a parallel program,
 * every processor would output the same error message -- leading to incredible
 * amounts of entangled error messages. This has now been resolved: every processor
 * still fails, but only processor 0 reports the error.
 * <br>
 * (Wolfgang Bangerth, 2014/06/06)
 *
 * <li> Fixed: When setting "Use years in output instead of seconds" the
 * velocity solution is now exported in m/year instead of m/s in visualization
 * files.
 * <br>
 * (Timo Heister, 2014/06/02)
 *
 * </ol>
 */
