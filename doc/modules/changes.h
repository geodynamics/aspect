/**
 * @page changes_current Changes after the latest release (v1.1)
 *
 * <p> This is the list of changes made after the release of Aspect version
 * 1.1. All entries are signed with the names of the author. </p>
 *
 *
 * <ol>
 *
 * <li> Changed: Updated maximum refinement function with the same
 * bugfix as the minimum refinement function. It now always declares
 * dim number of variables. When using the 'depth' coordinate system,
 * only the first variable is used. The others are set to zero.
 * <br>
 * (Jonathan Perry-Houts, 2014/08/04)
 *
 * <li> Fixed: The viscosity mesh refinement criterion did not ask
 * the material model for the viscosity, therefore it did not work
 * with many material models. Some material models calculate the
 * viscosity anyway, so it worked in that cases. Now the criterion 
 * asks correctly for the viscosity.
 * <br>
 * (Rene Gassmoeller, Juliane Dannberg 2014/08/04)
 *
 * <li> Changed: The minimum refinement function now always declares
 * three variables. In the 'depth' option only the first one is
 * set to the depth, the other variables are set to zero.
 * <br>
 * (Rene Gassmoeller, 2014/08/04)
 *
 * <li> New: There are now function parser plugins for the
 * gravity model and boundary temperature model. Gravity can
 * be prescribed in dependence of position and boundary temperature
 * in dependence of position and time.
 * <br>
 * (Sanja Panovska, Rene Gassmoeller, 2014/08/01)
 *
 * <li> New: There is now a postprocessor that outputs multiple 
 * material properties with just a single call to the material
 * model. This is more efficient, but only matters for complex
 * material models.
 * <br>
 * (Rene Gassmoeller, 2014/08/01)
 * 
 * <li> New: There is now a section in the manual detailing the
 * "Burstedde" benchmark and its verification.
 * <br>
 * (Iris van Zelst, Wolfgang Bangerth, 2014/07/07)
 *
 * <li> New: Added "maximum refinement function" plugin to limit the depth of
 * the mesh octree.
 * <br>
 * (Jonathan Perry-Houts, 2014/07/03)
 *
 * <li> Changed: The viscosity in the steinberger material model can now be
 * calculated by taking either the lateral averaged temperature or the
 * adiabatic temperature in this depth as reference.
 * <br>
 * (Rene Gassmoeller, 2014/06/20)
 *
 * <li> Changed: The minimum refinement function is now evaluated at every
 * cell corner instead of the cell midpoint. This especially helps in refining
 * thin boundary layers.
 * <br>
 * (Rene Gassmoeller, 2014/06/20)
 *
 * <li> Changed: The internal heating statistics postprocessor now calculates
 * the integral and averaged internal heating by mass averaging instead
 * of volume averaging.
 * <br>
 * (Rene Gassmoeller, 2014/06/20)
 *
 * <li> New: Added a visualization processor to output the local gravity vector.
 * <br>
 * (Timo Heister, 2014/06/17)
 *
 * <li> Fixed: The direct solver now works with nonlinear schemes like
 * iterated IMPES.
 * <br>
 * (Timo Heister, 2014/06/11)
 *
 * <li> Changed: The maximum number of solver iterations for all linear solvers
 * (including the inner solvers) is now set to 1000.
 * <br>
 * (Timo Heister, 2014/06/08)
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
