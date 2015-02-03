/**
 * @page changes_between_1.1_and_1.2 Changes between version 1.1 and version 1.2
 *
 * <p> This is the list of changes made after the release of Aspect version
 * 1.1 for version 1.2. All entries are signed with the names of the author.
 * </p>
 *
 * <ol>
 * <li> The default for the number of zones in "depth average" got reduced
 * from 100 to 10. A large number of zones causes erroneous output on coarse
 * meshes.
 * <br>
 * (Timo Heister, 2015/01/19)
 *
 * <li> New: ASPECT now uses meshes for the spherical shell where each node is
 * placed on concentric shells. Previously, upon mesh refinement, new vertices
 * were placed at an averaged location of the vertices of the mother cell,
 * leading to a certain degree of distortion of cells. It also led to a
 * situation where not all cells are equally shaped due to this distortion.
 * The new mesh, in contrast, is completely symmetric. Addtionally, a higher
 * order mapping to represent curved faces is now also used for the interior
 * cells using the manifold description of the shell/sphere.
 * <br>
 * (Wolfgang Bangerth, Timo Heister, 2015/01/16)
 *
 * <li> New: GeometryModel::has_curved_elements() will allow for optimizations
 * if all boundaries (and manifold description) are straight.
 * <br>
 * (Timo Heister, 2015/01/16)
 *
 * <li> Fixed: the inner solvers for the Stokes block preconditioner now start
 * with a zero initial guess. This should make the solvers more robust.
 * <br>
 * (Timo Heister, 2015/01/16)
 *
 * <li> ASPECT now requires deal.II version 8.2 or later.
 * <br>
 * (Timo Heister, 2015/01/13)
 *
 * <li> New: The parameters for the inner Algebraic Multigrid preconditioner
 * used in the Stokes solve have been tuned to achieve a speedup of the solver
 * phase by a factor between 1.2 and 2, depending on application.
 * <br>
 * (Wolfgang Bangerth, 2014/12/05)
 *
 * <li> New:  There is now the option to specify an origin when using the box
 * goemetry model.
 * <br>
 * (Ryan Grove, 2014/12/09)
 *
 * <li> Changed: The behaviour when one changed the visualization
 * output_interval during a checkpoint restart was previously undefined, and
 * working in slightly unexpected ways like never writing output for the first
 * timestep after the restart. This now works as one would expect, e.g. every
 * timestep that ends more than output_interval after the last output time
 * step will produce output. Old checkpoint files will continue to work, only
 * with a possible short gap in visualization output right after restart.
 * <br>
 * (Rene Gassmoeller, 2014/12/03)
 *
 * <li> New: There is now the possibility to interpolate the visualization
 * output to a refined output mesh. This accounts for the fact that most
 * visualization software only offers linear interpolation between grid points
 * and therefore the output file is a very coarse representation of the actual
 * solution field. The new output is still a linear interpolation, but at
 * least on a finer grid. The visualization output quality and used disk space
 * will increase when activating this option.
 * <br>
 * (Rene Gassmoeller, 2014/09/14)
 *
 * <li> Fixed: There was a race condition in writing output files: we write to
 * them in a temporary location and then move them to their final location. We
 * used to do the second step before actually closing the file we just wrote
 * to. Under some circumstances, this could lead to incomplete or empty output
 * files. This is now fixed.
 * <br>
 * (Wolfgang Bangerth, 2014/09/03)
 *
 * <li> Fixed: Running in release mode and with user-defined libraries loaded
 * that were compiled against the debug version of deal.II, or the other way
 * around, likely produces effects that are undesirable. ASPECT now tries to
 * detect these situations and abort the program.
 * <br>
 * (Wolfgang Bangerth, 2014/09/02)
 *
 * <li> New: There is now a new section in the manual that documents running a
 * free surface computation with a crust as a stagnant lid overlying a
 * convecting mantle.
 * <br>
 * (William Durkin IV, Wolfgang Bangerth, 2014/09/01)
 *
 * <li> New: There is now a new section in the manual that documents running
 * the benchmarks proposed by Davies et al.
 * <br>
 * (William Durkin IV, Wolfgang Bangerth, 2014/08/29)
 *
 * <li> Fixed: The initial condition 'adiabatic' had a minor bug in the case
 * of using a model without adiabatic heating and a prescribed bottom boundary
 * layer. In this case the model should create a constant temperature profile
 * with boundary layers to the boundary temperature. Due to a bug the bottom
 * boundary layer amplitude was calculated against an adiabatic profile
 * instead of a constant temperature. This is fixed now.
 * <br>
 * (Rene Gassmoeller, 2014/08/19)
 * <li> Changed: The GPlates velocity boundary plugin prescribed only the
 * normal vector of the model plane of 2D models. Therefore it was very hard
 * to rotate the model into its actual orientation, when comparing it with
 * other datasets. Now the first of the user input points is always rotated to
 * a known position (0,1,0), and the plugin outputs the according rotation
 * angles to rotate the model into its proper orientation and the inverse
 * angles to rotate other datasets into the model plane.
 * <br>
 * (Rene Gassmoeller, 2014/08/15)
 *
 * <li> New: There are numerous places in the input file where one can or has
 * to input boundary indicators. These indicators identify individual parts of
 * the boundary of the domain, such as the inner or outer surfaces of a
 * spherical shell. However, remembering these numerical indicators when
 * writing or reading input files is a hassle and error prone.
 * <br>
 * Given this problem, there is now functionality in Aspect so that geometry
 * models can (and, for the existing models, do) provide symbolic names for
 * each of the boundary parts. These symbolic names can then be used in the
 * input files wherever it was previously necessary to use numbers.
 * <br>
 * (Wolfgang Bangerth, 2014/07/21)
 *
 * <li> New: The <code>spherical hexagonal perturbation</code> initial
 * temperature model has been generalized to allow other modes than just the
 * hexagonal one.
 * <br>
 * (Joey Durkin, 2014/08/04)
 *
 * <li> New: The <code>spherical shell</code> geometry model in 2d now takes
 * additional parameters indicating the number of circumferential cells.
 * <br>
 * (Joey Durkin, 2014/08/04)
 *
 * <li> Changed: Updated maximum refinement function with the same bugfix as
 * the minimum refinement function. It now always declares dim number of
 * variables. When using the 'depth' coordinate system, only the first
 * variable is used. The others are set to zero.
 * <br>
 * (Jonathan Perry-Houts, 2014/08/04)
 *
 * <li> Fixed: The viscosity mesh refinement criterion did not ask the
 * material model for the viscosity, therefore it did not work with many
 * material models. Some material models calculate the viscosity anyway, so it
 * worked in that cases. Now the criterion asks correctly for the viscosity.
 * <br>
 * (Rene Gassmoeller, Juliane Dannberg 2014/08/04)
 *
 * <li> Changed: The minimum refinement function now always declares three
 * variables. In the 'depth' option only the first one is set to the depth,
 * the other variables are set to zero.
 * <br>
 * (Rene Gassmoeller, 2014/08/04)
 *
 * <li> New: There are now function parser plugins for the gravity model and
 * boundary temperature model. Gravity can be prescribed in dependence of
 * position and boundary temperature in dependence of position and time.
 * <br>
 * (Sanja Panovska, Rene Gassmoeller, 2014/08/01)
 *
 * <li> New: There is now a postprocessor that outputs multiple material
 * properties with just a single call to the material model. This is more
 * efficient, but only matters for complex material models.
 * <br>
 * (Rene Gassmoeller, 2014/08/01)
 *
 * <li> New: There is now a section in the manual detailing the "Burstedde"
 * benchmark and its verification.
 * <br>
 * (Iris van Zelst, Wolfgang Bangerth, 2014/07/07)
 *
 * <li> New: Added "maximum refinement function" plugin to limit the depth of
 * the mesh octree.
 * <br>
 * (Jonathan Perry-Houts, 2014/07/03)
 *
 * <li> Changed: The viscosity in the Steinberger material model can now be
 * calculated by taking either the laterally averaged temperature or the
 * adiabatic temperature at this depth as reference. (Rene Gassmoeller,
 * 2014/06/20)
 *
 * <li> Changed: The minimum refinement function is now evaluated at every
 * cell corner instead of the cell midpoint. This especially helps in refining
 * thin boundary layers.
 * <br>
 * (Rene Gassmoeller, 2014/06/20)
 *
 * <li> Changed: The internal heating statistics postprocessor now calculates
 * the integral and averaged internal heating by mass averaging instead of
 * volume averaging.
 * <br>
 * (Rene Gassmoeller, 2014/06/20)
 *
 * <li> New: Added a visualization processor to output the local gravity
 * vector.
 * <br>
 * (Timo Heister, 2014/06/17)
 *
 * <li> Fixed: The direct solver now works with nonlinear schemes like
 * iterated IMPES.
 * <br>
 * (Timo Heister, 2014/06/11)
 *
 * <li> Changed: The maximum number of solver iterations for all linear
 * solvers (including the inner solvers) is now set to 1000.
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
 * every processor would output the same error message -- leading to
 * incredible amounts of entangled error messages. This has now been resolved:
 * every processor still fails, but only processor 0 reports the error.
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
