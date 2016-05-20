/**
 * @page changes_current Changes after the latest release (v1.4.0)
 *
 * <p> This is the list of changes made after the release of Aspect version
<<<<<<< be6cd0b0e203518b84f6c5f12ec914becece0b1e
 * 1.4.0. All entries are signed with the names of the author. </p>
 * 
 * <li> Changed: The traction boundary conditions now use a new interface
 * that includes the boundary indicator (analogous to the velocity boundary 
 * conditions).
=======
 * 1.3. All entries are signed with the names of the author. </p>
 *
 * <ol>
 *
 * <li> New: A material model plugin for visco-plastic rheologies,
 * which combines diffusion, dislocation or composite viscous 
 * creep with a Drucker Prager yield criterion.
 * <br>
 * (John Naliboff, 2016/05/19)
 *
 * <li> New: The variable $ASPECT_SOURCE_DIR can now be used inside
 * .prm files in various locations (including "include" statements).
 * <br>
 * (Timo Heister, 2016/05/02)
 *
 * <li> Improved: Mesh refinement plugins can now update member variables
 * at the beginning of each time step. This can be used to make things like
 * time-dependent mesh refinement criteria, or fancy models that depend on
 * calling an external program. The most immediate added functionality is
 * the ability to define time-dependent minimum/maximum refinement functions.
 * <br>
 * (Jonathan Perry-Houts, 2016/04/13)
 *
 * <li> New: There is a new function Utilities::orthogonal_vectors() that
 * allows computing one (in 2d) or two (in 3d) orthogonal vectors to a given
 * input vector.
 * <br>
 * (Wolfgang Bangerth, 2016/04/13)
 *
 * <li> Improved: The particle sorting algorithm performance was considerably
 * improved. It now sorts neighbor cells by the distance between the particle
 * and the face center that is between the old cell and the neighbor and then
 * checks if the particle is in this list of neighbor cells. Only if it is
 * not found there, a search over all local cells is performed.
 * <br>
 * (Rene Gassmoeller, Wolfgang Bangerth, 2016/04/12)
 *
 * <li> New: There is now a postprocessor "point values" that allows evaluating
 * the solution at a number of user-defined evaluation points.
 * <br>
 * (Wolfgang Bangerth, 2016/03/24)
 *
 * <li> Fixed: Combining particles with initial adaptive refinement steps
 * used to create a multiple of the selected number of particles.
 * This is fixed now.
 * <br>
 * (Rene Gassmoeller, 2016/03/18)
 *
 * <li> Improved: The Introspection class has now a new base class
 * FEVariableCollection that allows flexible modification of the finite
 * element variables involved in a computation, even inside a plugin.
 * <br>
 * (Timo Heister, 2016/03/08)
 *
 * <li> Fixed: The uniform radial and uniform box particle generators now
 * produce globally unique particle IDs.
 * <br>
 * (Harsha Lokavarapu, Gerry Puckett, 2016/03/04)
 *
 * <li> Fixed: The 'Simpler' material model produced floating point exceptions
 * in models with compositional fields. This is fixed now.
 * <br>
 * (Lev Karatun, Rene Gassmoeller, 2016/02/26)
 *
 * <li> New: The advection systems (for temperature and compositions) can now
 * be discretized using the symmetric interior penalty discontinuous Galerkin
 * method. This can be useful to explore solution without adding artificial
 * smoothing. This is controlled by two new input parameters in
 * 'Discretization': use_discontinuous_temperature_discretization and
 * use_discontinuous_composition_discretization.
 * <br>
 * (Sam Cox, 2016/02/22)
 *
 * <li> Changed: ASPECT by default wrote one output file per MPI process that
 * was written in a background thread to a temporary location first and then
 * moved to its final location. If any of the steps failed it tried again by
 * writing directly to the output location. This approach needed complicated
 * logic and did not succeed on all systems. In order to increase stability
 * the new default behaviour is to write straight to the output folder. This
 * might decrease performance on clusters with slow network file systems.
 * The old behaviour can be recovered by setting 'Write in background thread'
 * to true and set a temporary storage location by 'set Temporary output
 * location'. Note that this functionality was and is only available if
 * 'Number of grouped files' is set to its default value of 0, and therefore
 * MPI-IO is not used for parallel output. For larger models with hundreds of
 * parallel processes using MPI-IO is recommended.
 * <br>
 * (Rene Gassmoeller, 2016/02/14)
 *
 * <li> New: Added 'command' postprocessor for executing arbitrary commands.
 * <br>
 * (Jonathan Perry-Houts, 2016/02/11)
 *
 * <li> Improved: The option to increase the output resolution by linear
 * interpolation of the quadratic elements now correctly uses the mapping of
 * curved geometries to interpolate cells. This increases output accuracy for
 * models that use curved geometries and use 'Set Interpolate output = true'.
 * The simulation itself is not affected.
 * <br>
 * (Rene Gassmoeller, 2016/02/08)
 *
 * <li> Changed: The GPlates plugin is restructured in the style of the
 * AsciiData Plugin. The major difference is that the interpolation is now
 * performed in spherical coordinates instead of Cartesian coordinates. Note
 * that some input parameters have changed: "Time step" is now called "Data
 * file time step", "Velocity file start time" is now called "First data file
 * model time", "Interpolation width" does not exist any more, but there are
 * three new parameters called "First data file number", "Decreasing file
 * order" and "Lithosphere thickness".
 * <br>
 * (Eva Bredow, Rene Gassmoeller, 2016/02/04)
 *
 * <li> New: ASPECT no longer relies on the availability of a command-
 * processor (terminal) at run-time, by providing fallbacks to C
 * commands. This adds support for architectures that do not offer a
 * terminal on compute nodes (like IBM BlueGene/Q).
 * <br>
 * (Rene Gassmoeller, 2016/02/03)
 *
 * <li> Changed: The 'depth' function of the 'box' geometry model and
 * the 'two merged boxes' geometry model previously threw an exception
 * when asked for the depth of a point outside of the initial model domain.
 * This is not longer appropriate for models with free surfaces and therefore
 * the behaviour was changed to the behaviour of the 'spherical shell' geometry
 * model, which is a cutoff of the depth to the range (0,maximal_depth).
 * <br>
 * (Rene Gassmoeller, Sascha Brune, 2016/01/11)
 *
 * <li> New: There is now a parameter called 'Additional tangential
 * mesh velocity boundary indicators' that allows to specify boundaries
 * which elements are allowed to deform tangential to the boundary.
 * This can be useful in models with free surface and a prescribed
 * material in-/outflow at the sides. Previously in this case the
 * uppermost element became distorted over time, now the whole
 * boundary mesh adjusts according to the deformation. This change
 * also fixes the handling of traction boundary conditions in models
 * with free surface.
 * <br>
 * (Anne Glerum, Rene Gassmoeller, Ian Rose, 2016/01/11)
 *
 * <li> Changed: The interfaces of the boundary composition and boundary
 * temperature plugins have been deprecated. Their replacements not longer
 * contain references to the geometry model, which was a leftover from an
 * earlier development stage. Users should derive their plugins from
 * SimulatorAccess if they need access to the geometry model. The
 * deprecated functions will be removed in a future ASPECT release.
 * <br>
 * (Rene Gassmoeller, 2016/01/04)
 *
 * <li> New: A new mesh refinement plugin was added that refines cells
 * according to the density of particles in that cell.
 * <br>
 * (Rene Gassmoeller, 2015/12/19)
 *
 * <li> Changed: The boundary_velocity(const Point<dim> &position) const
 * function has now been deprecated in favor of the new function
 * boundary_velocity (const types::boundary_id boundary_indicator,
 * const Point<dim> &position) const.
 * <br>
 * (Menno Fraters, 2015/12/16)
 *
 * <li> New: Visualization postprocessors for thermal conductivity and
 * thermal diffusivity.
 * <br>
 * (Anne Glerum, 2015/12/03)
 *
 * <li> New: The tracer architecture has been completely overhauled. It is now
 * more flexible and allows for easier modification. Additionally tracers
 * now carry properties with them, which allows for a variety of new use cases.
 * A number of bugs related to curved cells in spherical models with tracers
 * have been resolved.
 * <br>
 * (Rene Gassmöller, 2015/12/01)
 *
 * <li> Fixed: Whenever the base models used by either the "depth dependent"
 * or "averaging" material models depended on anything that requires accessing
 * the simulator, then this led to segmentation faults. This is now fixed.
 * <br>
 * (Wolfgang Bangerth, Shangxin Liu, 2015/11/24)
 *
 * <li> New: The tolerance of the preconditioners of the A and S block
 * are now available as parameters in the prm file. There is now also
 * a section added to the manual on how to use these parameters to
 * make ASPECT in certain situation faster.
 * <br>
 * (Menno Fraters, 2015/11/08)
 *
 * <li> Changed: The DynamicTopography postprocessor and visualization
 * plugins now use the more accurate Gaussian quadrature rule for evaluating
 * cell averages of the surface stress.
 * <br>
 * (Ian Rose, 2015/11/04)
 *
 * <li> New: Add depth postprocessor which visually outputs the
 * depth for all points inside the domain, as determined by the
 * geometry model.
 * <br>
 * (Menno Fraters, 2015/10/15)
 *
 * <li> New: The "Stokes residual" postprocessor will output the convergence
 * of the Stokes residual of each linear solve.
 * <br>
 * (Timo Heister, 2015/10/21)
 *
 * <li> New: The DepthAverage postprocessor can now output as a plain
 * ascii text table.
 * <br>
 * (Ian Rose, 2015/10/05)
 *
 * <li> Fixed: Free surface computations now work with checkpointing.
 * <br>
 * (Ian Rose, 2015/09/29)
 *
 * <li> Fixed: certain combinations of boundary conditions with a free surface
 * could result in accessing nonexistent matrix entries. This is fixed.
 * <br>
 * (Ian Rose, 2015/09/28)
 *
 * <li> New: Ability to automatically resume models from a checkpoint
 * if previous checkpoint files exist (Resume computation = auto).
 * <br>
 * (Jonathan Perry-Houts, 2015/09/25)
 *
 * <li> Changed: the user can now select a subset of the laterally-
 * averaged quantities to be computed in the DepthAverage postprocessor.
 * <br>
 * (Ian Rose, 2015/09/04)
 *
 * <li> New: The history of the Stokes solver residuals is saved and can
 * be accessed using the post_stokes_solver signal and will be written to
 * a file automatically in case the solver doesn't converge.
 * <br>
 * (Timo Heister, 2015/09/02)
 *
 * <li> Fixed: The laterally averaged sinking velocity and velocity
 * magnitude calculations did not check whether the user selected m/s
 * or m/yr for output values.  Now they do.
 * <br>
 * (Ian Rose, 2015/08/28)
 *
 * <li> Fixed: The lateral averaging of the velocity magnitude was
 * mistakenly calculating the square of the velocity. Now it calculates
 * the magnitude.
 * <br>
 * (Ian Rose, 2015/08/28)
 *
 * <li> Changed: The interface of material models no longer declares
 * property_depends_on() functions. The dependencies of parameters
 * on solution variables are instead handled by a structure in the
 * base class. All included material models have been updated. User
 * written material models will continue to work as before. Since
 * there is no solver yet that utilizes the other dependencies,
 * the impact of this change is limited.
 * <br>
 * (Kimberly Moore, Rene Gassmoeller, Wolfgang Bangerth, 2015/08/26)
 *
 * <li> New: The DepthAverage postprocessor now can calculate the laterally
 * averaged heat flux in the interior of the simulation.
 * <br>
 * (Ian Rose, 2015/08/24)
 *
 * <li> New: There is now a new initial condition in which the temperature field is perturbed
 * following the SAVANI shear wave velocity model by Auer et al., 2014. The data were
 * downloaded from http://n.ethz.ch/~auerl/research.html .
 * <br>
 * (Shangxin Liu, 2015/08/20)
 *
 * <li> New: A box Geometry Model plugin with additional boundary indicators
 * for the upper part of the box and corresponding Boundary Temperature and
 * Composition Model plugins. With this plugin, different boundary conditions
 * can be prescribed on the upper and lower part of the vertical domain boundaries.
 * <br>
 * (Anne Glerum, 2015/08/14)
 *
 * <li> New: Plugin for visualizing the boundary indicators used by the
 * Geometry Model.
 * <br>
 * (Anne Glerum, 2015/08/14)
 *
 * <li> New: There is a new visualization postprocessor which displays
 * the heat flux in the vertical direction, where upwards heat flux
 * is positive.
 * <br>
 * (Ian Rose, 2015/08/12)
 *
 * <li> New: A new material averaging option using logarithms is added.
 * This is combined with the existing averaging schemes. Taking the viscosity for example,
 * the log averaging will average 10^23 and 10^21 to 10^22.
 * <br>
 * (Shangxin Liu, 2015/08/09)
 *
 * <li> New: A material model plugin for Drucker-Prager plasticity.
 * <br>
 * (Anne Glerum, 2015/08/03)
 *
 * <li> Fixed: Quasi-implicit stabilization of a free surface had used the
 * reference density instead of the density as evaluated by the material
 * model.  Now it uses the actual density of at the surface.  This should
 * not change much unless the density at the surface is significantly
 * different from the reference density.
 * <br>
 * (Ian Rose, 2015/07/22)
 *
 * <li> New: Plugin for visualizing groups of compositional fields as vectors.
 * <br>
 * (Jonathan Perry-Houts, 2015/07/12)
 *
 * <li> New: For free surface computations there is an option to advect the
 * mesh vertically (in the direction of gravity),  in addition to the old
 * formulation which advects it in the direction normal to the surface.
 * This can be enabled by setting "Surface velocity projection" to "vertical"
 * in the "Free surface" section of a parameter file.
 * This scheme can maintain better mesh regularity properties for computations
 * where there is a large deformation, or large curvature.
 * <br>
 * (Ian Rose, 2015/07/10)
 *
 * <li> New: There is now an option in the Visualization postprocessor for
 * outputting the mesh velocity in free surface runs.
 * <br>
 * (Ian Rose, 2015/06/16)
 *
 * <li> New: There are now parameter files and a section in the manual for
 * reproducing the benchmarks for free surface computations from Crameri et
 * al. (2012).
 * <br>
 * (Ian Rose, 2015/06/14)
 *
 * <li> New: One can now prescribe the traction on a boundary instead of
 * supplying velocity boundary conditions.
 * This is done in a similar way as for the prescribed velocity boundary conditions:
 * For a given boundary indicator, one can prescribe all or a selection of the
 * traction components.
 * <br>
 * (Wolfgang Bangerth, Anne Glerum, 2015/05/29)
 *
 * <li> Changed: We now use the exact formulation with the
 * compressible strain rate instead of an approximation
 * using the right hand side of the mass conservation equation
 * to calculate the shear heating. This is more accurate in
 * compressible models.
 * <br>
 * (Juliane Dannberg, 2015/05/29)
 *
 * <li> New: There is now a new geometry model called chunk, which
 * takes radius, longitude (and latitude) pairs and creates a regional
 * chunk of a spherical shell. Spherical boundary and initial conditions
 * have been updated to accept this new model. The conversion conventions
 * between [longitude, latitude and radius], [phi, theta and radius] and
 * Cartesian [x, y, z] are consistent with mathematical convention and
 * other models in dealii/ASPECT.
 * This model was based largely on work on deal.ii by D. Sarah Stamps,
 * Wolfgang Bangerth, and in ASPECT by Menno Fraters.
 * <br>
 * (Bob Myhill, 2015/05/29)
 *
 * <li> New: Added the ability to prescribe internal velocities with an ascii
 * file.
>>>>>>> removed statistics and depth average plots for visco_plastic test case
 * <br>
 * (Juliane Dannberg, 2016/05/19)
 *
 * <li> New: There is a new visualization postprocessor "artificial viscosity
 * composition" to visualize the artificial viscosity of a compositional
 * field.
 * <br>
 * (Juliane Dannberg, Timo Heister, 2016/05/03)
 *
 * <li> New: Mesh refinement strategies based on artificial viscosity,
 * composition gradients, or composition threshold.
 * <br>
 * (Juliane Dannberg, Timo Heister, 2016/05/03)
 *
 * <ol>
 *
 * </ol>
 */
