/**
 * @page changes_current Changes after the latest release (v1.3)
 *
 * <p> This is the list of changes made after the release of Aspect version
 * 1.3. All entries are signed with the names of the author. </p>
 *
 * <ol>
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
 * <br>
 * (Scott Tarlow, 2015/05/29)
 *
 * <li> New: Three new material averaging schemes are added. These are combined
 * with the existing averaging schemes, except for the Q1 averaging schemes into
 * a compositing material model. The new averaging schemes are normalised weighted
 * distance versions of the arithmetic, harmonic and geometric averages.
 * <br>
 * (Menno Fraters, 2015/05/28)
 *
 * <li> New: There are now postprocessors BoundaryDensities and
 * BoundaryPressures which calculate laterally averaged densities
 * and pressures at the top and bottom of the domain.
 * <br>
 * (Ian Rose, 2015/05/28)
 *
 * <li> New: Postprocessor and visualization postprocessor plugins can
 * now state that they require other postprocessors to run as well,
 * for example because they want to query information computed
 * by these other postprocessors in computing their own information.
 * This is done using the
 * aspect::Postprocess::Interface::requires_other_postprocessors()
 * function.
 * <br>
 * (Wolfgang Bangerth, 2015/05/28)
 *
 * <li> New: Added cookbook to prescribe initial condition from shear
 * wave velocity model.
 * <br>
 * (Jacqueline Austermann, 2015/05/28)
 *
 * <li> Changed: The heating models have a new straucture now:
 * Instead of the implementation in the assembly, there is a heating
 * plugin for each model now that can be used both in the assembly
 * and the postprocessors, and a new heating model manager that
 * combines the plugins by adding the individual heating terms.
 * <br>
 * (Juliane Dannberg, 2015/05/27)
 *
 * <li> Fixed: There was a bug in the make pressure rhs compatibility
 * function that caused the linear solver to fail in models with a
 * significant in- or outflux of material. This is fixed now.
 * <br>
 * (Juliane Dannberg, 2015/05/27)
 *
 * <li> New: Added cookbook for prescribed internal velocity values.
 * <br>
 * (Jonathan Perry-Houts, 2015/05/26)
 *
 * <li> New: There is now a Material model called diffusion dislocation that
 * implements diffusion and dislocation creep viscosity.
 * <br>
 * (Bob Myhill, 2015/05/26)
 *
 * <li> Changed: Modified S40RTS initial condition file to incorporate
 * the option to zero out heterogeneities within a given depth.
 * <br>
 * (Jacqueline Austermann, 2015/05/26)
 *
 * <li> New: There is now a refinement plugin based on strain rate which
 * will come handy to capture shear bands when combined with plasticity.
 * <br>
 * (Cedric Thieulot, 2015/05/26)
 *
 * <li> New: There is now a Material model called Depth dependent that implements
 * depth-dependent viscosity through the modification of any other Material model.
 * <br>
 * (Jacqueline Austermann and Max Rudolph, 2015/05/24)
 *
 * <li> The manual.pdf is no longer part of the git repository but you can
 * find it online at http://aspect.dealii.org or you can build it yourself.
 * <br>
 * (Timo Heister, 2015/05/23)
 *
 * <li> New: There are now a set of global constants defined for physical
 * properties and for radius and gravity parameters relevant to
 * Earth and Mars.
 * <br>
 * (Bob Myhill, 2015/05/22)
 *
 * <li> New: There is now a SimulatorAccess::get_statistics_object() function
 * that allows all users of this class to record information in the statistics
 * file, whether they are postprocessors or not.
 * <br>
 * (Wolfgang Bangerth, 2015/05/22)
 *
 * <li> New: ASPECT now also provides a signals mechanism to attach user-defined
 * functions to certain events over the course of a simulation. This allows more
 * fine-grained observation and intervention to user additions. A new section
 * in the manual explains how to extend ASPECT this way.
 * <br>
 * (Wolfgang Bangerth, 2015/05/21)
 *
 * <li> Changed: The free surface handler now detaches internal manifolds
 * for cases when the domain has them, since they are not necessarily a
 * good description of the geometry when there has been large mesh deformation.
 * <br>
 * (Ian Rose, 2015/05/21)
 *
 * <li> Changed: The documentation for nullspace removal is now more
 * descriptive of what Aspect is actually doing.
 * <br>
 * (Ian Rose, 2015/05/21)
 *
 * <li> New: There is now a mesh refinement criterion called Slope, intended
 * for use with a free surface, which refines where the topography has a
 * large slope.
 * <br>
 * (Ian Rose, 2015/05/21)
 *
 * <li> Changed: The specific heating plugin has a new interface now; it gets
 * the material model inputs and outputs and fills a vector with heating
 * model outputs for the whole cell.
 * <br>
 * (Juliane Dannberg, 2015/05/20)
 *
 * <li> New: There is now an (3D) ellipsoidal chunk geometry model where two
 * of the axis have the same length. The ellipsoidal chunk can be non-coordinate
 * parallel part of the ellipsoid.
 * This plugin is a joined effort of Menno Fraters, D Sarah Stamps and Wolfgang
 * Bangerth
 * <br>
 * (Menno Fraters, 2015/08/28)
 * </ol>
 */
