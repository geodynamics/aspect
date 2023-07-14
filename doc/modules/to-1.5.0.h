/**
 * @page changes_between_1.4.0_and_1.5.0 Changes between version 1.4.0 and version 1.5.0
 *
 * <p> This is the list of changes made after the release of ASPECT version
 * 1.4.0 for version 1.5.0. All entries are signed with the names of the author.
 * </p>
 *
 * <ol>
 *
 * <li> New: Blankenbach benchmarks were added.
 * <br>
 * (Timo Heister, 2017/02/16)
 *
 * <li> New: The particle property interpolation interface can now be
 * queried for only a part (or only one) of the particle properties.
 * This makes interpolation more efficient, because the compositional
 * fields are solved one-by-one and therefore oftentimes only one
 * of the properties is actually needed. Interpolating all needed
 * properties in one call is still more efficient in the general case.
 * <br>
 * (Rene Gassmoeller, 2017/02/15)
 *
 * <li> Removed: The material model interface contained the functions
 * 'reference_density', 'reference_thermal_expansion_coefficient',
 * 'viscosity_ratio', and 'thermodynamic_phase'. Additionally many
 * material models implemented functions like 'reference_cp', and
 * 'reference_thermal_diffusivity'. These functions were not used in
 * any place and could be inconsistent with the model properties,
 * because they had no input to figure out the real model conditions.
 * Therefore they were removed. Additionally, the postprocessors
 * visualizing the 'viscosity_ratio' and 'thermodynamic_phase' functions
 * were removed, because there is no material model in the current version
 * that actually implements these functions, and they are not tested.
 * The new recommended way to determine a reference
 * material property for a given model is to evaluate the material model
 * at reference conditions (given by the adiabatic profile) as implemented
 * in the 'basic statistics' postprocessor. Existing user plugins that
 * implement these functions will continue to work, as long as the
 * functions are only used within that plugin. User plugins that called these
 * interface functions from outside can either cast the material model
 * reference to a particular material model that implements these functions,
 * or use the recommended approach to determine reference properties described
 * above.
 * <br>
 * (Rene Gassmoeller, 2017/02/09)
 *
 * <li> Changed: The 'basic statistics' postprocessor now uses a better
 * defined reference state to compute properties like the Rayleigh number.
 * In particular, it computes material properties for the prescribed
 * adiabatic temperature and pressure conditions at the surface, instead
 * of using the deprecated reference_property() functions of the material
 * model interface. This way it can also be used for other material models
 * than the 'simple' model, and also highlights possible inconsistencies
 * between the reference profile and reference temperatures assumed in the
 * matherial model. During the rework the multiplicative output while doing
 * initial adaptive refinement steps was also removed.
 * <br>
 * (Rene Gassmoeller, 2017/02/08)
 *
 * <li> Removed: The 'steinberger' material model contained an option to
 * use it for incompressible models. This option included a somewhat
 * complicated density scaling, and relied on the availability of a
 * compressible adiabatic profile in incompressible models. This profile
 * is not longer available with the recent changes to the adiabatic
 * reference profile, and the feature is too obscure to justify a fix.
 * Therefore, it was removed and the 'steinberger' material model now
 * only works for compressible formulations.
 * <br>
 * (Rene Gassmoeller, 2017/02/07)
 *
 * <li> New: A cookbook for continental extension was added, which uses
 * uses the visco plastic material to simulate the thermal-mechanical
 * evolution of a continental plate extending at a constant velocity.
 * <br>
 * (John Naliboff, 2017/01/20)
 *
 * <li> Changed: The finite strain cookbook used a way to compute
 * the finite strain that was only appropriate for a small amount of
 * strain. This was fixed to a proper finite strain computation.
 * Added tests and benchmarks for the new formulation, and updated
 * the according particle property in the same way. Also updated
 * and extended the cookbook description in the manual.
 * <br>
 * (Rene Gassmoeller, 2017/01/19)
 *
 * <li> New: Added a "heat flux densities" postprocessor.
 * <br>
 * (Timo Heister, 2016/12/26)
 *
 * <li> New: ASPECT now allows querying its own version number as well as
 * the version of all of the underlying libraries using the
 * <code>--version</code> or <code>-v</code> command line flags.
 * <br>
 * Similarly, using <code>--help</code> or <code>-h</code> allows
 * querying command line usage of ASPECT.
 * <br>
 * (Wolfgang Bangerth, 2016/12/19-31)
 *
 * <li> New: A material model for incompressible (using the Boussinesq
 * approximation) and compressible computations (with ALA or TALA) for a
 * nondimensionalized problem. This can be used for several benchmark problems
 * like Blankenbach, King, etc..
 * <br>
 * (Timo Heister, Juliane Dannberg, Rene Gassmoeller, 2016/12/18)
 *
 * <li> New: ASPECT now supports the choice between different formulations for
 * the governing equations including Boussinesq and anelastic liquid
 * approximation. For this, the adiabatic conditions have been extended to
 * provide values and gradients of the reference density. Several benchmarks
 * for these formulations have been added.
 * <br>
 * (Juliane Dannberg, Rene Gassmoeller, Timo Heister, 2016/12/14)
 *
 * <li> New: ASPECT now also generates a <code>output/particles.visit</code>
 * file that allows Visit to read in all files from all processors and
 * all time steps that contain particle data.
 * <br>
 * (Wolfgang Bangerth, 2016/12/02)
 *
 * <li> Changed: The 'Stokes only' nonlinear solver scheme now uses
 * the nonlinear solver tolerance parameter (previously it was hard-
 * coded to 1e-8), and it computes the nonlinear residual as current
 * residual divided by initial residual, consistent with the 'iterated
 * Stokes' solver scheme.
 * <br>
 * (Juliane Dannberg, Rene Gassmoeller, 2016/11/22)
 *
 * <li> Changed: The adiabatic profile now contains a reference density
 * profile, and the derivative of this reference density profile. The
 * InitialProfile adiabatic profile now relies on the adiabatic heating
 * to decide if the temperature increases with depth. Additionally, an
 * off-by-one bug was fixed in InitialProfile leading to minor changes
 * in the adiabatic profiles (relative change of ~1e-4), changing models
 * that rely on the profile. Models that rely on the adiabatic profile
 * might be changed by this PR if they are compressible, but do not
 * use adiabatic heating, or vice versa.
 * <br>
 * (Timo Heister, Juliane Dannberg, Rene Gassmoeller, 2016/11/20)
 *
 * <li> New: The visco plastic material model now includes an option for
 * strain-weakening of cohesion and the internal angle of friction.
 * Strain-weakeing of these properties is commonly used to help localize
 * deformation in regions undergoing plastic failure. Note that using this
 * option requires tracking the finite strain tensor through particles or
 * compositional fields.
 * <br>
 * (John Naliboff, 2016/11/10)
 *
 * <li> New: Two particle generators were added. One, generates particles
 * at the quadrature points for each active cell in the triangulation.
 * Two, generates a uniform distribution of particles in the unit cell
 * and transforms each of the particles back to real region in the model
 * domain for each active cell in the triangulation.
 * <br>
 * (Harsha Lokavarapu, 2016/11/10)
 *
 * <li> Changed: The exchange of ghost particles that was introduced lately
 * can be quite expensive for models with many particles,
 * and is often unnecessary if the particles are used as passive tracers.
 * Therefore, a new input parameter 'Update ghost particles' controls this
 * exchange, and its default is set to 'false'. Model parameter files using
 * active particles will need to be changed accordingly.
 * <br>
 * (Rene Gassmoeller, 2016/10/18)
 *
 * <li> Improved: The matrix assembly of Stokes and Advection systems has been
 * optimized, by assembling less (only the relevant) DoFs, and by optimizing
 * calls to deal.II functions. The overall speedup for box models is between
 * 20 and 40% of the assembly time, likely somewhat less for curved geometries.
 * This change will require changes in user written assembler plugins, because
 * the Stokes system assembly now only loops over Stokes degrees of freedom.
 * <br>
 * (Rene Gassmoeller, 2016/10/17)
 *
 * <li> Improved: Box models without deformed mesh now use a MappingCartesian,
 * which assumes all mesh cells are aligned with cartesian coordinate axes.
 * Matrix assembly and particle transport in such mappings is around 20 % faster
 * compared to a general MappingQ1 for other box models.
 * <br>
 * (Rene Gassmoeller, 2016/10/14)
 *
 * <li> Changed: HDF5 particle output files are now named 'particles-...'
 * instead of 'particle-...' to be consistent with the vtu output. Also
 * particle properties with more than one component are now correctly split
 * into scalar fields in the output files, if they have more or less components
 * than the number of spatial dimensions in the model.
 * <br>
 * (Rene Gassmoeller, 2016/09/20)
 *
 * </li>
 * <li> New: Multiple particle properties can be initialized by specifying
 * multiple particle property function components as opposed to one particle
 * property.
 * <br>
 * (Harsha Lokavarapu, Gerry Puckett, 2016/09/20)
 *
 * <li> Changed: The timestep entry in the statistics file has been moved to
 * column 3 and is now the timestep used for the timestep corresponding to the
 * current row.
 * <br>
 * (Jonathan Robey, 2016/09/16)
 * </li>
 *
 * <li> Changed: The 'cell average' particle interpolator is now more
 * tolerant against cells without particles by interpolating properties
 * from neighboring cells. This is necessary, because during refinement
 * even children of cells with a reasonable number of particles can be
 * void of particles.
 * <br>
 * (Rene Gassmoeller, Jonathan Perry-Houts, 2016/08/31)
 *
 * <li> Changed: Particle properties should now declare which solution
 * properties they need to update themselves. The particle world then
 * only computes values and gradients of the solution at
 * the particle positions if necessary, which can reduce the computational
 * cost of the particle update for simple particle properties.
 * <br>
 * (Rene Gassmoeller, 2016/08/30)
 *
 * <li> New: .visit output files now also contain information about
 * the model time, as long as ASPECT was build with at least
 * deal.II 8.5.0.pre. Previously, this information was only available
 * in the Paraview .pvd files.
 * <br>
 * (Rene Gassmoeller, Juliane Dannberg, 2016/08/24)
 *
 * <li> New: There is now an initial topography plugin that returns
 * initial topography values based on an ascii data file.
 * <br>
 * (Anne Glerum, 2016/08/22)
 *
 * <li> Fixed: The point value postprocessor forgot to take into
 * account the mapping we use when describing curved boundaries.
 * <br>
 * (Rene Gassmoeller, Wolfgang Bangerth, 2016/08/16)
 *
 * <li> Changed: Particles now also store their location in the
 * coordinate system of their current cell. This decreases the
 * number of times this location has to be computed by inverting
 * the mapping for the current cell, which is expensive.
 * On average this change will save 40-50% of the overall
 * particle computing time, while increasing the particle
 * memory footprint (which is usually small compared to the
 * system matrix).
 * <br>
 * (Rene Gassmoeller, 2016/08/12)
 *
 * <li> Fixed: The chunk geometry pull back function now returns
 * a corrected longitude value when 180 hemisphere is crossed.
 * <br>
 * (Anne Glerum, 2016/08/09)
 *
 * <li> Changed: It is now possible to read in ascii data files of
 * which the coordinates are not equally spaced.
 * <br>
 * (Anne Glerum, 2016/08/05)
 *
 * <li> New: It is now possible to create compositional fields that are
 * not advected by a field method, but interpolated from particle properties.
 * This simplifies the process of using the particles as 'active' particles
 * that carry information that influences the solution. E.g. the material
 * model can access the compositional field that is interpolated from particle
 * properties and use this as 'composition' information.
 * <br>
 * (Rene Gassmoeller, 2016/08/02)
 *
 * <li> New: There is now an initial topography plugin which reads
 * from the prm file polygon definitions and set the initial topography
 * to be constant within those polygons.
 * <br>
 * (Menno Fraters, 2016/07/26)
 *
 * <li> Changed: Particle initialization no longer routinely computes
 * the solution at the particle positions, since it is usually not needed
 * and complicates the initialization process. Instead it evaluates the
 * initial conditions at the particle positions. It is still possible to
 * access the solution by evaluating it manually inside of particle
 * property plugins. Additionally the 'initial composition' property
 * now utilizes the user-provided names of the compositional fields
 * to identify its particle properties (they are now named
 * 'initial [field_name]', where '[field_name]'
 * is replaced by the user provided name).
 * <br>
 * (Rene Gassmoeller, 2016/07/18)
 *
 * <li> New: Added parameter “Adapt by fraction of cells” to switch between
 * refining a certain number of cells based on the fraction of total error
 * (default behavior) and the fraction of total number of cells
 * <br>
 * (Lev Karatun, 2016/07/20)
 *
 * <li> Changed: It is now possible to set the gravity to a negative
 * value in order to calculate backward advection.
 * <br>
 * (Jacky Austermann, 2016/07/13)
 *
 * <li> New: There is a new postprocessor that outputs statistics to
 * the screen about the memory usage and nonzero entries of matrices.
 * It can be called with the name 'matrix statistics'.
 * <br>
 * (Sam Cox, 2016/07/01)
 *
 * <li> New: There is now a plugin structure to add initial topography
 * to geometry models.
 * <br>
 * (Menno Fraters and Anne Glerum, 2016/07/01)
 *
 * <li> New: There is a new postprocessor that outputs statistics about
 * the memory usage at each timestep. It can be called with the name
 * 'memory statistics'. This replaces the helper function
 * 'output_program_stats()', which has been removed, along with the
 * variable 'output_parallel_statistics'.
 * <br>
 * (Sam Cox, 2016/06/30)
 *
 * <li> Changed: In input files, lines that end in a backslash now require
 * a space before the backslash if the two lines should not be conjoined
 * immediately, because leading spaces on the continuing line is ignored.
 * See the Section "The structure of parameter files" in the manual.
 * <br>
 * (Jonathan Robey, 2016/06/30)
 *
 * <li> Changed: The default for "Initial adaptive refinement" cycles is
 * now 0.
 * <br>
 * (Timo Heister, 2016/06/28)
 *
 * <li> New: There is a new parameter gamma in the entropy viscosity
 * stabilization that allows to scale the stabilization with the strain rate
 * (in addition to the velocity).
 * <br>
 * (Juliane Dannberg, 2016/06/28)
 *
 * <li> New: There is now a postprocessor that outputs the heatflux
 * density at each boundary face into a text file and a
 * postprocessor that outputs the heatflux density at each boundary
 * face for visualization.
 * <br>
 * (Jacky Austermann, 2016/06/28)
 *
 * <li> New: There is a new function Utilities::real_spherical_harmonic
 * which calculates the values of a fully normalized real spherical harmonic.
 * <br>
 * (Ian Rose, 2016/06/28)
 *
 * <li> Changed: The files handling the free surface implementation
 * have been renamed to free_surface.h and free_surface.cc for better
 * consistency with the rest of the file names.
 * <br>
 * (Ian Rose, 2016/06/27)
 *
 * <li> Changed: The previously-deprecated functions 'composition' and
 * 'temperature' have now been removed. Their replacements are
 * 'boundary_composition' and 'boundary_temperature'.
 * <br>
 * (Sam Cox, 2016/06/27)
 *
 * <li> Changed: The free surface implementation now represents mesh
 * deformation using a vector of displacements for the mesh vertices,
 * using a MappingQ1Eulerian to transform from the reference cell to
 * the physical cell.
 * <br>
 * (Ian Rose, 2016/06/26)
 *
 * <li> Changed: The visualization postprocessor now writes every initial
 * refinement stage if 'Run postprocessors on initial refinement' is
 * set to true.
 * <br>
 * (Rene Gassmoeller, 2016/06/25)
 *
 * <li> New: There is now a postprocessor that outputs statistics about
 * the number of particles per cell in the model domain. It outputs
 * global maximum, minimum and average number of particles per cell.
 * <br>
 * (Rene Gassmoeller, 2016/06/24)
 *
 * <li> New: There is now the option to model melt transport (two-phase
 * flow). This core of the implementations includes additional
 * variables in the solution vector, a new assembler with an additional
 * equation that will be solved and a modified advection equation for the
 * porosity field, a new preconditioner for models with melt transport, and
 * additional melt outputs for the material model.
 * <br>
 * (Juliane Dannberg, Timo Heister, 2016/06/24)
 *
 * <li> New: Particles can now carry the integrated strain they have
 * experienced over the course of the model. They store all components
 * of the symmetric strain tensor, which can be converted into the
 * invariants or investigated individually.
 * <br>
 * (Rene Gassmoeller, 2016/06/10)
 *
 * <li> New: There is a new optional feature for the discontinuous temperature
 * and compositional solutions. After solving the advection equation,
 * a "bound preserving limiter" working as a correction procedure is applied
 * to the discontinuous advection fields. The limiter will stabilize the
 * discontinuous advection solutions and keep it in the range of user defined
 * global maximum/minimum values. Whether or not the limiter is used is
 * determined by an entry to the parameter file.
 * <br>
 * (Ying He, 2016/06/02)
 *
 * <li> New: Tests can now be marked that they are expected to fail by the
 * EXPECT FAILURE keyword in the .prm.
 * <br>
 * (Timo Heister, 2016/06/01)
 *
 * <li> New: There is a new boundary traction model "ascii data"
 * that prescribes the boundary traction according to pressure values
 * read from an ascii data file.
 * <br>
 * (Juliane Dannberg, 2016/05/24)
 *
 * <li> New: A material model plugin for visco-plastic rheologies,
 * which combines diffusion, dislocation or composite viscous
 * creep with a Drucker Prager yield criterion.
 * <br>
 * (John Naliboff, 2016/05/19)
 *
 * <li> Changed: The traction boundary conditions now use a new interface
 * that includes the boundary indicator (analogous to the velocity boundary
 * conditions).
 * <br>
 * (Juliane Dannberg, 2016/05/19)
 *
 * <li> New: There is now a visualization plugin to visualize the maximum
 * horizontal component of the compressive stress.
 * <br>
 * (Wolfgang Bangerth, D. Sarah Stamps, 2016/05/12)
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
 * </ol>
 */
