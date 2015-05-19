/**
 * @page changes_current Changes after the latest release (v1.3)
 *
 * <p> This is the list of changes made after the release of Aspect version
 * 1.3. All entries are signed with the names of the author. </p>
 *
 *
 * <ol>
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
 * </ol>
 */
