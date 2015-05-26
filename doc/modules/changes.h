/**
 * @page changes_current Changes after the latest release (v1.3)
 *
 * <p> This is the list of changes made after the release of Aspect version
 * 1.3. All entries are signed with the names of the author. </p>
 *
 *
 * <ol>
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
 * <li> Changed: The specific heating plugin has a new interface now; it gets
 * the material model inputs and outputs and fills a vector with heating
 * model outputs for the whole cell.
 * <br>
 * (Juliane Dannberg, 2015/05/20)
 *
 * <li> New: There is now a mesh refinement criterion called Slope, intended
 * for use with a free surface, which refines where the topography has a 
 * large slope.
 * <br>
 * (Ian Rose, 2015/05/21)
 *
 * </ol>
 */
