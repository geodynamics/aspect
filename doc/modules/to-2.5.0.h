/**
 * @page changes_between_2.4.0_and_2.5.0 Changes between version 2.4.0 and version 2.5.0
 *
 * <p> This is the list of changes made after the release of ASPECT version
 * 2.4.0 for version 2.5.0. All entries are signed with the names of the author.
 * </p>
 *
 * <ol>
 *
 * <li> Improved: ASPECT now includes version 0.5.0 of the Geodynamic World Builder (https://doi.org/10.5281/zenodo.7998525).
 * <br>
 * (Menno Fraters and other contributors, 2023/06/06)
 *
 * <li> Changed: ASPECT's manual has been converted from LaTeX to Markdown to be
 * hosted as a website on https://aspect-documentation.readthedocs.io.
 * This change was made to make the manual more accessible to users and
 * developers, and to make it easier to update the manual.
 * <br>
 * (Chris Mills, Mack Gregory, Timo Heister, Wolfgang Bangerth, Rene Gassmoeller, and many others, 2023/05/23)
 *
 * <li> New: ASPECT now requires deal.II 9.4 or newer.
 * <br>
 * (Rene Gassmoeller, Timo Heister, 2023/04/29)
 *
 * <li> Improved: Significantly reduced the number of duplicated file system
 * accesses across ASPECT. This should increase the stability and
 * potentially the speed for models with large number of MPI ranks
 * and systems with slow parallel file systems.
 * <br>
 * (Rene Gassmoeller, 2023/04/14)
 *
 * <li> Changed: ASPECT now computes the pressure scaling factor that scales the
 * pressure solution to a similar magnitude as the velocity solution in
 * every Stokes solve instead of once per timestep. This reduces
 * the necessary linear solver iterations in nonlinear iterations
 * after the first if the viscosity changes significantly.
 * The cost of every nonlinear iteration increases
 * slightly because of this change.
 * <br>
 * (Rene Gassmoeller, 2023/03/22)
 *
 * <li> New: ASPECT now has a cmake option ASPECT_INSTALL_EXAMPLES that allows
 * to build and install all cookbooks and benchmarks. This will include
 * building all shared libraries during 'make' and installing all files
 * to the installation location on 'make install'. This is helpful for
 * installations that are used for teaching and tutorials.
 * <br>
 * (Rene Gassmoeller, 2023/03/13)
 *
 * <li> New: Plugins in shared libraries now build libmy_plugin.debug.so
 * or libmy_plugin.release.so depending on the current build type.
 * The "Additional shared libraries" parameter value does not need
 * to be adjusted, because ASPECT automatically loads the correct
 * library when "./libmy_plugin.so" is specified.
 * <br>
 * (Timo Heister, 2023/03/09)
 *
 * <li> New: ASPECT now supports a DebugRelease build type that creates
 * a debug build and a release build of ASPECT at the same time. It
 * can be enabled by setting the CMake option CMAKE_BUILD_TYPE to
 * DebugRelease or by typing "make debugrelease".
 * <br>
 * (Timo Heister, 2023/03/03)
 *
 * <li> New: There is a new initial composition plugin that reads in
 * a slab model (e.g. Slab2) and converts it to a compositional field.
 * <br>
 * (Arushi Saxena, Rene Gassmoeller, Juliane Dannberg, 2023/03/03)
 *
 * <li> New: The signal start_timestep can be used to get a
 * callback at the beginning of each timestep.
 * <br>
 * (Timo Heister, 2023/02/22)
 *
 * <li> Fixed: The depth average postprocessor would crash when trying to compute the
 * quantities "adiabatic_temperature", "adiabatic_pressure", "adiabatic_density",
 * and "adiabatic_density_derivative" for a single depth slice. It would also
 * compute these quantities at equidistant positions along a depth profile from
 * the surface to the bottom of the model, even if user defined depth bounds were
 * specified in the input file. This was fixed so that adiabatic quantities are
 * computed like all other depth averaged quantities by iterating over all cells
 * in the model domain and averaging the properties spatially. As a consequence
 * some functions in the adiabatic conditions interface that contained the bugs
 * and were only used by the depth average postprocessor have been deprecated.
 * <br>
 * (Rene Gassmoeller, 2023/02/21)
 *
 * <li> New: ASPECT detects when refinement does not actually change the mesh
 * and skips the redundant computation.
 * <br>
 * (Timo Heister, 2023/02/13)
 *
 * <li> Fixed: The Allow fixed temperature on outflow boundaries = false
 * option (and the equivalent option for composition) now also work
 * in time steps where the mesh is refined. In addition, if this
 * parameter is set to false, fixed temperature/composition are now
 * no longer applied to the boundary in cases where the flow normal
 * to the boundary is zero (as is the case in time step 0).
 * <br>
 * (Juliane Dannberg, Rene Gassmoeller, 2023/02/03)
 *
 * <li> New: The PhaseFunction class now includes helper functions
 * n_phases_for_each_composition() and n_phases().
 * <br>
 * (Bob Myhill, 2023/01/22)
 *
 * <li> New: ASPECT has new introspection functions
 * get_indices_for_fields_of_type and get_names_for_fields_of_type.
 * There is also a new utilities function called
 * compute_only_composition_fractions that takes a vector of
 * compositional field values and a list of indices that designates
 * which of those fields correspond to chemical compositions.
 * It then calculates the missing fraction corresponding to the
 * background field, and returns the full vector of fractions.
 * <br>
 * (Bob Myhill, 2023/01/19)
 *
 * <li> New: The "Compositional fields / Types of fields" parameter in ASPECT
 * has a new type called "strain".
 * <br>
 * (Bob Myhill, 2023/01/19)
 *
 * <li> New: There is now a new mesh refinement plugin called
 * "Nonadiabatic temperature threshold" that uses the temperature
 * anomaly (compared to the adiabat) to mark additional cells for
 * refinement.
 * <br>
 * (Juliane Dannberg, 2023/01/17)
 *
 * <li> New: There is now a new additional material model output called
 * ShearHeatingOutputs that contains the fraction of the deformation
 * work that is released as shear heating rather than being
 * converted into other forms of energy, to be filled by the
 * material model. The boundary_area_change_work_fractions variable
 * (which provides the same information) has been removed from the
 * DislocationViscosityOutputs in the grain size material model.
 * <br>
 * (Juliane Dannberg, 2023/01/12)
 *
 * <li> Added: The 'adiabatic' initial temperature plugin can now use a spatially
 * variable top boundary layer thickness read from a data file or specified as a
 * function in the input file. Additionally, the boundary layer temperature can
 * now also be computed following the plate cooling model instead of the
 * half-space cooling model.
 * <br>
 * (Daniel Douglas, John Naliboff, Juliane Dannberg, Rene Gassmoeller, 2022/12/06)
 *
 * <li> New: Added a particle property 'grain size' that tracks grain size
 * evolution on particles using the 'grain size' material model.
 * <br>
 * (Juliane Dannberg, Rene Gassmoeller, 2022/12/06)
 *
 * <li> Added: New transient annulus advection benchmark with analytical solution.
 * <br>
 * (Rene Gassmoeller, 2022/12/02)
 *
 * <li> Improved: The interpolation from particles to fields used to be done field-by-field, creating a lot of overhead in the interpolation algorithms. Now all particle properties are interpolated to fields at the same time, significantly reducing the time required for interpolation for models with many particle properties.
 * <br>
 * (Rene Gassmoeller, 2022/10/10)
 *
 * <li> New: We now support tangential velocity boundary conditions with GMG for more geometries,
 * such as 2D and 3D chunks.
 * <br>
 * (Timo Heister, Haoyuan Li, Jiaqi Zhang, 2022/09/13)
 *
 * <li> Changed: Moved composition velocity statistics postprocessor from cookbook to main code.
 * Compatibility of the postprocessor's input parameter "Names of slab compositional fields" is broken by renaming it to the more general "Names of selected compositional fields".
 * <br>
 * (Elodie Kendall and Anne Glerum, 2022/07/29)
 *
 * <li> Changed: Models that use a strain-dependent rheology
 * can now use an interactive Advection scheme, which was
 * confirmed not to cause unrealistic accumulation of
 * strain in compositional fields during nonlinear iterations.
 * <br>
 * (John Naliboff, 2022/07/22)
 *
 * <li> New: The graphical output of a computation with mesh deformation can now be optionally done on the undeformed reference mesh ("set Output undeformed mesh = true") and the displacement can be written as a separate vector field ("set Output mesh displacement = true").
 * <br>
 * (Timo Heister, 2022/07/16)
 *
 * <li> Changed: In our plugin systems, we have many factory functions
 * that create objects and returned raw pointers; in many but not
 * all places, these returned raw pointers were then converted to
 * `std::unique_ptr`s. This has been changed now
 * by (hopefully) consistently letting factory functions
 * return `std::unique_ptr`s right away, and then dealing with them at
 * the receiving end.
 * <br>
 * (Wolfgang Bangerth, 2022/07/07)
 *
 * <li> New: There is now a postprocessor that outputs the
 * total volume of the computational domain.
 * <br>
 * (Anne Glerum, 2022/06/28)
 *
 * <li> New: Phase transitions can now be deactivated outside a given temperature range specified by
 * upper and lower temperature limits for each phase transition. This allows to implement
 * complex phase diagrams with transitions that intersect in pressure-temperature space.
 *
 * <br>
 * (Haoyuan Li, 2022/05/11)
 *
 * </ol>
 */
