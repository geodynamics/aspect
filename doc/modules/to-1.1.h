/**
 * @page changes_between_1.0_and_1.1 Changes between version 1.0 and version 1.1
 *
 * <p> This is the list of changes made after the release of Aspect version
 * 1.0 for version 1.1. All entries are signed with the names of the author.
 * </p>
 *
 * <ol>
 *
 * <li> Changed: the Stokes linear solver tolerance now includes the current
 * pressure instead of being only based on the rhs. This makes the tolerance
 * less problem dependent. Please note that in many cases you can loosen the
 * tolerance to something between 1e-4 to 1e-7.
 * <br>
 * (Timo Heister, 2014/05/30)
 *
 * <li> New: ASPECT now logs all screen output to output/log.txt.
 * <br>
 * (Timo Heister, 2014/05/29)
 *
 * <li> Changed: The adiabatic reference profile is now defined in a plugin
 * architecture as well. Currently, only the old funtionality of having a
 * constant profile, which is computed before the first timestep, is included.
 * However, now there is the possibility for more simple or more realistic
 * user-written plugins (like an analytic function, or a profile that is
 * updated over time).
 * <br>
 * (Rene Gassmoeller, 2014/05/29)
 *
 * <li> Changed: The initial conditions do not longer use pointers to other
 * parts of the program (like the geometry model), instead they now should be
 * derived from the SimulatorAccess class as all other plugins.
 * <br>
 * (Rene Gassmoeller, 2014/05/29)
 *
 * <li> Changed: For some modules, e.g. the selection of which postprocessor
 * to run, "all" was an option. This may have made sense when we had 3 or 4
 * such postprocessors, but this is no longer the case. Furthermore, not all
 * postprocessors can always be run for a given geometry, input parameters,
 * etc. Consequently, the option of specifying "all" has been removed. Where
 * this applied, we have also removed the default value "all".
 *
 * <li> Changed: If the input file is invalid, for example because a
 * parameter's value did not satisfy its constraints or because it tried to
 * define a parameter that did not exist, all processors in a parallel
 * computation produced the same error message -- leading to massive amounts
 * of output that were barely readable. This has been changed now: only
 * processor zero now generates the output, but all processors abort
 * nonetheless.
 * <br>
 * (Wolfgang Bangerth, 2014/05/29)
 *
 * <li> Changed: Previously, we aborted the program if the specified output
 * directory did not exist. This is inconvenient for jobs that may have been
 * sitting in a queue for a long time and for which forgetting to create the
 * output directory may have just been an oversight. We now simply generate
 * this directory if it is missing.
 * <br>
 * (Wolfgang Bangerth, 2014/05/28)
 *
 * <li> New: Given the proliferation of plugins of all kinds, plugin names are
 * now listed in alphabetical order in the documentation.
 * <br>
 * (Wolfgang Bangerth, 2014/05/27)
 *
 * <li> New: PETSc support can now be enabled using 'cmake -D
 * ASPECT_USE_PETSC=ON'.
 * <br>
 * (Timo Heister, 2014/05/27)
 *
 * <li> Fixed: The GPlates plugin now correctly handles meshes created by
 * GPlates 1.4 and later. Previous Aspect versions may only read in files
 * created by GPlates 1.3.
 * <br>
 * (Rene Gassmoeller, 2014/05/23)
 *
 * <li> Add a find_postprocessor() function in Simulator_access class to pass
 * a pointer of a certain postprocessor, if it is not found, return a NULL
 * pointer.
 * <br>
 * (Siqi Zhang, 2014/05/23)
 *
 *
 * <li> Mesh refinement plugins can now interact with initial global mesh
 * refinement and unflag certain cells if desired.
 * <br>
 * (Timo Heister, 2014/05/22)
 *
 * <li> Changed: The functionality that outputs some basic statistics values
 * like the Rayleigh number in case of a simple material model was moved from
 * the velocity statistics postprocessor to a new postprocessor called basic
 * statistics.
 * <br>
 * (Rene Gassmoeller, 2014/05/21)
 *
 * <li> Fixed: We accidentally evaluated the viscosity of the material model
 * from the place where we compute the adiabatic conditions, but this required
 * information that we didn't have. This is also unnecessary since we don't
 * need the viscosity in this place. This is now fixed.
 * <br>
 * (Wolfgang Bangerth, 2014/05/21)
 *
 * <li> Fixed: Temperature and compositional boundary conditions were
 * previously evaluated only once at the beginning of the simulation, but this
 * did not allow for time dependent boundary conditions for these variables.
 * This is now fixed.
 * <br>
 * (Wolfgang Bangerth, 2014/05/21)
 *
 * <li> Fixed: We forgot to set the initial time before we evaluate the
 * temperature boundary conditions, so temperature boundary conditions could
 * not use <code>get_time()</code> -- they just got NaN. This is now fixed.
 * <br>
 * (Wolfgang Bangerth, 2014/05/21)
 *
 * <li> New: There is now functionality for running models with a free surface
 * with an Arbitrary Lagrangian-Eulerian framework. The user specifies which
 * boundary indicators are to be free surface boundaries, as well as a
 * stabilization parameter that keeps the free surface from exhibiting a
 * sloshing instability. There is an associated postprocessor that calculates
 * mininum/maximum topography, as well as a new section and cookbook in the
 * manual.  Finally, there are two tests for the free surface formulation.
 * <br>
 * (Ian Rose, 2014/05/21)
 *
 * <li>New: There is now a new postprocessor "velocity boundary statistics"
 * that computes statistics about the minimal and maximal velocities on each
 * part of the boundary.
 * <br>
 * (Anne Glerum, 2014/05/21)
 *
 * <li>New: Compositional fields have names now. If no names are specified in
 * the input file, the previously used names C_1, C_2, ... are used as a
 * default.
 * <br>
 * (Juliane Dannberg, 2014/05/21)
 *
 * <li> New: There is now a simple compressible material model with constant
 * compressibility resulting in an exponential dependeny of density on
 * pressure and a linear dependence on temperature deviation from the
 * adiabatic profile. All other material properties are constant.
 * <br>
 * (Rene Gassmoeller, 2014/05/20)
 *
 * <li>New: There is now a new initial temperature condition plugin that can
 * read in solidus temperatures from a file and add perturbations to it. There
 * is also a corresponding data file for initial conditions for Mars from
 * Permentier et al. in <code>data/initial-temperature/solidus.Mars</code>.
 * <br>
 * (Siqi Zhang, 2014/05/20)
 *
 * <li>New: There is now a new postprocessor "spherical velocity statistics"
 * that computes statistics about the radial and tangential velocity field.
 * <br>
 * (Anne Glerum, 2014/05/20)
 *
 * <li>New: There is now a mesh refinement criterium that checks the mesh
 * always has a minimum refinement level determined by a function given in the
 * input file.
 * <br>
 * (Juliane Dannberg, 2014/05/20)
 *
 * <li> New: Input files can now contain lines ending in backslashes to
 * concatenate subsequent lines.
 * <br>
 * (Menno Fraters, 2014/05/19)
 *
 * <li> New: There is a new material model called multicomponent for having an
 * arbitrary number of compositional fields with different material
 * properties, where each compositional field represents a rock type. Within
 * each rock type the material properties are assumed to be constant.  When
 * more than one rock type is present, the material model averages their
 * properties with a weighted arithmetic average. An exception is viscosity,
 * where the user may specify a harmonic, geometric, or arithmetic averaging,
 * or selecting the viscosity of the maximum composition.
 * <br>
 * (Ian Rose, 2014/05/19)
 *
 * <li>New: The dynamic topography postprocessors now take into account (i)
 * the dynamic pressure, (ii) effects due to compressibility, (iii) they now
 * have the option to subtract the mean topography to ensure that the
 * topography is, on average, zero.
 * <br>
 * (Jacqueline Austermann, 2014/05/19)
 *
 * <li>Fixed: The GPlates cookbook in Aspect-1.0 contained a bug that
 * prevented it from using the GPlates plugin. This is fixed now.
 * <br>
 * (Rene Gassmoeller, 2014/05/18)
 *
 * <li>New: There is a new plugin architecture for heating models that allows
 * for plugins defining the radiogenic heating rate in dependency of
 * temperature, pressure, composition, position and time. This introduces an
 * incompatibility of input files to older versions, but the old functionality
 * is preserved and extended.
 * <br>
 * (Rene Gassmoeller, 2014/05/16)
 *
 * <li>New: There is now a material model with reactions between compositional
 * fields and a cookbook describing the usage of this material model.
 * <br>
 * (Juliane Dannberg, Rene Gassmoeller, 2014/05/16)
 *
 * <li>Fixed: Times associated with visualization output are now correctly
 * modified when "Use years in output instead of seconds" is true.
 * <br>
 * (Eric Heien, 2014/05/09)
 *
 * </ol>
 *
 *
 */
