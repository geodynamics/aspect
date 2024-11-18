/**
 * @page changes_between_2.5.0_and_3.0.0 Changes between version 2.5.0 and version 3.0.0
 *
 * <p> This is the list of changes made after the release of ASPECT version
 * 2.5.0 for version 3.0.0. All entries are signed with the names of the author.
 * </p>
 *
 * <ol>
 *
 * <li>  Fixed: The function to determine whether a boundary is an inflow
 * or outflow boundary (for the purpose of deciding whether to apply
 * boundary conditions for compositional fields) did not take into
 * account the mesh deformation. This caused incorrectly applied
 * surface boundary conditions in models with free surface and is
 * fixed now.
 * <br>
 * (Anne Glerum, 2024/11/11)
 *
 * <li> Fixed: Several fixes for the Newton solver that would slow down convergence or
 * rarely lead to wrong results. The Newton solver and defect correction Picard
 * solver computed a wrong pressure scaling for the linear solver. This has been
 * fixed, leading to a reduction in linear solver iterations. The solution
 * variable temporarily contained the update to the solution when using the Newton
 * solver or defect correction Picard solver, leading to bugs in models that
 * assumed the vector contains the values. This has been fixed by always storing
 * the values of the solution in the variable and never the update. The
 * Eisenstat-Walker method to compute linear solver tolerances for the Newton
 * solver would compute too coarse linear residuals leading to poor nonlinear
 * convergence. This has been fixed.
 * <br>
 * (Rene Gassmoeller, Menno Fraters, 2024/11/08)
 *
 * <li> Added: When using the Geodynamic World Builder, ASPECT will create a file
 * called 'original.wb' in the output directory which contains the exact
 * worldbuilder file used to run the ASPECT model.
 * <br>
 * (Daniel Douglas, 2024/11/06)
 *
 * <li> New: The new input parameter "Output directory LFS stripe count" allows for configuring
 * the ASPECT output directory for better performance on Lustre file systems.
 * <br>
 * (Wolfgang Bangerth, 2024/10/22)
 *
 * <li> Changed: The time dimension of the surface condition function of the 'compute profile'
 * adiabatic conditions plugin now respects the value of the parameter 'Use years in output instead of seconds'.
 * <br>
 * (Ranpeng Li, 2024/10/21)
 *
 * <li> Changed: The 'random uniform' particle generator no longer uses
 * the input parameters in the subsection 'Probability density function'.
 * Instead it uses its own subsection 'Random uniform'. Additionally,
 * the default parameters of the 'probability density function' generator
 * have been fixed to allow not specifying a probability density function.
 * If no function is specified the 'probability density function' generator
 * now behaves identically to the 'random uniform' generator.
 * <br>
 * (Rene Gassmoeller, Qianyi Lu, 2024/10/21)
 *
 * <li> Removed: A number of deprecated source code functions and input parameters
 * have been removed. In particular the deprecated option to specify
 * individual material properties in the parameter
 * 'Postprocess/Visualization/List of output variables' has been removed.
 * Input files will be automatically fixed and updated
 * by the update scripts in the directory contrib/utilities.
 * <br>
 * (Rene Gassmoeller, 2024/10/16)
 *
 * <li> New: ASPECT now has a Visual Studio Code extension named
 * ASPECT (https://marketplace.visualstudio.com/items?itemName=zhikui.vscode-aspect).
 * This extension provides syntax highlighting and auto-completion for input parameter
 * files, and can display the documentation of input parameters while writing parameter
 * files in VS Code. It is now one of the recommended extensions for VS Code when using
 * ASPECT.
 * <br>
 * (Zhikui Guo, Timo Heister, Rene Gassmoeller, 2024/10/09)
 *
 * <li> Changed: The class Particle::World has been renamed to
 * Particle::Manager. All references to world have been replaced
 * by references to manager in the source code. The existing header file "world.h"
 * is now deprecated and user plugins should instead include
 * "manager.h". In input and output options we now refer to
 * different groups of particles as "particle systems".
 * <br>
 * (Rene Gassmoeller, 2024/10/08)
 *
 * <li> New: The interface for particle property plugins has been updated.
 * The new interface can update particles cell-wise and is more efficient.
 * In addition, the interface is more extensible for future improvements.
 * All particle property plugins in ASPECT have been updated, user plugins
 * will have to be updated accordingly before the support for the now deprecated
 * old interface is dropped.
 * <br>
 * (Rene Gassmoeller, 2024/09/25)
 *
 * <li> Added: A python script which uses Paraview's pvpython to extract velocities
 * from a global model along user-specified boundaries and saves them to ASCII
 * files which can then be applied as velocity boundary conditions in regional
 * chunk models.
 * <br>
 * (Daniel Douglas, 2024/09/05)
 *
 * <li> New: The grain size evolution functionality has been extracted from the grain
 * size material model into the new grain size evolution reaction model. This
 * simplifies using the grain size evolution equations in other material models.
 * <br>
 * (Rene Gassmoeller, 2024/07/26)
 *
 * <li> New: The CompositeViscoPlastic rheology in ASPECT now
 * allows the user to select isostress (Reuss) or
 * isostrain (Voigt) averaging of viscosities for
 * multicomponent materials via the parameter
 * "Viscosity averaging scheme".
 * <br>
 * (Bob Myhill, 2024/07/25)
 *
 * <li> Added: The functions in the CompositeViscoPlastic
 * rheology model now have grain_size as an argument.
 * <br>
 * (Bob Myhill, 2024/07/17)
 *
 * <li> Changed: The CompositeViscoPlastic rheology in ASPECT
 * now uses a log-stress Newton scheme, which requires
 * fewer iterations to converge.
 * <br>
 * (Bob Myhill, 2024/07/17)
 *
 * <li> Changed: The CompositeViscoPlastic rheology in ASPECT is now
 * required to use the approximate Peierls flow rheology if
 * Peierls flow is switched on.
 * <br>
 * (Bob Myhill, 2024/07/17)
 *
 * <li> Changed: The CompositeViscoPlastic rheology in ASPECT now uses
 * the DruckerPragerPower rheology.
 * <br>
 * (Bob Myhill, 2024/07/17)
 *
 * <li> New: ASPECT now has a DruckerPragerPower rheology. This rheology
 * is based on the DruckerPrager rheology, but the strain-rate-independent
 * yield stress is replaced with a power-law rheology. The parameters
 * "Reference plastic strain rate" and "Plastic stress exponent" control
 * how the yield stress varies with strain rate.
 * <br>
 * (Bob Myhill, 2024/07/17)
 *
 * <li> Added: The DiffusionCreep rheology module can now
 * calculate viscosities and strain rates with grain size
 * as an variable argument.
 * <br>
 * (Bob Myhill, 2024/07/04)
 *
 * <li> Changed: The particle interpolator plugin 'bilinear least squares'
 * now activates a limiter by default, which prevents unwanted over-
 * and undershoots by limiting the returned properties to the bounds
 * of the particle properties in the cell. The limiter can be
 * deactivated by the input parameter 'Use linear least squares limiter'.
 * <br>
 * (Rene Gassmoeller, 2024/06/25)
 *
 * <li> Removed: The input parameter 'Update ghost particles' has been deprecated.
 * Ghost particles are now always updated. Existing parameter files can
 * be updated by the script in 'contrib/utilities/update_prm_files.sh'.
 * <br>
 * (Rene Gassmoeller, 2024/06/25)
 *
 * <li> Added: The Frank Kamenetskii module now has user inputs for reference
 * temperature and reference pressure. Adiabatic surface temperature and pressure
 * are used as the defaults for those values.
 * <br>
 * (Grant Block, 2024/06/24)
 *
 * <li> New: ASPECT now supports compositional fields with
 * different discretizations (continuous or discontinuous)
 * and different polynomial degrees at the same time.
 * <br>
 * (Timo Heister, 2024/06/18)
 *
 * <li> Changed: The particle plugin functions declare_parameters and parse_parameters
 * now assume that they are already in the particles subsection. User plugins will
 * have to change accordingly.
 * <br>
 * (Menno Fraters, 2024/06/15)
 *
 * <li> Changed: The operator splitting scheme can now use SUNDIALs
 * ARKode to compute the reactions (rather than a hand-written
 * forward Euler scheme). ARKode is the new default, but the
 * original scheme can still be used by setting the input
 * parameter Reaction solver type to fixed step.
 * <br>
 * (Juliane Dannberg, 2024/06/14)
 *
 * <li> Added: A postprocessor 'darcy velocity' that visualizes the fluid
 * velocity in models with a compositional field named 'porosity'.
 * This postprocessor is independent of the method the porosity is actually
 * advected with, i.e. it could be static, advected with the Darcy velocity,
 * or advected according to the two-phase flow equations.
 * <br>
 * (Daniel Douglas, 2024/06/13)
 *
 * <li> New: There is now a postprocessor that computes statistics
 * about how many iterations ARKode needs to solve ODEs.
 * <br>
 * (Juliane Dannberg, 2024/06/12)
 *
 * <li> Fixed: The grain size material model now computes the
 * stress that is used to calculate the grain size reduction
 * rate with the viscosity that has been limited between its
 * minimum and maximum given in the input file, rather than
 * with the unlimited viscosity. This is more consistent,
 * because the limited viscosity is the one we use to compute
 * the strain rate (and stress is 2*viscosity_strain_rate).
 * This updated implementation prevents artificial grain size
 * reduction in regions with large viscosity.
 * <br>
 * (Juliane Dannberg, 2024/06/10)
 *
 * <li> New: The CPO particle property plugin can now use the Geodynamic
 * World Builder to set initial grain-sizes and grain orientations
 * <br>
 * (Menno Fraters, 2024/06/09)
 *
 * <li> Added: A new particle interpolator plugin 'distance weighted average' that
 * computes the interpolated properties as an average of particle properties
 * weighted by a distance function. This plugin is similarly accurate
 * as the 'cell average' interpolator, but provides smoother results and is
 * less susceptible to jumps in interpolated properties if particles move
 * small distances.
 * <br>
 * (Rene Gassmoeller, 2024/06/09)
 *
 * <li> Added: Functionality to limit the timestep based on the
 * Darcy velocity.
 * <br>
 * (Daniel Douglas, 2024/06/09)
 *
 * <li> New: There is now a new cookbook which models compressible
 * mantle convection in an annulus.
 * <br>
 * (Cedric Thieulot, 2024/06/06)
 *
 * <li> Changed: The grain size material model no longer has the
 * option to advect the logarithm of the grain size. This
 * option used to be helpful for models with large gradients
 * in grain size, but now that ASPECT has particles it is
 * more accurate to advect the grain size on particles
 * instead.
 * <br>
 * (Juliane Dannberg, 2024/06/07)
 *
 * <li> Changed: The grain size material model no longer has the
 * option to scale the grain size in the lower mantle. This
 * option used to be helpful in models where the grain size
 * is very different between upper and lower mantle, but
 * now that ASPECT has particles, it is more accurate to
 * advect the grain size on particles instead.
 * <br>
 * (Juliane Dannberg, 2024/06/07)
 *
 * <li> New: The CPO particle property plugin set now contains a plugin
 * which can compute the symmetry decompositions of the elastic
 * tensor. There are now also tests which test the whole CPO
 * particle property plugin set.
 * <br>
 * (Menno Fraters, 2024/06/06)
 *
 * <li> New: Added a Rayleigh Taylor instability benchmark with a free surface.
 * The benchmark was described by Kaus et al. (2010) and reproduced with
 * ASPECT with different stabilization schemes for the free surface by
 * Rose et al. (2017). It can be used to show that the stabilization
 * scheme implemented by I. Rose allows for larger timesteps, and that
 * the projection of velocities onto the surface leads to significantly
 * different topography.
 * <br>
 * (Anne Glerum, 2024/06/06)
 *
 * <li> Changed: The grain size material model now uses SUNDIALs
 * to compute the grain size change (rather than a hand-
 * written forward Euler scheme).
 * <br>
 * (Juliane Dannberg, 2024/06/06)
 *
 * <li> New: Added fluid reaction and properties calculations to the
 * Katz2003MantleMelting model in the ReactionModels module. These calculations
 * were taken from MeltSimple and the model now calls Katz2003MantleMelting. Additionally,
 * the katz 2003 option in ReactiveFluidTransport calls these calculations, allowing for
 * melting and melt transport to be composited with any material model via ReactiveFluidTransport.
 * <br>
 * (Grant Block, 2024/06/06)
 *
 * <li> New: There is now a python script in the contrib folder that
 * creates netcdf files from an input ascii file
 * formatted following structured ascii data file convention.
 * <br>
 * (Arushi Saxena, 2024/06/06)
 *
 * <li> Changed: We no longer support "make test" to run the
 * testsuite. Please use "ctest" directly.
 * <br>
 * (Timo Heister, 2024/06/05)
 *
 * <li> New: There is now a new cookbook based on the
 * setup of Allken et al, G3, 2012 which models
 * rift interaction in 3D brittle-ductile coupled systems.
 * <br>
 * (Cedric Thieulot, 2024/06/05)
 *
 * <li> Added: A cookbook using a model with simple 2D annulus convection
 * to showcase the capabilities of the pyvista Python package in plotting
 * and conducting simple mesh operations using numpy arrays. The cookbook
 * contains a parameter file, documentation, and figures produces by
 * manipulating data on pyvista and plotting new fields. The relevant script
 * is located in [cookbooks/twoD_annulus_visualization/2D_annulus_example.prm].
 * <br>
 * (Madeleine Kerr, 2024/06/05)
 *
 * <li> Added: Cookbook which showcases the use of the tian approximation
 * in the Reactive Fluid Transport material model. The cookbook features
 * a kinematically subducting slab with an initial hydration state
 * dehydrating as it advects through a hot mantle wedge.
 * <br>
 * (Daniel Douglas, 2024/06/05)
 *
 * <li> Fixed: There was a bug in spherical harmonic data reading of S40RTS
 * and SAVANI initial temperature models when a maximum degree lower
 * than the maximum degree of the spherical harmonic coefficients of
 * the data file is specified by users to set up temperature field
 * from spherical harmonic summation. This bug is now fixed.
 * <br>
 * (Shangxin Liu, 2024/06/04)
 *
 * <li> Changed: In S40RTS and SAVANI initial temperature models, the input parameter,
 * 'Specify a lower maximum order' is changed to 'Specify a lower maximum degree',
 * and the input parameter, 'Maximum order' is changed to 'Maximum degree'.
 * This avoids confusion because the function of these two parameters is to
 * allow users to specify a maximum degree lower than the maximum degree
 * of the spherical harmonic data file to calculate the temperature field.
 * <br>
 * (Shangxin Liu, 2024/06/04)
 *
 * <li> Changed: The QT based Parameter GUI has been removed from the ASPECT repository
 * and the documentation as it was no longer maintained and no longer working
 * correctly. Instead, use the parameter list linked from the website or edit
 * directly in Visual Studio Code with the ASPECT plugin.
 * <br>
 * (Timo Heister, 2024/06/02)
 *
 * <li> Added: The Function expression plug-in for the Gravity model can now evaluate
 * functions in multiple coordinate systems. The plug-in was amended to read in entries
 * to the "Coordinate system" parameter in the Function subsection and Gravity Model
 * section of the parameter file. There are three added tests for the Gravity model
 * function plug-in to ensure the cartesian, spherical, and depth coordinate systems
 * evaluate functions as expected in a 2D cylindrical annulus.
 * <br>
 * (Madeleine Kerr, 2024/06/02)
 *
 * <li> Added: The heating function plugin now supports input in different coordinate
 * systems (cartesian, spherical, depth).
 * <br>
 * (Rene Gassmoeller, 2024/06/02)
 *
 * <li> New: Benchmark from Gassmoller et al. 2018 for particle integration schemes.
 * <br>
 * (Gabriel Johnston, 2024/06/01)
 *
 * <li> Improved: The limiter for the DG solution can now be selected
 * for each compositional field individually, instead of using
 * a single setting that applies to all compositional fields.
 * If only a single setting is found it will still apply to all fields.
 * <br>
 * (Rene Gassmoeller, 2024/06/01)
 *
 * <li> Added: Functionality for multiplying the dislocation and/or diffusion creep
 * viscosities by prefactors which depend on the model composition. Created a
 * specific use-case for multiplying the viscosity by a factor that accounts for
 * bound H2O from Hirth & Kohlstaedt 2004 10.1029/138GM06
 * <br>
 * (Daniel Douglas, 2024/06/01)
 *
 * <li> New: The grain size material model now has the option
 * to apply plastic yielding (Drucker-Prager yield
 * criterion) to limit the viscosity. This can be switched
 * on through a new input parameter. At the moment, grain
 * size evolution is not affected by plastic yielding.
 * <br>
 * (Juliane Dannberg, 2024/06/01)
 *
 * <li> Changed: Melt fraction calculated by the Katz2003 model
 * in the ReactionModels module was added to the reactive fluid transport model.
 * A test was added compositing the reactive fluid transport and visco plastic models,
 * benchmarking the melt fraction results with Katz et. al. 2003.
 * <br>
 * (Grant Block, 2024/06/01)
 *
 * <li> New: ASPECT now outputs the physical units of quantities into .pvtu files.
 * <br>
 * (Wolfgang Bangerth, 2024/06/01)
 *
 * <li> Changed: ASPECT's default Stokes preconditioner has been changed from
 * a block AMG preconditioner to the geometric multigrid (GMG) preconditioner
 * described in Clevenger and Heister, 2021 (https://doi.org/10.1002/nla.2375).
 * The GMG preconditioner generally performs better, but requires that
 * the viscosity is averaged using one of the material model averaging functions.
 * If models use features that are not supported by the GMG preconditioner
 * and the GMG preconditioner has not been explicitly set ASPECT will
 * fall back to the AMG preconditioner.
 * <br>
 * (Rene Gassmoeller, 2024/06/09)
 *
 * <li> Changed: The grain size material model now uses the
 * phase function class to compute where phase transitions
 * occur. This means that phase transition temperature,
 * depth and Clapeyron slope now have to be specified in
 * a different format in the input file (rather than a
 * comma-separated list, they can now be specified as
 * key:value1|value2|...). This is in line with other
 * material models. Note that the width of phase
 * transitions is still zero in the grain size model.
 * <br>
 * (Juliane Dannberg, 2024/05/31)
 *
 * <li> Changed: The reset of the grain size at phase transitions
 * is now determined only by the velocity of the material
 * moving through a transition rather than by a somewhat
 * arbitrary 'transition width' parameter.
 * <br>
 * (Juliane Dannberg, 2024/05/30)
 *
 * <li> New: ASPECT now has a DiffusionDislocation rheology. This rheology
 * is identical to the one that was originally in the Material Model
 * of the same name.
 * <br>
 * (Bob Myhill, 2024/05/29)
 *
 * <li> Changed: A pressure term was added to the Frank Kamenetskii viscous flow law
 * in the Visco Plastic material model. The new FK function is given by
 * viscosity = A * exp(E * 0.5 * (1.0-(T/ref_T)) + F * (P-ref_P)/(ref_rho*g*h)),
 * where F * (P-ref_P)/(ref_rho*g*h) are the new terms added here.
 * <br>
 * (Grant Block, 2024/05/29)
 *
 * <li> Changed: The A block solver in the expensive Stokes iterations
 * has been changed for free surface models from a CG solver to a BiCGStab
 * solver. This should improve the stability of the inner (top left)
 * preconditioner in free surface models. The new solver can also be forced
 * by an input parameter, which can be a useful option for models that
 * crash when solving the top left (A) block of the Stokes system.
 * <br>
 * (Rene Gassmoeller, 2024/03/31)
 *
 * <li> Changed: The solver used during the expensive iterations of
 * the AMG preconditioner has been changed from the trilinos implementation
 * to the one provided by deal.II. The new solver should be faster and more
 * stable, but can produce marginally different results.
 * <br>
 * (Rene Gassmoeller, 2024/03/26)
 *
 * <li> New: It is now possible to use the "initial topography" plugin system
 * for the spherical shell geometry, including using the ascii data
 * system that can then be used to import a text description of the
 * initial topography. As an example, we also include the file
 * `data/geometry-model/initial-topography-model/ascii-data/global_1deg.txt.gz`
 * that describes the topography of Earth on a (relatively coarse) 1
 * degree mesh.
 * <br>
 * (Wolfgang Bangerth, 2024/03/21)
 *
 * <li> Changed: The postprocessor 'heat flux map' now outputs the heat flux values in
 * the same format as the 'dynamic topography' postprocessor does for dynamic
 * topography. This in particular means instead of writing all values into a
 * single file 'heat flux map' now writes one file per boundary. Both
 * postprocessors now take care to not include spaces in the data file
 * headers so that separating the header columns by spaces works correctly.
 * <br>
 * (Rene Gassmoeller, 2024/03/15)
 *
 * <li> Removed: The outdated input option for the geoid postprocessor named
 * 'Include the contributon from dynamic topography' has been removed.
 * Use 'Include surface topography contribution' or 'Include CMB topography
 * contribution' instead.
 * <br>
 * (Rene Gassmoeller, 2024/03/15)
 *
 * <li> Implemented a new fluid-rock interaction scheme for the reactive fluid
 * transport material model based on a parametrization by Tian et al., 2019 G3,
 * https://doi.org/10.1029/2019GC008488.
 * <br>
 *  (Daniel Douglas, 2024/02/24)
 *
 * <li> Incompatible: When using the "AsciiDataBoundary" class to read initial
 * topographies for the chunk geometry, the input file needed to provide
 * topographies spaced on a grid in latitude and longitude where the
 * longitude was to be provided from -pi (180 degrees west) to +pi (180
 * degrees east). This is contrary to the rest of ASPECT, which generally
 * uses 0 to 2*pi (0 to 360 degrees east) as a convention. This has now
 * been changed: The initial topography plugin will be asked about
 * topographies in lat-long with longitudes between 0 and 2*pi. This also
 * means that input files have to be changed correspondingly.
 * <br>
 * (Wolfgang Bangerth, 2024/02/16)
 *
 * <li> New: There is now a new cookbook that demonstrates how
 * to use the Geodynamic World Builder to set up initial
 * conditions for a mid-ocean ridge model with a transform
 * fault following the setup of Behn et al., 2007.
 * <br>
 * (Juliane Dannberg, 2024/02/13)
 *
 * <li> New: There is now a new visualization postprocessor called "surface
 * elevation" that can be used to plot the initial topography or the
 * effects of a free surface.
 * <br>
 * (Wolfgang Bangerth, 2024/02/04)
 *
 * <li>  New: Solution variables can now also be outputted
 *  on the surface mesh of the model domain.
 *  <br>
 *  (Anne Glerum, 2024/05/30)
 *
 * <li>  Improved: The geometry function point_is_in_domain
 *  now also works for meshes that include initial topography
 *  and/or mesh deformation for box and chunk geometries.
 *  <br>
 *  (Anne Glerum, 2023/12/14)
 *
 * <li> Removed: The GeometryModel::Chunk and GeometryModel::TwoMergedChunks
 * classes had member functions `depth_wrt_topo()`, but these were not
 * used anywhere and did not work with the mesh deformation
 * framework. They have consequently been removed.
 * <br>
 * (Wolfgang Bangerth, 2023/12/07)
 *
 * <li> Changed: The interface
 * Postprocess::VisualizationPostprocessors::CellDataVectorCreator::execute()
 * used to return an object of type
 * `std::pair<std::string, Vector<float>*>` where the second part is just a
 * raw pointer. This has been changed (incompatibly) to
 * `std::pair<std::string, std::unique_ptr<Vector<float>>>` to avoid the use
 * of raw pointers and ensure that memory de-allocation happens automatically.
 * <br>
 * (Wolfgang Bangerth, 2023/11/28)
 *
 * <li> Fixed: The stress and shear stress visualization postprocessors
 * now output the correct stresses when elasticity is enabled.
 * <br>
 * (Bob Myhill, Rebecca Fildes, 2023/11/10)
 *
 * <li> Fixed: The 'principal stress' postprocessor used the wrong
 * sign for stresses in models with elastic deformation, leading
 * to principal stress directions that were rotated by 90 degrees.
 * This is fixed now.
 * <br>
 * (Rene Gassmoeller, Rebecca Fildes, 2023/11/09)
 *
 * <li> Changed: ASPECT now considers boundaries with no normal flow as
 * boundaries with inflow for the purposes of the parameter
 * 'Allow fixed temperature on outflow boundaries' and the
 * corresponding parameter for composition. This was the default
 * behavior up to ASPECT 2.4.0. This behavior was changed in ASPECT 2.5.0,
 * in which boundaries with no normal flow are treated like outflow
 * boundaries. The new behavior caused unintended side effects, therefore
 * it is reverted back to the original behavior.
 * The reason for the initial change was a bugfix for boundary conditions
 * in the first timestep that is now implemented differently.
 * <br>
 * (Rene Gassmoeller, 2023/11/06)
 *
 * <li> Fixed: The RK4 particle interpolation scheme computed wrong
 * particle locations for particles that crossed periodic boundaries
 * in a box geometry. This is fixed now.
 * <br>
 * (Rene Gassmoeller, 2023/11/03)
 *
 * <li> New: There is now a new cookbook that demonstrates how
 * to use the grain size material model and how to choose
 * the particle parameters in a model with a complex
 * nonlinear rheology.
 * <br>
 * (Juliane Dannberg, 2023/10/07)
 *
 * <li> New: The amount of shear heating in a computation can
 * now be limited by setting a maximum stress to be used,
 * based on a Drucker-Prager yield stress with a user-
 * specified cohesion and friction angle. This is useful
 * in models with unrealistically high stresses.
 * <br>
 * (Juliane Dannberg, 2023/10/01)
 *
 * <li> Added: There is now a prescribed viscosity material model plugin, which can
 * overwrite the viscosity of a specified location with a viscosity prescribed
 * by a function.
 * <br>
 * (Menno Fraters, 2023/09/27)
 *
 * <li> Changed: ASPECT has been renamed from "Advanced Solver for Problems in
 * Earth's ConvecTion" to "Advanced Solver for Planetary Evolution,
 * Convection, and Tectonics" to reflect that the scope of ASPECT has
 * grown beyond mantle convection.
 * <br>
 * (Timo Heister, 2023/09/17)
 *
 * <li>  New: There is now an entropy table lookup model
 * in the initial composition model.
 * This model takes an initial temperature field and
 * converts it to an initial entropy field using
 * a lookup table.
 * <br>
 * (Haoyuan Li, 2023/08/11)
 *
 * <li> Changed: The input subsection 'Particles' has been moved
 * from the 'Postprocess' subsection into its own subsection
 * at the top level. This is to make it more clear that
 * particles are not simply a postprocessing feature.
 * All included input files have been adjusted and the reformatting
 * script has been adjusted. Users will have to apply the reformatting
 * script to their input files manually.
 * <br>
 * (Rene Gassmoeller, 2023/09/10)
 *
 * <li> Changed: ASPECT now requires CMake 3.13.4 (just like deal.II v9.5.0)
 * and no longer supports GeodynamicWorldBuilder older than 0.5.0.
 * <br>
 * (Rene Gassmoeller, 2023/09/09)
 *
 * <li> Changed: The benchmarks 'blankenbach' and 'king2dcompressible' have been
 * reworked and extended. Changes include: (1) higher resolution, (2) more
 * reproducible convergence results, (3) additional figures, (4) a comparison
 * between gradient based and CBF heat flux calculation, (5) reproducible
 * Richardson extrapolation of ASPECT reference results.
 * <br>
 * (Rene Gassmoeller, 2023/08/29)
 *
 * <li> New: There is now a new viscosity profile data file that
 * is more consistent with how the profile is computed in
 * the original Steinberger and Calderwood (2006) paper.
 * <br>
 * (Juliane Dannberg, Bernhard Steinberger, 2023/08/22)
 *
 * <li> New: The steinberger material model now allows the
 * viscosities of different materials to differ from
 * the reference viscosity by a constant factor
 * via the parameter Composition viscosity prefactors.
 * <br>
 * (Poulami Roy and Bob Myhill, 2023/08/14)
 *
 * <li>  Improved: The 'spherical constant' boundary composition and boundary
 * temperature plugins have been updated. In particular the boundary
 * composition plugin can now prescribe different values for different
 * compositional fields.
 * <br>
 * (Rene Gassmoeller, 2023/08/02)
 *
 * <li> Improved: The iteration scheme in the diffusion dislocation material
 * model now uses the derivative of the logarithm of the strain rate to the
 * logarithm of the stress. This improves efficiency and stability of the
 * iterative scheme. The parameter "Strain rate residual tolerance" now
 * corresponds to the residual for the logarithm of the strain rate (which
 * close to the solution corresponds to the relative residual).
 * <br>
 * (Haoyuan Li and Bob Myhill, 2023/07/18)
 *
 * <li>  New: There is now a strain rate tensor postprocessor
 * that is only outputted on the surface of the model
 * domain.
 * <br>
 * (Anne Glerum, 2023/07/17)
 *
 * <li>  Added: A sea level postprocessor that computes the
 * sea level for glacial isostatic adjustment modeling.
 * It computes the sea level based on the free
 * surface topography, ocean basin, ice melt from
 * boundary traction, and perturbed gravitational
 * potential of the Earth model from the geoid
 * postprocessor.
 * <br>
 * (Maaike Weerdesteijn, 2023/07/14)
 *
 * <li> New: There is now a new material model ('reactive fluid
 * transport') that is designed to advect fluids and
 * compute fluid release and absorption based on
 * different models for fluid-rock interaction. The
 * properties of the solid are taken from a
 * base material model.
 * <br>
 * (John Naliboff, 2023/07/14)
 *
 * <li> Fixed: The strain rheology now uses the correct 'old strain'
 * when computing the reaction term updates.
 * <br>
 * (John Naliboff 2023/07/14)
 *
 * <li> New: There is now a new cookbook which models
 * the (Poiseuille) flow of the lower crust around a rigid obstacle.
 * <br>
 * (Cedric Thieulot, 2023/07/13)
 *
 * <li>  Improved: The initial lithostatic pressure plugin for
 *  boundary tractions now includes the initial topography
 *  of the reference point into its pressure profile
 *  instead of the maximum topography within the domain.
 *  <br>
 *  (Anne Glerum, 2023/07/13)
 *
 * <li> Changed: Boundary traction models are now organized in
 * a manager class that allows to assign multiple
 * boundary traction plugins to each boundary. Existing
 * user plugins that use the function
 * this->get_boundary_traction() will need to be modified
 * to use this->get_boundary_traction_manager().
 * <br>
 * (Chameera Silva, Rene Gassmoeller, 2023/07/13)
 *
 * <li> New: Added a cookbook on how to set up global instantaneous models
 * based on recent geophysical constraints with a heterogeneous density and
 * viscosity distribution and weak plate boundaries prescribed using
 * different plate boundary configurations.
 * <br>
 * (Arushi Saxena, Juliane Dannberg, and Rene Gassmoeller 2023/07/13)
 *
 * <li> New: There is now a new cookbook which is a very
 * simple subduction initiation model as published
 * by Matsumoto and Tomoda (1983).
 * <br>
 * (Cedric Thieulot, 2023/07/12)
 *
 * <li> Changed: The entropy plugin now includes a plasticity rheology, which
 * can be set by cohesion and frictional angle. It can also use a viscosity
 *  prefactor profile, set the lateral viscosity variation limit, and choose
 * p-T dependent conductivity.
 * <br>
 * (Ranpeng Li & Juliane Dannberg, 2023/07/13)
 *
 * <li> Changed: Material models are now always provided with the current
 * strain rate. Checks if the strain rate is provided are not longer necessary.
 * MaterialModelInputs can no longer be initialized without a strain rate.
 * <br>
 * (Rene Gassmoeller, 2023/07/12)
 *
 * <li> New: A cookbook which replicates the model setup of the
 * van Keken et al., 2008 2D corner flow subduction benchmark.
 * <br>
 * (Daniel Douglas, Cedric Thieulot, Wolfgang Bangerth, Max Rudolph, 2023/07/12)
 *
 * <li> New: There is now a new cookbook which models
 * the deformation of elliptical and rectangular inclusions
 * in simple shear and pure shear.
 * <br>
 * (Cedric Thieulot, 2023/07/11)
 *
 * <li> Fixed:
 * When using the 'artificial viscosity' visualization postprocessor in parallel,
 * we used to put NaNs into the output vector (for ghost and artificial
 * cells) that then triggered floating point exceptions. This is now fixed.
 * <br>
 * (Wolfgang Bangerth, 2023/07/10)
 *
 * <li> New: The base code for the computation and output of
 * Crystal Preferred Orientation (CPO) has been added to
 * ASPECT.
 * <br>
 * (Menno Fraters, 2023/07/09)
 *
 * <li>  Changed: The input parameter 'Include viscoelasticity' has been
 * removed from the visco_plastic material model. When
 * 'Enable elasticity' in the Formulation subsection is
 * switched on, elasticity will automatically be included
 * in the visco_plastic material model.
 * <br>
 * (Anne Glerum, 2023/07/09)
 *
 * <li> New: ASPECT now allows it to select a different list of
 * assemblers for each advection field.
 * <br>
 * (Juliane Dannberg, 2023/07/09)
 *
 * <li> Improved: The iteration scheme in the Peierls creep rheology
 * now uses the derivative of the logarithm of the strain rate to
 * the logarithm of the stress.
 * <br>
 * (Haoyuan Li and Bob Myhill, 2023/07/08)
 *
 * <li>  New: The cutoff stress for the Peierls stress
 * could now be applied strictly as the lower
 * limit of the stress in the Peierls rheology
 * <br>
 * (Haoyuan Li, Bob Myhill, 2023/07/08)
 *
 * <li> Fixed: The computation of the upper box origin in the 'box
 * with lithosphere boundary indicators' geometry model did not
 * take into account a nonzero origin. This is fixed now.
 * <br>
 * (Juliane Dannberg, 2023/07/08)
 *
 * <li> Changed: Default field types are now defined by ASPECT
 * based on Field Names if the types are not defined by the user.
 * The logic is as follows: if the name contains the substring stress,
 * it has type stress. If the name is equal to
 * s11, s12, s22, s13, s23, s33, it has type strain.
 * If the name is equal to grain_size, porosity, density_field,
 * or entropy it has the types grain size, porosity, density or entropy
 * respectively. Otherwise it has the type chemical composition.
 * <br>
 * (Bob Myhill, 2023/07/08)
 *
 * <li> New: There is now a new additional benchmark called solubility
 * that demonstrates the mass of water is conserved as water is
 * released, migrates and is being reabsorbed.
 * <br>
 * (Juliane Dannberg, 2023/07/08)
 *
 * <li>  Improved: ASPECT is now by default compiled in DebugRelease
 * mode which compiles both a debug and a release (optimized)
 * version of the executable. Individual build types can still
 * be selected using 'make debug' or 'make release' or by
 * setting the cmake variable CMAKE_BUILD_TYPE to Debug or
 * Release.
 * <br>
 * (Rene Gassmoeller, Timo Heister, 2023/07/07)
 *
 * <li> New: The maximum relative increase in the time step length is now
 * bounded by a factor of 1.91, where before it was unbounded. This value
 * is obtained from theoretical considerations about parabolic problems.
 * <br>
 * (Laila Busaleh, 2022/09/15)
 *
 * <li> New: Mesh deformation plugin that uses the landscape
 * evolution code FastScape to deform the surface through
 * erosion and sediment deposition.
 *
 * Citations:
 * Neuharth, D., Brune, S., Wrona, T., Glerum, A., Braun, J., & Yuan, X. (2022). Evolution of rift systems and their fault networks in response to surface processes. Tectonics, 41(3), e2021TC007166.
 *
 * Neuharth, D., Brune, S., Glerum, A., Morley, C. K., Yuan, X., & Braun, J. (2022). Flexural strike-slip basins. Geology, 50(3), 361-365.
 *
 * <br>
 * (Derek Neuharth, Anne Glerum, Sascha Brune, Esther Heckenbach 2023/03/01)
 *
 * <li> New: There is now a new initial temperature plugin that adds
 * a fixed number of Gaussian perturbations placed at random
 * locations to the initial temperature.
 * <br>
 * (Juliane Dannberg, 2023/02/10)
 *
 * <li> New: The 'grain size' material model now supports the grain size evolution
 * equations of Mulyukova & Bercovici 2018 in addition to the existing
 * formulations. For this purpose the old input parameter 'Use paleowattmeter'
 * has been deprecated and replaced with 'Grain size evolution formulation'. The
 * default behavior of the material model has not changed. Benchmarks and tests
 * for the new feature have been added.
 * <br>
 * (Arushi Saxena, Ranpeng Li, Juliane Dannberg, Rene Gassmoeller, 2023/08/23)
 *
 * </ol>
 */
