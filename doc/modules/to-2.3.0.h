/**
 * @page changes_between_2.2.0_and_2.3.0 Changes between version 2.2.0 and version 2.3.0
 *
 * <p> This is the list of changes made after the release of ASPECT version
 * 2.2.0 for version 2.3.0. All entries are signed with the names of the author.
 * </p>
 *
 * <ol>
 *
 * <li> <li> New: A mesh refinement plugin that allows to set regions of minimum and
 * maximum refinement based on specific values of variables (of e.g. temperature) between two
 * isosurfaces of that value. For example the minimum refinement is 2 and maximum refinement
 * is 4 if the temperature is between 273 K and 1600 K.
 * <br>
 * (Menno Fraters and Haoyuan Li, 2021/07/09)
 *
 * <li> Fixed: The Steinberger material model now uses the correct profile for
 * lateral viscosity variations (in dependence of temperature). Before
 * this change, it had artificial jumps at the transition from upper to
 * lower mantle and used n=3 instead of n=3.5.
 * <br>
 * (Juliane Dannberg, 2021/05/20)
 *
 * <li> Fixed: The 'random uniform' and 'probability density function' particle
 * generator plugins crashed when one MPI rank had no cells at all (e.g. when
 * running a very coarse model on a large number of processors). This has been
 * fixed.
 * <br>
 * (Rene Gassmoeller, 2021/02/01)
 *
 * <li> New: Mesh deformation now also works in combination with particles.
 * Instead of at the end of the timestep, particles are now advected
 * before the solving of the compositional field advection equations.
 * The values of fields that are advected through the particle method
 * are now in sync with the values carried by the particles.
 * In iterative advection schemes, the particle location is restored
 * before each iteration.
 * A new postprocessor outputs the reference location of the particles
 * in their respective cells.
 * <br>
 * (Anne Glerum, Rene Gassmoeller, Robert Citron, 2021/02/04)
 *
 * <li> Changed: ASPECT now uses a mapping caching strategy provided by deal.II for
 * models with curved cells (spherical shells, chunk, sphere) as long as there is
 * no mesh deformation. The new mapping drastically decreases the time for sorting
 * particles, but has little influence (as far as we know) on other timings, and
 * has no influence on the accuracy.
 * <br>
 * (Rene Gassmoeller, 2021/02/01)
 *
 * <li> Changed: The multicomponent_compressible material model
 * now treats the compositional field values as mass fractions, and use
 * these mass fractions to calculate material properties. The specific heat
 * is arithmetically averaged according to mass fraction, while the density,
 * compressibility, thermal expansivity, and thermal conductivity are
 * all arithmetically averaged according to volume. With the exception of
 * thermal conductivity, these choices are dictated by mass conservation
 * and thermodynamics. Viscosity is averaged according to volume fractions,
 * with the scheme chosen by the user (this functionality has not changed).
 * These changes will probably not be noticeable for model runs using
 * compositional fields with values of either 1 or 0 or materials with
 * the same or very similar densities, but significant differences may be
 * seen for models using distinct material properties. If users have used the
 * compute_volume_fractions function in utilities, they should be aware that
 * the function is now called compute_composition_fractions.
 * <br>
 * (Bob Myhill, 2021/01/27)
 *
 * <li> Changed: The internal pressure scaling of ASPECT no longer uses the 'reference
 * viscosity' provided by the material model as a representative viscosity.
 * Instead it computes the representative viscosity as an average of the logarithm
 * of the actual viscosities present in the model.  This prevents accidental
 * misconfigurations when the reference viscosity is vastly different from the
 * actual viscosities or when viscosities change drastically during a model run.
 * The 'reference_viscosity' function in the material models is now largely
 * ignored and will be removed in a future release.
 * <br>
 * (Rene Gassmoeller, 2021/01/18)
 *
 * <li> New: There is a new nullspace removal option called 'net surface rotation',
 * which removes the net rotation of the surface. This can be useful to compare
 * plate velocities on a sphere with plate reconstructions in a no net rotation
 * reference frame. Additionally, the 'rotation statistics' postprocessor was
 * extended to also compute the net surface rotation.
 * <br>
 * (Rene Gassmoeller, 2021/06/20)
 *
 * <li> Changed: The AsciiDataLookup class has been renamed to StructuredDataLookup to
 * reflect that other file formats will be supported in the future.
 * <br>
 * (Timo Heister, 2021/01/05)
 *
 * <li> New: Particle advection can now be used in combination with
 * the repetition of timesteps. Before each repetition the particles
 * are restored to their previous position.
 * <br>
 * (Anne Glerum, 2020/12/18)
 *
 * <li> Fixed: The current viscous stress that is compared to the yield stress
 * in the visco_plastic material model does now include the effects
 * of viscous strain weakening.
 * <br>
 * (Anne Glerum, 2020/12/15)
 *
 * <li> New: ASPECT's ViscoPlastic rheology has been migrated from the
 * ViscoPlastic material model into its own rheology module. It is therefore
 * accessible by not only the ViscoPlastic material model, but by other
 * user-generated material models as well, without needing a large amount
 * of code duplication.
 * <br>
 * (Bob Myhill, 2020/11/30)
 *
 * <li> Fixed: The Multicomponent Incompressible equation
 * of state (EOS) module now uses the correct reference
 * temperature for density when adiabatic heating is used.
 * <br>
 * (John Naliboff and Anne Glerum, 2020/11/28)
 *
 * <li> New: The Drucker Prager rheology module now has an option
 * to include a plastic damper, which acts to stabilize the
 * plasticity formulation. At sufficient resolutions for a given
 * plastic damper viscosity, the plastic shear band characteristics
 * will be resolution independent.
 * <br>
 * (John Naliboff and Cedric Thieulot, 2020/11/19)
 *
 * <li> Fixed: The visualization postprocessors that output stresses now use the
 * correct sign for the shear stress term (under the convention of positive
 * compressive stress), i.e. -2*eta*strain_rate+pressure*I instead of
 * 2*eta*strain_rate+pressure*I.
 * <br>
 * (Anne Glerum, Sungho Lee, 2020/11/19)
 *
 * <li> New: The CompositeViscoPlastic rheology module in ASPECT now provides the
 * option to include a damped Drucker-Prager plastic element.
 * <br>
 * (Bob Myhill, 2020/11/13)
 *
 * <li> New: Particles now supported for a spherical shell geometry with periodicity
 * in polar angle direction for a 2D quarter shell (90 degree opening angle).
 * <br>
 * (Kiran Chotalia, Rene Gassmoeller, 2020/11/10)
 *
 * <li> New: There is now a new property in the depth average postprocessor
 * that averages the mass of a compositional field (rather than its volume,
 * as the "composition" property).
 * <br>
 * (Juliane Dannberg, 2020/11/09)
 *
 * <li> New: ASPECT now includes a CompositeViscoPlastic rheology module.
 * This rheology is an isostress composite of diffusion, dislocation and
 * Peierls creep rheologies. The composite viscosity and partial strain rates
 * are calculated using an internal Newton-Raphson iterative scheme.
 * For improved stability, the composite rheology is arranged in parallel with a
 * constant viscosity strain rate limiter and the whole package is then arranged in
 * series with a constant high viscosity limiter. This ensures that the viscosity
 * function is continuous, differentiable and bounded. The new rheology can be used
 * with multiple compositional fields.
 * <br>
 * (Bob Myhill, 2020/11/05)
 *
 * <li> New: There is a new visualization postprocessor called 'principal stress',
 * which outputs the principal stress values and directions at every point in the
 * model. The postprocessor can either operate on the full stress tensor
 * (including pressure) or on the trace-free deviatoric stress tensor.
 * <br>
 * (Rene Gassmoeller, 2020/10/30)
 *
 * <li> New: Added the functionality to compute averages in user defined depth
 * layers (e.g. lithosphere, asthenosphere, transition zone, lower mantle)
 * to the depth average postprocessor and the lateral averaging plugin.
 * <br>
 * (Rene Gassmoeller, 2020/10/22)
 *
 * <li> Fix: Fixed heatflow boundaries for the Newton solver and refactor the
 * code so it is more maintainable in the future.
 * <br>
 * (Menno Fraters, 2020/10/21)
 *
 * <li> New: The 'spherical shell' geometry model now supports periodic boundary conditions
 * in polar angle direction for a 2D quarter shell (90 degree opening angle).
 * <br>
 * (Kiran Chotalia, Timo Heister, Rene Gassmoeller, 2020/10/08)
 *
 * <li> New: Added a rheology model that computes viscosity as depth-dependent values
 * taken from an input ascii data file. The existing depth-dependent material
 * model is modified to utilize the new rheology.
 * <br>
 * (Arushi Saxena, 2020/09/21)
 *
 * <li> New: The Peierls creep rheology can now accept different
 * parameter values for distinct phases in each compositional field.
 * This unifies the functionality of the Peierls creep rheology
 * with that of diffusion creep and dislocation creep.
 * <br>
 * (Bob Myhill, 2020/09/04)
 *
 * <li> New: Added a new benchmark that examines shortening of a
 * visco-plastic or viscoelastic-plastic block in the absence of
 * gravity. The benchmark is modified from Exercise 13.2 in
 * Gerya 2019 (Introduction to Numerical Geodynamic Modeling).
 * <br>
 * (John Naliboff and Cedric Thieulot, 2020/09/03)
 *
 * <li> New: Added a visualization object in postprocessing that computes the velocity residual at the top boundary of the model domain. The residual is computed between the modeled and the input velocities that can either be in an ascii data or gplates model format.
 * <br>
 * (Arushi Saxena, 2020/09/02)
 *
 * <li> New: The rheology module for Peierls creep now includes a formulation
 * to compute the exact Peierls viscosity, using an internal Newton-Raphson
 * iterative scheme. The module now also contains a function to calculate the
 * strain rate and derivative as a function of stress
 * for both the exact and approximate formulations.
 * <br>
 * (Bob Myhill, John Naliboff and Magali Billen, 2020/08/29)
 *
 * <li> New: The current time step size (dt) is now printed to the screen.
 * <br>
 * (Timo Heister, 2020/08/28)
 *
 * <li> Changed: The steinberger material model now performs
 * consistent averaging of material properties when multiple
 * material lookups are used. Specifically, the density,
 * thermal expansivity, compressibility and phase volume fractions
 * are averaged according to volume fractions, the specific heat
 * is averaged according to mass fractions and the
 * seismic velocities are Voigt-Reuss-Hill averaged.
 * <br>
 * (Bob Myhill, 2020/08/28)
 *
 * <li> New: ASPECT now has a new, reproducible logo.
 * <br>
 * (Rene Gassmoeller, 2020/08/20)
 *
 * <li> New: A new particle interpolator based upon quadratic least squares has
 * been added.
 * <br>
 * (Mack Gregory, 2020/08/19)
 *
 * <li> New: A new class TimeStepping::Manager to control time stepping with a
 * plugin architecture has been added.
 * <br>
 * (Timo Heister, 2020/08/14)
 *
 * <li> New: The PerpleX Lookup module can now accept phase volume
 * fractions, for use in the material models or simply to visualize
 * after the model runs. A python file is also provided in the
 * contributions folder to allow users to easily create input files
 * in the correct format.
 * <br>
 * (Bob Myhill, 2020/08/14)
 *
 * <li> New: The diffusion and dislocation creep rheology modules
 * now have a new function to calculate the strain rate and
 * first stress-derivative of the strain rate to facilitate
 * the development of composite rheologies.
 * <br>
 * (Bob Myhill, John Naliboff, Cedric Theulot, 2020/08/14)
 *
 * <li> New: Added a new input parameter "Stabilization time scale factor" to Rheology::Elasticity.
 *
 * This variable is a stabilization factor for the elastic stresses that influences how fast
 *  elastic stresses adjust to deformation.
 *
 * <br>
 * (Daniel Douglas, Rene Gassmöller, Esther Heckenbach 2020/08/13)
 *
 * <li> Improved: Data loading in AsciiDataLookup has been refactored to be more
 * robust (errors in incorrect coordinate values are now being reported when
 * detected) and data can be set directly without going through a text file by
 * calling reinit() directly. This adds the possibility for other file formats in
 * the future.
 * <br>
 * (Timo Heister, 2020/08/12)
 *
 * <li> New: Added a new particle property (ElasticStress) to store and advect
 * elastic stresses. This particle property allows replacing compositional fields
 * with particles in any model using viscoelasticity.
 * <br>
 * (John Naliboff, 2020/08/12)
 *
 * <li> New: The 'mesh deformation' plugins can now prescribe an initial
 * deformation to the mesh that acts as a starting point for deformations during
 * the model run. This can be used to prescribe an initial topography and will
 * replace the less general 'initial topography' plugin system in a future
 * version.
 * <br>
 * (Rene Gassmoeller, 2020/08/12)
 *
 * <li> New: Added calculation for temperature-dependent strain healing in the strain
 * dependent rheology module.
 * <br>
 * (Erin Heilman, 2020/08/12)
 *
 * <li> Changed: The 'stress', 'shear stress', and 'strain rate tensor'
 * visualization postprocessors now output their information in the form
 * of actual tensors that (with deal.II versions post-9.2) can be
 * visualized as tensors and not just as individual scalar
 * components. This also implies that now *all* components of these
 * tensors, not just the ones on and above the diagonal, are output, even
 * though these tensors are all symmetric.
 * <br>
 * (Wolfgang Bangerth, 2020/08/11)
 *
 * <li> New: A rheology module for Peierls creep has been added
 * (namespace MaterialModel::Rheology::PeierlsCreep). This module
 * is used to define parameters for Peierls creep and compute an approximate
 * viscosity, which is integrated in the visco_plastic material model.
 * <br>
 * (John Naliboff and Magali Billen, 2020/08/10)
 *
 * <li> New: The gravity model 'function' does now evaluate the time at each timestep and the function expression can therefore now be time dependent.
 * <br>
 * (Jacky Austermann, 2020/08/10)
 *
 * <li> Fixed: The Newton solver schemes now correctly postprocess
 * nonlinear iterations.
 * <br>
 * (Juliane Dannberg, 2020/08/07)
 *
 * <li> New: Added new rheology module, namespace MaterialModel::Rheology::FrankKamenetskii.
 * This new rheology computes the temperature dependent Frank Kamenetskii viscosity approximation.
 * <br>
 * (Erin Heilman, 2020/08/07)
 *
 * <li> New: Change the implementations in MaterialModel::Rheology to allow computation of viscosity on phases.
 * Also add optional entry to MaterialModel::ViscoPlastic::calculate_isostrain_viscosities in order to pass in values of phase functions.
 * <br>
 * (Haoyuan Li, 2020/08/06) *
 * <li> Fixed: The latent heat and latent heat melt material models now
 * use the adiabatic temperature as the reference temperature for the
 * viscosity if adiabatic heating is switched on (and hence, the
 * geotherm is an adiabat rather than a constant temperature).
 * <br>
 * (Juliane Dannberg, 2020/08/05)
 *
 * <li> New: There is now a signal called post_nonlinear_solver that is triggered
 * after the nonlinear solver has finished in a timestep and can be used to query
 * information like the number of nonlinear iterations, the initial and final
 * nonlinear residuals, and whether the solver succeeded.
 * <br>
 * (Rene Gassmoeller, 2020/08/04)
 *
 * <li> New: ASPECT now supports the creation of visualization
 * postprocessors that only output data on the surface of a model. An
 * example is the "surface stress" visualization postprocessor.
 * <br>
 * (Wolfgang Bangerth, 2020/08/03)
 *
 * <li> New: Additional compiler flags can be set using the CMake variable
 * ASPECT_ADDITIONAL_CXX_FLAGS.
 * <br>
 * (Timo Heister, 2020/08/02)
 *
 * <li> New: ASPECT now requires deal.II 9.2.0 or newer.
 * <br>
 * (Timo Heister, 2020/08/01)
 *
 * <li> New: ASPECT can now optionally link to the NetCDF library.
 * <br>
 * (Timo Heister, 2020/07/27)
 *
 * <li> New: We now support the CMake option ASPECT_UNITY_BUILD (ON by default)
 * and a rewritten ASPECT_PRECOMPILE_HEADERS option (ON by default) if you use
 * CMake 3.16 or newer. Both options combined can speed up compilation by more
 * than 3x.
 * <br>
 * (Rene Gassmoeller, Timo Heister, 2020/07/06)
 *
 * <li> New: Added a new function PhaseUtilities::phase_average_value to the namespace MaterialModel::MaterialUtilities.
 * This new function averages parameters(e.g density, viscosity) across all phases for a composition.
 * <br>
 * (Haoyuan Li, 2020/07/01)
 *
 * <li> New: There is now a mesh deformation plugin "diffusion" that can be
 * used to diffuse away surface topography in box geometry models. The
 * plugin can be used alone or in combination with the free surface
 * to approximate surface processes and increase numerical stability
 * by limiting mesh distortion.
 * <br>
 * (Anne Glerum, 2020/06/19)
 *
 * </ol>
 */
