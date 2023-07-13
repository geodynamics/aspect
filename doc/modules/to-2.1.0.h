/**
 * @page changes_between_2.0.0_and_2.1.0 Changes between version 2.0.0 and version 2.1.0
 *
 * <p> This is the list of changes made after the release of ASPECT version
 * 2.0.0 for version 2.1.0. All entries are signed with the names of the author.
 * </p>
 *
 * <ol>
 *
 * <li> New: ASPECT can optionally compile with the Geodynamic World Builder
 * (https://github.com/GeodynamicWorldBuilder/WorldBuilder/). There
 * are new initial conditions plugins for temperature and composition which use the
 * Geodynamic World Builder to construct initial conditions.
 *
 * (Menno Fraters, 2019/04/24)
 *
 * <li> New: There is a new postprocessor 'entropy viscosity statistics' that
 * computes the maximum and average entropy viscosity stabilization for
 * the temperature equation for diagnostic purposes.
 *
 * (Rene Gassmoeller, 2019/03/22)
 *
 * <li> New: There is now an option to output visualization data as higher order
 * polynomials. This is an improvement in accuracy and requires less disk space
 * than the 'Interpolate output' option that was available before (the new option
 * requires the old option to be set, which it is by default). However the new
 * output can only be read by ParaView version 5.5 and newer and is therefore
 * disabled by default.
 * <br>
 * (Rene Gassmoeller, 2019/03/13)
 *
 * <li> Changed: The default option of the 'Interpolate output' input parameter of the
 * 'Visualization' postprocessor was changed to 'true'. This will by default
 * increase the output resolution compared to the resolution of the internal mesh
 * by a factor of 2 in each dimension. This counteracts the effect that
 * visualization programs usually interpolate linearly between neighboring
 * vertices, while ASPECT's solution usually consists of higher order polynomials.
 * The output now resembles the internal solution more closely at the cost of
 * additional disk space.
 * <br>
 * (Rene Gassmoeller, 2019/03/13)
 *
 * <li> New: Added basic support for a volume-of-fluid interface tracking advection
 * method in 2D incompressible box models. The VoF method is an efficient method
 * to track a distinct compositional field without artificial diffusion.
 *
 * (Jonathan Robey, 2019/01/28)
 *
 * <li> New: ASPECT now outputs a file named original.prm in the output directory with the exact content of the parameter it got started with.
 * <br>
 * (Timo Heister, 2019/01/17)
 *
 * <li> Fixed: If using an initial topography for box geometries with a non-default
 * box origin the topography was computed incorrectly. This is fixed now.
 * <br>
 * (Rene Gassmoeller, Juliane Dannberg, Eva Bredow, 2019/01/08)
 *
 * <li> New: ASPECT now calculates gravity anomalies in addition to the geoid. This can be outputted in ascii format during postprocessing.
 * <br>
 * (Jacky Austermann, 2018/11/15)
 *
 * <li> Changed: If the 'prescribed field' compositional field method is used
 * for a field, this field will now be solved in the order it is listed
 * in in the input file (as is done for all compositional fields) instead
 * of being updated after all of the other fields are solved (the
 * previous behavior).
 * <br>
 * (Juliane Dannberg, 2018/10/31)
 *
 * <li> New: There is now a new method for compositional fields that allows
 * to prescribe them to a value that is computed in the material model
 * and then, in a second step, to diffuse them based on a length scale
 * given as an input parameter. This option can be used by selecting the
 * 'prescribed field with diffusion' compositional field method in the
 * input file, and allows it to smooth quantities that are computed in
 * the material model.
 * <br>
 * (Juliane Dannberg, 2018/10/29)
 *
 * <li> Changed: The heat flux through the boundary cells that is computed
 * in several postprocessors is now computed using the consistent boundary
 * flux method as described in
 *
 * Gresho, et al. (1987). The consistent Galerkin FEM for computing derived
 * boundary quantities in thermal and or fluids problems. International Journal
 * for Numerical Methods in Fluids, 7(4), 371-394.
 *
 * This leads to significantly different (and more accurate) heat flux
 * output. In particular the following changes are expected:
 *   - Heat flux through reflecting boundaries is exactly 0,
 *   - Heat flux through boundaries with Neumann boundary conditions is
 *     computed correctly (previously ignored),
 *   - Heat flux through boundaries with Dirichlet boundary conditions is
 *     much more accurate than before,
 *   - Heat flux due to artificial stabilization terms is correctly included (fixing
 *     apparent heat flux imbalances for steady-state models in spherical geometry),
 *   - Heat flux now includes advective heat flux (i.e. the term velocity times
 *     temperature times density times specific heat capacity), not purely
 *     conductive heat flux.
 *
 * Note that all of these changes purely affect the postprocessing. The
 * accuracy of the solution has not changed.
 *
 * Benchmark and test results have been updated and therefore now differ from
 * publications that described their initial results.
 * <br>
 * (Rene Gassmoeller, 2018/10/18)
 *
 * <li> Fixed: The 'named additional outputs' postprocessor used to only support
 * a single named output from the material model due to an indexing bug. This
 * is fixed now.
 * <br>
 * (Rene Gassmoeller, 2018/10/04)
 *
 * <li> New: Compositional fields can now be prescribed to a value that is
 * computed in the material model as an additional output at every time
 * step and then interpolated and copied into the solution vector (as
 * opposed to advecting the compositional field). This option can be
 * used by selecting the 'prescribed field' compositional field method
 * in the input file.
 * <br>
 * (Juliane Dannberg, 2018/10/01)
 *
 * <li> Improved: The artificial diffusion term that is added in the entropy viscosity
 * method to the temperature and composition equations is now computed as the
 * maximum of the physical diffusion and entropy viscosity instead of the sum.
 * This is reasonable, as the entropy viscosity is the smallest diffusion that is
 * necessary to stabilize the equation, and the presence of physical diffusion
 * stabilizes the equation exactly as an artificial diffusion would. This change
 * likely changes all benchmark results, but leads to more accurate, less
 * diffusive solutions. We verified the new method in selected benchmarks, and
 * updated test results but did not rerun and update all of the contained
 * benchmark cases.
 * <br>
 * (Wolfgang Bangerth, Rene Gassmoeller, 2018/09/26)
 *
 * <li> Fixed: In case of combining reference profile formulations of the equations
 * (ALA, Boussinesq) with discontinuous temperature elements the face terms
 * in the assembly were incorrectly assembled with the full, instead of the
 * reference density. The volume terms were correct. The result were temperature
 * overshoots and possible crashes if the face terms dominate and create densities
 * smaller than zero. This is fixed now, the face terms are now correctly computed
 * with the reference density.
 * <br>
 * (Rene Gassmoeller, 2018/09/15)
 *
 * <li> Fixed: The global entropy variation that is used in the artificial
 * viscosity computation was always computed with a quadrature of the
 * polynomial degree of the temperature (even for compositional fields).
 * This is fixed now, and causes minor changes in the stabilization for
 * compositional fields if they use a different element order than the
 * temperature (usually not the case).
 * <br>
 * (Rene Gassmoeller, 2018/08/31)
 *
 * <li> New: There is now a 'material statistics' postprocessor that computes the
 * average density, average viscosity, and total mass in a
 * model.
 * <br>
 * (Rene Gassmoeller, 2018/08/22)
 *
 * <li> New: ASPECT has two visualization postprocessors which calculate and output the grain lag angle and the infinite strain axis (ISA) rotation timescale, respectively.
 * These two quantities can be used to calculate the grain orientation lag parameter of Kaminski and Ribe (G3, 2002)
 * <br>
 * (Hannah Mark, 2018/08/12)
 *
 * <li> Fixed: Fixed a bug in the 'Named additional outputs' postprocessor, which
 * would cause a segmentation fault if more than one named additional material
 * output object was active in the same model.
 * <br>
 * (Rene Gassmoeller, 2018/08/06)
 *
 * <li> Fixed: Fixed a bug in the 'Steinberger' material model that caused ASPECT
 * to use wrong compressibilities when only a single compositional field, and
 * a single material data file were prescribed. The correct behavior would have
 * been to ignore the compositional field and use the data file as background
 * material, but instead the compositional field was interpreted as volume
 * fraction and multiplied with the compressibility calculated from the
 * data file. This is fixed now.
 * <br>
 * (Rene Gassmoeller, 2018/07/29)
 *
 * <li> New: ASPECT now outputs a dynamically generated URL based on used features to ask people to cite appropriate papers.
 * <br>
 * (Timo Heister, 2018/07/25)
 *
 * <li> New: Added a postprocessor for the temperature anomaly,
 * defined as temperature - (depth average of temperature).
 * <br>
 * (Max Rudolph, 2018/07/11)
 *
 * <li> New: ASPECT now has a Patch on S40RTS initial temperature model which combines S40RTS and
 * an ascii grid containing an upper mantle tomography. Data is input in Vs and converted to
 * temperature within the function (using the same method as S40RTS_perturbation).
 * <br>
 * (Sophie Coulson, 2018/07/09)
 *
 * <li> New: Analytical benchmark for 2D Cartesian incompressible Stokes flow
 * following the FEM book by Donea and Huerta.
 * <br>
 * (Cedric Thieulot, 2018/07/02)
 *
 *
 *
 * <li> New: The benchmarks now include a viscoplastic thermal convection benchmark as
 * described and performed in Tosi et al. (G3 2015).
 * <br>
 * (Anne Glerum, 2018/06/27)
 *
 * <li> Improved: The assembly speed for the Newton solver has been significantly improved
 * by caching various variables.
 * <br>
 * (Menno Fraters and Rene Gassmoeller, 2018/06/26)
 *
 * <li> New: Created a new plugin system that allows it to prescribe a fixed
 * heat flux (instead of prescribing the temperature) at the model boundaries.
 * <br>
 * (Juliane Dannberg and Rene Gassmoeller, 2018/06/25)
 *
 * <li> Fixed: The default latitude and longitude limits for the 'uniform radial' particle generator were set to (0, pi) even though they are measured in degrees.
 * The new default limits correspond to the entire surface of the sphere, and the manual clarifies that latitude is in (0, 180) and not (-90, 90).
 * <br>
 * (Bart Niday, 2018/06/25)
 *
 * <li> New: ASPECT can now call PerpleX to calculate material properties, phase
 * amounts and compositions on-the-fly. A simple material model (PerpleX lookup)
 * is provided which uses PerpleX to provide densities, compressibilities,
 * thermal expansivities and heat capacities on a cell-by-cell basis. This
 * model is provided as a proof-of-concept; more efficient procedures
 * will be required for production runs.
 * <br>
 * (Bob Myhill, 2018/06/25)
 *
 * <li> New: We now have unit_tests in the unit_test/ folder using the Catch testing
 * framework. They should be used instead of writing a test in tests/ if they
 * don't require an ASPECT run.
 * <br>
 * (Timo Heister, 2018/06/24)
 *
 * <li> New: The ascii data and function boundary velocity plugins now allow velocities to be specified along spherical (up, east, north) unit vectors.
 * <br>
 * (Bart Niday, 2018/06/24)
 * <br>
 * New: Added a visualization plugin that directly outputs the strain rate tensor.
 * <br>
 * (Bart Niday, 2018/06/24)
 *
 * <li> New: ASPECT can now read in a depth-dependent initial temperature from ASCII data. This is intended to provide a background temperature profile on which to add perturbations for composite initial temperature models.
 * <br>
 * (Marie Kajan, 2018/06/24)
 *
 * <li> New: ASPECT can now read in a depth-dependent vs to density conversion file, which
 * can be used for S40RTS and SAVANI.
 * <br>
 * (Jacky Austermann, 2018/06/24)
 *
 * <li> New: Function GeometryModel::Interface::cartesian_to_other_coordinates converts between defined coordinate systems. This might break some plugins, which can be fixed by including aspect/geometry_model/interface.h
 * <br>
 * (Bart Niday, 2018/06/23)
 *
 * <li> Changed: We now only set time step number and timestep sizes after
 * initialization is finished.
 *
 * Right now, there is no way to find out, for example in the hook that
 * is called from the function that builds constraints, whether we
 * are already in the process of time stepping, or only in the
 * initialization phase. So move the initialization of the time
 * step number and time step variables to the very end of the
 * initialization phase. Before this, these variables have invalid
 * values.
 * <br>
 * (Wolfgang Bangerth, 2018/06/20)
 *
 * <li> New: There is now a visualization postprocessor that outputs
 * the compaction length, the characteristic length scale of
 * melt transport.
 * <br>
 * (Joe Schools, 2018/06/19)
 *
 * <li> Fixed: When gnuplot was used as the output format for the 'depth average'
 * postprocessor plugin, the first two columns were labeled 'x' 'y', even though
 * 'depth' and 'time' are always the first two variables. This is now fixed
 * provided deal.II 9.0 or newer is used.
 * <br>
 * (Bob Myhill, 2018/06/19)
 *
 * <li> New: Compositional fields can now be advected with the melt velocity
 * (as opposed to the solid velocity), and this option can be used by
 * selecting the 'melt field' advection method in the input file.
 * <br>
 * (Juliane Dannberg, 2018/06/19)
 *
 * <li> New: ASPECT now does some basic checks that essential parameters (such
 * as the number of compositional fields) are the same before creating a
 * checkpoint and after restarting from it.
 * <br>
 * (Wolfgang Bangerth, 2018/06/18)
 *
 * <li> Fixed: The 'compositional heating' heating plugin had a parameter
 * 'Use compositional field for heat production averaging' that was used
 * inconsistently with its description. Its first entry did not correspond
 * to the background field, but to the first compositional field, and the
 * last value was ignored. This is fixed now, the first entry is used for
 * the background field, and all following values determine whether to include
 * the compositional fields.
 * <br>
 * (Rene Gassmoeller, 2018/05/21)
 *
 * <li> Fixed: The nullspace removal had a minor bug that caused a small
 * remaining rotation in three dimensional models even if the net rotation
 * or angular momentum should have been removed. The bug was fixed, and
 * a new postprocessor called 'rotation statistics' was added that can be
 * used to calculate net rotation and angular momentum.
 * <br>
 * (Rene Gassmoeller, Shangxin Liu, 2018/05/03)
 *
 * <li> New: ASPECT now stores the steps for the continuous integration setup
 * inside the repository in form of a Jenkins pipeline file. This file
 * documents each step in our testing procedure, and allows to easily setup
 * new Jenkins server that reproduce the official ASPECT tester.
 * <br>
 * (Rene Gassmoeller, Timo Heister, Tyler Esser, 2018/05/02)
 *
 * </ol>
 */
