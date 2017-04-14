/**
 * @page changes_current Changes after the latest release (v1.5.0)
 *
 * <p> This is the list of changes made after the release of Aspect version
 * 1.5.0. All entries are signed with the names of the author. </p>
 *
 *
 * <ol>
 *
 * <li> Changed: The global statistics information (time, timestep, degrees of
 * freedom, solver iterations), is now written at one point in a new
 * postprocessor 'global_statistics' instead of wherever they happened
 * to be available. This allows to control the order and format of the
 * columns in the statistics file, and fixes a bug in the statistics file,
 * in which data from later timesteps was written into lines
 * of previous timesteps in case of nonlinear solver schemes. Additionally,
 * the postprocessor by default now outputs only a single line per 
 * timestep also for nonlinear solver schemes. To restore the previous
 * behavior (one line per nonlinear iteration) there is a new input
 * parameter "Postprocess/Global statistics/Write statistics for all nonlinear
 * iterations".
 * <br>
 * (Juliane Dannberg, Rene Gassmoeller, 2017/04/14)
 *
 * <li> New: aspect now provides as script that -- when used together with
 * the deal.II parameter GUI program -- allows for a graphical creation and 
 * modification of input parameter files. All available parameters are listed,
 * including their documentation, type, allowed range, and default value. The
 * parameter file written by the GUI will only contain values that are 
 * different from the default values, to keep the file easily readable.
 * <br>
 * (Timo Heister, Rene Gassmoeller, Juliane Dannberg, 2017/03/28)
 *
 * <li> New: aspect now supports the --output-xml flag to generate .xml parameter
 * files that can be edited using the deal.II parameter GUI.
 * <br>
 * (Timo Heister, 2017/03/24)
 *
 * <li> New: Compressible computations now work with discontinuous pressure
 * elements.
 * <br>
 * (Timo Heister, 2017/03/23)
 *
 * </ol>
 */
