/**
 * @page changes_current Changes after the latest release (v1.5.0)
 *
 * <p> This is the list of changes made after the release of Aspect version
 * 1.5.0. All entries are signed with the names of the author. </p>
 *
 *
 * <ol>
 *
 * <li> Changed: The visco_plastic material model now uses a single 
 * compositional field for strain weakening, which can be calculated
 * with the integrated strain invariant particle property or an 
 * equivalent compositional field plugin located in the benchmark 
 * buiter_et_al_2008_jgr.
 * <br>
 * (John Naliboff, 2017/04/14)
 *
 * <li> Changed: aspect now outputs only relative nonlinear residuals of
 * nonlinear iterations in the screen output. In addition, the number of
 * the nonlinear iteration is included in the screen output.
 * <br>
 * (Juliane Dannberg, 2017/04/13)
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
