/**
 * @page changes_after_1_0 Changes after Version 1.0
 *
 * <p> This is the list of changes made after the release of Aspect version
 * 1.0. All entries are signed with the names of the author. </p>
 *
 * <ol>
 * <li>New: There is now a mesh refinement criterium that checks the mesh
 * always has a minimum refinement level determined by a function given in the
 * input file. 
 * <br>
 * (Juliane Dannberg, 2014/05/19)
 *
 * <li>New: There is a new plugin architecture for heating models that
 * allows for plugins defining the radiogenic heating rate in dependency
 * of temperature, pressure, composition, position and time. This introduces
 * an incompatibility of input files to older versions, but the old functionality
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
