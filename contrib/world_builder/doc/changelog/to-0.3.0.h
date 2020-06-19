/**
 * @page changes_between_0.2.0_and_0.3.0 Changes between version 0.2.0 and version 0.3.0
 * 
 * This is the list of changes made after the release of ASPECT version
 * 0.2.0 and before 0.3.0.
 * 
 * 
 * <ol>
 * <li> Change: Fortran compiler is now optional
 *
 * <li> Change: The world builder no longer requires or uses OpenMP. It was not used in the core 
 * library, but used to speed up the Visualization program. That uses its own threads now. 
 * 
 * <li> New: A Python interface has been added. It is now possible to use the World Builder from a
 * Python program. Integration with numpy has not yet been added.
 * 
 * <li> New: An experimental grains interface containing grain sizes and orientations has been added.
 * 
 * <li> Change: The usage of cmake has been moderized allowing for more flexibility
 * 
 * <li> Fixed: The code has been refactored and cleaned exstinively and new compiler warnings have been 
 * enabled and fixed.
 * </ol>
 * 
 */