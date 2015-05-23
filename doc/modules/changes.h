/**
 * @page changes_current Changes after the latest release (v1.3)
 *
 * <p> This is the list of changes made after the release of Aspect version
 * 1.3. All entries are signed with the names of the author. </p>
 *
 *
 * <ol>
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
 * </ol>
 */
