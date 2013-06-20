/**
 * @page changes_after_0_3 Changes after Version 0.3

<p>
This is the list of changes made after the release of
Aspect version 0.3.
All entries are signed with the names of the author.
</p>

<ol>
  <li>New: Nonlinear solver scheme = 'Stokes only' only
  solves the Stokes system and ignores temperature and
  compositional fields.
  <br>
  (Timo Heister, 2013/06/20)

  <li>New: If the underlying deal.II version supports this,
  then VTK or VTU files now contain information about the time
  and time step number corresponding to each file, and this
  is then displayed when using Visit as the visualization
  program.
  <br>
  (Wolfgang Bangerth, 2013/06/16)

  <li>New: In order to implement extensions, in particular
  new plugins for material models, geometries, etc, it used to
  be necessary to put the new files into the Aspect source
  directories and re-compile all of Aspect. This is now no
  longer necessary: You can just compile your additional
  plugins into a shared library and tell Aspect via the
  parameter file to load this shared library at start-up.
  Details on this process are provided in the manual in
  the section "How to write a plugin".
  <br>
  (Wolfgang Bangerth, 2013/06/16)
</ol>


*/
