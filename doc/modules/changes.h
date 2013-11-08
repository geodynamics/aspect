/**
 * @page changes_after_0_3 Changes after Version 0.3

<p>
This is the list of changes made after the release of
Aspect version 0.3.
All entries are signed with the names of the author.
</p>

<ol>
  <li>New: It is now possible to prescribe the velocity only for certain
  components in the 'Prescribed velocity boundary indicators' parameter.
  <br>
  (Timo Heister 2013/11/08)

  <li>New: the "iterated Stokes" nonlinear solver will now stop iterating
  if the residual is smaller than the new "Nonlinear solver tolerance".
  <br>
  (Timo Heister 2013/11/02)

  <li>New: add a visual postprocessor that outputs the artificial
  viscosity parameter for the temperature equation on each cell.
  <br>
  (Timo Heister 2013/10/28)

  <li>Fixed: moved particle generation to a class, changed particle 
  integration and generation to be factory patterned classes. There
  should be no effect on the user but this will allow for easier
  extension of particle functionality in the future.
  <br>
  (Eric Heien 2013/10/14)

  <li>New: Aspect now not only generates a <code>solution-NNNNN.visit</code>
  file for each time step but also a global <code>solution.visit</code> file
  that Visit can use to visualize the entire time dependent solution. (Both
  of these work with versions of Visit that support this, including Visit 2.5.0.
  Unfortunately, versions of Visit between 2.5.1 and the version current at
  the time of writing this, 2.6.3, have a bug that prevents this.)
  <br>
  (Wolfgang Bangerth 2013/10/08)

  <li>Fixed: Performance of matrix assembly has been improved significantly,
  especially in 3d: assembly of the temperature system is up to three
  times faster, assembly of the Stokes system up to 50%.
  <br>
  (Timo Heister, Thomas Geenen, Wolfgang Bangerth 2013/10/08)

  <li>New: HDF5/XDMF will only output mesh data when the mesh changes,
  reducing total data output significantly. XDMF serialization is also
  properly implemented.
  <br>
  (Eric Heien, 2013/09/27)

  <li>New: HDF5 output now uses DataOutFilter to remove redundant
  vertices/values in output data.
  <br>
  (Eric Heien, 2013/09/27)

  <li>New: Aspect now supports periodic domains (a recent development version of deal.II
  is required).
  <br>
  (Ian Rose, Timo Heister, 2013/09/11)

  <li>New: The manual now has a new section discussing ways to make computations
  go faster.
  <br>
  (Wolfgang Bangerth, 2013/09/03)

  <li>Fixed: The assembly of the temperature and compositional linear systems
  in each time step used an unnecessarily large number of quadrature points.
  This is now fixed.
  <br>
  (Thomas Geenen, 2013/08/30)

  <li>Extended: The ability to compute with tracers has been completely
  overhauled. In particular, there is now a cookbook in the manual that
  describes how to use them.
  <br>
  (Wolfgang Bangerth, 2013/08/13)

  <li>Added: ability to terminate simulation after
  specified number of steps. This is implemented as one of the
  terminator modules that can be selected from the input file.
  <br>
  (Ted Studley, 2013/07/03)

  <li>Added: new initial conditions for the Box geometry
  to test different shaped inclusions.
  <br>
  (Ted Studley, 2013/06/20)

  <li>Removed: the parameter "Nonlinear iteration" was
  not used anywhere so it got removed. You might need
  to delete this line from your prm files.
  <br>
  (Timo Heister, 2013/06/20)

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
