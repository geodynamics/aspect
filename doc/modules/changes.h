/**
 * @page changes_after_0_1 Changes after Version 0.1

<p>
This is the list of changes made after the release of
Aspect version 0.1.
All entries are signed with the names of the author.
</p>



<ol>
<li>
New: The number of space dimensions in which a simulation happens is now
a parameter that is set in the input parameter filer, rather than
statically at compile time as before.
<br>
(Wolfgang Bangerth, 2012/04/03)

<li>
New: A postprocessor module for particles was added. This module creates
a random uniform distribution of particles in the mesh and moves them based
on the velocity field using an RK2 scheme. There is still a bug that occurs
in some setups, once this is fixed this entry will be updated.
<br>
(Eric Heien, 2012/04/02)
</ol>


*/
