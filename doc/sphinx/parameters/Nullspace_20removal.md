(parameters:Nullspace_20removal)=
# Nullspace removal


## **Subsection:** Nullspace removal


(parameters:Nullspace_20removal/Remove_20nullspace)=
### __Parameter name:__ Remove nullspace
**Default value:**

**Pattern:** [MultipleSelection net rotation|angular momentum|net surface rotation|net translation|linear momentum|net x translation|net y translation|net z translation|linear x momentum|linear y momentum|linear z momentum ]

**Documentation:** Choose none, one or several from

\begin{itemize} \item net rotation \item angular momentum \item net translation \item net surface rotation\item linear momentum \item net x translation \item net y translation \item net z translation \item linear x momentum \item linear y momentum \item linear z momentum \end{itemize}

These are a selection of operations to remove certain parts of the nullspace from the velocity after solving. For some geometries and certain boundary conditions the velocity field is not uniquely determined but contains free translations and/or rotations. Depending on what you specify here, these non-determined modes will be removed from the velocity field at the end of the Stokes solve step.


The &ldquo;angular momentum&rdquo; option removes a rotation such that the net angular momentum is zero. The &ldquo;linear * momentum&rdquo; options remove translations such that the net momentum in the relevant direction is zero.  The &ldquo;net rotation&rdquo; option removes the net rotation of the whole domain, the &ldquo;net surface rotation&rdquo; option removes the net rotation of the top surface, and the &ldquo;net * translation&rdquo; options remove the net translations in the relevant directions.  For most problems there should not be a significant difference between the momentum and rotation/translation versions of nullspace removal, although the momentum versions are more physically motivated. They are equivalent for constant density simulations, and approximately equivalent when the density variations are small.

Note that while more than one operation can be selected it only makes sense to pick one rotational and one translational operation.
