(parameters:Postprocess)=
# Postprocess


## **Subsection:** Postprocess


(parameters:Postprocess/List_20of_20postprocessors)=
### __Parameter name:__ List of postprocessors
**Default value:**

**Pattern:** [MultipleSelection Stokes residual|basic statistics|boundary densities|boundary pressures|boundary strain rate residual statistics|boundary velocity residual statistics|command|composition statistics|core statistics|depth average|dynamic topography|entropy viscosity statistics|geoid|global statistics|gravity calculation|heat flux densities|heat flux map|heat flux statistics|heating statistics|load balance statistics|mass flux statistics|material statistics|matrix statistics|maximum depth of field|melt statistics|memory statistics|mobility statistics|particle count statistics|particles|point values|pressure statistics|rotation statistics|spherical velocity statistics|temperature statistics|topography|velocity boundary statistics|velocity statistics|viscous dissipation statistics|visualization|volume of fluid statistics ]

**Documentation:** A comma separated list of postprocessor objects that should be run at the end of each time step. Some of these postprocessors will declare their own parameters which may, for example, include that they will actually do something only every so many time steps or years. Alternatively, the text &lsquo;all&rsquo; indicates that all available postprocessors should be run after each time step.

The following postprocessors are available:

&lsquo;Stokes residual&rsquo;: A postprocessor that outputs the Stokes residuals during the iterative solver algorithm into a file stokes_residuals.txt in the output directory.

&lsquo;basic statistics&rsquo;: A postprocessor that outputs some simplified statistics like the Rayleigh number and other quantities that only make sense in certain model setups. The output is written after completing initial adaptive refinement steps. The postprocessor assumes a point at the surface at the adiabatic surface temperature and pressure is a reasonable reference condition for computing these properties. Furthermore, the Rayleigh number is computed using the model depth (i.e. not the radius of the Earth), as we need a definition that is geometry independent. Take care when comparing these values to published studies and make sure they use the same definitions.

&lsquo;boundary densities&rsquo;: A postprocessor that computes the laterally averaged density at the top and bottom of the domain.

&lsquo;boundary pressures&rsquo;: A postprocessor that computes the laterally averaged pressure at the top and bottom of the domain.

&lsquo;boundary strain rate residual statistics&rsquo;: A postprocessor that computes some statistics about the surface strain rate residual along the top boundary. The residual is the difference between the second invariant of the model strain rate and the second strain rate invariant read from the input data file. Currently, the strain residual statistics, i.e., min, max and the rms magnitude, are computed at the top suface.

&lsquo;boundary velocity residual statistics&rsquo;: A postprocessor that computes some statistics about the velocity residual along the top boundary. The velocity residual is the difference between the model solution velocities and the input velocities (GPlates model or ascii data). Currently, the velocity residual statistics, i.e., min, max and the rms magnitude, is computed at the top suface.

&lsquo;command&rsquo;: A postprocessor that executes a command line process.

&lsquo;composition statistics&rsquo;: A postprocessor that computes some statistics about the compositional fields, if present in this simulation. In particular, it computes maximal and minimal values of each field, as well as the total mass contained in this field as defined by the integral $m_i(t) = \int_\Omega c_i(\mathbf x,t) \; \text{d}x$.

&lsquo;core statistics&rsquo;: A postprocessor that computes some statistics about the core evolution. (Working only with dynamic core boundary temperature plugin)

&lsquo;depth average&rsquo;: A postprocessor that computes depth averaged quantities and writes them into a file <depth_average.ext> in the output directory, where the extension of the file is determined by the output format you select. In addition to the output format, a number of other parameters also influence this postprocessor, and they can be set in the section `Postprocess/Depth average` in the input file.

In the output files, the $x$-value of each data point corresponds to the depth, whereas the $y$-value corresponds to the simulation time. The time is provided in seconds or, if the global &ldquo;Use years in output instead of seconds&rdquo; parameter is set, in years.

&lsquo;dynamic topography&rsquo;: A postprocessor that computes a measure of dynamic topography based on the stress at the surface and bottom. The data is written into text files named &lsquo;dynamic\_topography.NNNNN&rsquo; in the output directory, where NNNNN is the number of the time step.

The exact approach works as follows: At the centers of all cells that sit along the top surface, we evaluate the stress and evaluate the component of it in the direction in which gravity acts. In other words, we compute $\sigma_{rr}={\hat g}^T(2 \eta \varepsilon(\mathbf u)- \frac 13 (\textrm{div}\;\mathbf u)I)\hat g - p_d$ where $\hat g = \mathbf g/\|\mathbf g\|$ is the direction of the gravity vector $\mathbf g$ and $p_d=p-p_a$ is the dynamic pressure computed by subtracting the adiabatic pressure $p_a$ from the total pressure $p$ computed as part of the Stokes solve. From this, the dynamic topography is computed using the formula $h=\frac{\sigma_{rr}}{(\mathbf g \cdot \mathbf n)  \rho}$ where $\rho$ is the density at the cell center. For the bottom surface we chose the convection that positive values are up (out) and negative values are in (down), analogous to the deformation of the upper surface. Note that this implementation takes the direction of gravity into account, which means that reversing the flow in backward advection calculations will not reverse the instantaneous topography because the reverse flow will be divided by the reverse surface gravity.
The file format then consists of lines with Euclidean coordinates followed by the corresponding topography value.

(As a side note, the postprocessor chooses the cell center instead of the center of the cell face at the surface, where we really are interested in the quantity, since this often gives better accuracy. The results should in essence be the same, though.)

&lsquo;entropy viscosity statistics&rsquo;: A postprocessor that computes the maximum and volume averagedentropy viscosity stabilization for the temperature field.

&lsquo;geoid&rsquo;: A postprocessor that computes a representation of the geoid based on the density structure in the mantle, as well as the topography at the surface and core mantle boundary (CMB) if desired. The topography is based on the dynamic topography postprocessor in case of no free surface, and based on the real surface from the geometry model in case of a free surface. The geoid is computed from a spherical harmonic expansion, so the geometry of the domain must be a 3D spherical shell.

&lsquo;global statistics&rsquo;: A postprocessor that outputs all the global statistics information, e.g. the time of the simulation, the timestep number, number of degrees of freedom and solver iterations for each timestep. The postprocessor can output different formats, the first printing one line in the statistics file per nonlinear solver iteration (if a nonlinear solver scheme is selected). The second prints one line per timestep, summing the information about all nonlinear iterations in this line. Note that this postprocessor is always active independent on whether or not it is selected in the parameter file.

&lsquo;gravity calculation&rsquo;: A postprocessor that computes gravity, gravity anomalies, gravity potential and gravity gradients for a set of points (e.g. satellites) in or above the model surface for either a user-defined range of latitudes, longitudes and radius or a list of point coordinates.Spherical coordinates in the output file are radius, colatitude and colongitude. Gravity is here based on the density distribution from the material model (and non adiabatic). This means that the density may come directly from an ascii file. This postprocessor also computes theoretical gravity and its derivatives, which corresponds to the analytical solution of gravity in the same geometry but filled with a reference density. The reference density is also used to determine density anomalies for computing gravity anomalies. Thus one must carefully evaluate the meaning of the gravity anomaly output, because the solution may not reflect the actual gravity anomaly (due to differences in the assumed reference density). On way to guarantee correct gravity anomalies is to subtract gravity of a certain point from the average gravity on the map. Another way is to directly use density anomalies for this postprocessor.The average- minimum- and maximum gravity acceleration and potential are written into the statistics file.

&lsquo;heat flux densities&rsquo;: A postprocessor that computes some statistics about the heat flux density for each boundary id. The heat flux density across each boundary is computed in outward direction, i.e., from the domain to the outside. The heat flux is computed as sum of advective heat flux and conductive heat flux through Neumann boundaries, both computed as integral over the boundary area, and conductive heat flux through Dirichlet boundaries, which is computed using the consistent boundary flux method as described in &ldquo;Gresho, Lee, Sani, Maslanik, Eaton (1987). The consistent Galerkin FEM for computing derived boundary quantities in thermal and or fluids problems. International Journal for Numerical Methods in Fluids, 7(4), 371-394.&rdquo;

Note that the &ldquo;heat flux statistics&rdquo; postprocessor computes the same quantity as the one here, but not divided by the area of the surface. In other words, it computes the *total* heat flux through each boundary.

&lsquo;heat flux map&rsquo;: A postprocessor that computes the heat flux density across each boundary in outward direction, i.e., from the domain to the outside. The heat flux is computed as sum of advective heat flux and conductive heat flux through Neumann boundaries, both computed as integral over the boundary area, and conductive heat flux through Dirichlet boundaries, which is computed using the consistent boundary flux method as described in &ldquo;Gresho, Lee, Sani, Maslanik, Eaton (1987). The consistent Galerkin FEM for computing derived boundary quantities in thermal and or fluids problems. International Journal for Numerical Methods in Fluids, 7(4), 371-394.&rdquo;

&lsquo;heat flux statistics&rsquo;: A postprocessor that computes some statistics about the heat flux density across each boundary in outward direction, i.e., from the domain to the outside. The heat flux is computed as sum of advective heat flux and conductive heat flux through Neumann boundaries, both computed as integral over the boundary area, and conductive heat flux through Dirichlet boundaries, which is computed using the consistent boundary flux method as described in &ldquo;Gresho, Lee, Sani, Maslanik, Eaton (1987). The consistent Galerkin FEM for computing derived boundary quantities in thermal and or fluids problems. International Journal for Numerical Methods in Fluids, 7(4), 371-394.&rdquo;The point-wise heat flux can be obtained from the heat flux map postprocessor, which outputs the heat flux to a file, or the heat flux map visualization postprocessor, which outputs the heat flux for visualization.

As stated, this postprocessor computes the *outbound* heat flux. If you are interested in the opposite direction, for example from the core into the mantle when the domain describes the mantle, then you need to multiply the result by -1.

\note{In geodynamics, the term &ldquo;heat flux&rdquo; is often understood to be the quantity $- k \nabla T$, which is really a heat flux *density*, i.e., a vector-valued field. In contrast to this, the current postprocessor only computes the integrated flux over each part of the boundary. Consequently, the units of the quantity computed here are $W=\frac{J}{s}$.}

The &ldquo;heat flux densities&rdquo; postprocessor computes the same quantity as the one here, but divided by the area of the surface.

&lsquo;heating statistics&rsquo;: A postprocessor that computes some statistics about heating, averaged by volume.

&lsquo;load balance statistics&rsquo;: A postprocessor that computes statistics about the distribution of cells, and if present particles across subdomains. In particular, it computes maximal, average and minimal number of cells across all ranks. If there are particles it also computes the maximal, average, and minimum number of particles across all ranks, and maximal, average, and minimal ratio between local number of particles and local number of cells across all processes. All of these numbers can be useful to assess the load balance between different MPI ranks, as the difference between the mimimal and maximal load should be as small as possible.

&lsquo;mass flux statistics&rsquo;: A postprocessor that computes some statistics about the mass flux across boundaries. For each boundary indicator (see your geometry description for which boundary indicators are used), the mass flux is computed in outward direction, i.e., from the domain to the outside, using the formula $\int_{\Gamma_i} \rho \mathbf v \cdot \mathbf n$ where $\Gamma_i$ is the part of the boundary with indicator $i$, $\rho$ is the density as reported by the material model, $\mathbf v$ is the velocity, and $\mathbf n$ is the outward normal.

As stated, this postprocessor computes the *outbound* mass flux. If you are interested in the opposite direction, for example from the core into the mantle when the domain describes the mantle, then you need to multiply the result by -1.

\note{In geodynamics, the term &ldquo;mass flux&rdquo; is often understood to be the quantity $\rho \mathbf v$, which is really a mass flux *density*, i.e., a vector-valued field. In contrast to this, the current postprocessor only computes the integrated flux over each part of the boundary. Consequently, the units of the quantity computed here are $\frac{kg}{s}$.}

&lsquo;material statistics&rsquo;: A postprocessor that computes some statistics about the material properties. In particular, it computes the volume-averages of the density and viscosity, and the total mass in the model. Specifically, it implements the following formulas in each time step: $\left<\rho\right> = \frac{1}{|\Omega|} \int_\Omega \rho(\mathbf x) \, \text{d}x$, $\left<\eta\right> = \frac{1}{|\Omega|} \int_\Omega \eta(\mathbf x) \, \text{d}x$, $M = \int_\Omega \rho(\mathbf x) \, \text{d}x$, where $|\Omega|$ is the volume of the domain.

&lsquo;matrix statistics&rsquo;: A postprocessor that computes some statistics about the matrices. In particular, it outputs total memory consumption, total non-zero elements, and non-zero elements per block, for system matrix and system preconditioner matrix.

&lsquo;maximum depth of field&rsquo;: A postprocessor that for each compositional field outputs the largest depth at which a quadrature point is found where the field has a value of 0.5 or larger. For fields that do not represent materials, but for example track a certain quantity like strain, this value of 0.5 does not necessarily make sense.

&lsquo;melt statistics&rsquo;: A postprocessor that computes some statistics about the melt (volume) fraction. If the material model does not implement a melt fraction function, the output is set to zero.

&lsquo;memory statistics&rsquo;: A postprocessor that computes some statistics about the memory consumption. In particular, it computes the memory usage of the system matrix, triangulation, p4est, DoFHandler, current constraints, solution vector, and peak virtual memory usage, all in MB. It also outputs the memory usage of the system matrix to the screen.

&lsquo;mobility statistics&rsquo;: A postprocessor that computes some statistics about mobility following Tackley (2000) and Lourenco et al. (2020).

&lsquo;particle count statistics&rsquo;: A postprocessor that computes some statistics about the particle distribution, if present in this simulation. In particular, it computes minimal, average and maximal values of particles per cell in the global domain.

&lsquo;particles&rsquo;: A Postprocessor that creates particles that follow the velocity field of the simulation. The particles can be generated and propagated in various ways and they can carry a number of constant or time-varying properties. The postprocessor can write output positions and properties of all particles at chosen intervals, although this is not mandatory. It also allows other parts of the code to query the particles for information.

&lsquo;point values&rsquo;: A postprocessor that evaluates the solution (i.e., velocity, pressure, temperature, and compositional fields along with other fields that are treated as primary variables) at the end of every time step or after a user-specified time interval at a given set of points and then writes this data into the file <point\_values.txt> in the output directory. The points at which the solution should be evaluated are specified in the section `Postprocess/Point values` in the input file.

In the output file, data is organized as (i) time, (ii) the 2 or 3 coordinates of the evaluation points, and (iii) followed by the values of the solution vector at this point. The time is provided in seconds or, if the global &ldquo;Use years in output instead of seconds&rdquo; parameter is set, in years. In the latter case, the velocity is also converted to meters/year, instead of meters/second.

\note{Evaluating the solution of a finite element field at arbitrarily chosen points is an expensive process. Using this postprocessor will only be efficient if the number of evaluation points or output times is relatively small. If you need a very large number of evaluation points, you should consider extracting this information from the visualization program you use to display the output of the &lsquo;visualization&rsquo; postprocessor.}

&lsquo;pressure statistics&rsquo;: A postprocessor that computes some statistics about the pressure field.

&lsquo;rotation statistics&rsquo;: A postprocessor that computes some statistics about the rotational velocity of the model (i.e. integrated net rotation and angular momentum). In 2D we assume the model to be a cross-section through an infinite domain in z direction, with a zero z-velocity. Thus, the z-axis is the only possible rotation axis and both moment of inertia and angular momentum are scalar instead of tensor quantities.

&lsquo;spherical velocity statistics&rsquo;: A postprocessor that computes radial, tangential and total RMS velocity.

&lsquo;temperature statistics&rsquo;: A postprocessor that computes some statistics about the temperature field.

&lsquo;topography&rsquo;: A postprocessor intended for use with a deforming top surface.  After every step it loops over all the vertices on the top surface and determines the maximum and minimum topography relative to a reference datum (initial box height for a box geometry model or initial radius for a sphere/spherical shell geometry model). If &rsquo;Topography.Output to file&rsquo; is set to true, also outputs topography into text files named &lsquo;topography.NNNNN&rsquo; in the output directory, where NNNNN is the number of the time step.
The file format then consists of lines with Euclidean coordinates followed by the corresponding topography value.Topography is printed/written in meters.

&lsquo;velocity boundary statistics&rsquo;: A postprocessor that computes some statistics about the velocity along the boundaries. For each boundary indicator (see your geometry description for which boundary indicators are used), the min and max velocity magnitude is computed.

&lsquo;velocity statistics&rsquo;: A postprocessor that computes some statistics about the velocity field.

&lsquo;viscous dissipation statistics&rsquo;: A postprocessor that outputs the viscous rate of dissipation of energy for each compositional field (where the field has a value of 0.5 or more) as well as over the whole domain. When all the fields represent lithologies and there is no background field, the sum of the individual field&rsquo;s dissipation should equal that over the whole domain. The viscous dissipation is computed as: $\int_{V}\left(\sigma&rsquo; \dot{\epsilon}&rsquo; \right)$, where $\sigma&rsquo;$  is the deviatoric stress and $\dot{\epsilon}&rsquo;$ the deviatoric strain rate.Note then when shear heating is included in the temperature equation, it is better to use the &rsquo;heating statistics&rsquo; postprocessor.

&lsquo;visualization&rsquo;: A postprocessor that takes the solution and writes it into files that can be read by a graphical visualization program. Additional run time parameters are read from the parameter subsection &rsquo;Visualization&rsquo;.

&lsquo;volume of fluid statistics&rsquo;: A postprocessor that computes some statistics about the volume-of-fluid fields.

(parameters:Postprocess/Run_20postprocessors_20on_20nonlinear_20iterations)=
### __Parameter name:__ Run postprocessors on nonlinear iterations
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Whether or not the postprocessors should be executed after each of the nonlinear iterations done within one time step. As this is mainly an option for the purposes of debugging, it is not supported when the &rsquo;Time between graphical output&rsquo; is larger than zero, or when the postprocessor is not intended to be run more than once per timestep.

(parameters:Postprocess/Boundary_20strain_20rate_20residual_20statistics)=
## **Subsection:** Postprocess / Boundary strain rate residual statistics
(parameters:Postprocess/Boundary_20strain_20rate_20residual_20statistics/Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/postprocess/boundary-strain-rate-residual/

**Pattern:** [DirectoryName]

**Documentation:** The name of a directory that contains the ascii data. This path may either be absolute (if starting with a &lsquo;/&rsquo;) or relative to the current directory. The path may also include the special text &lsquo;$ASPECT_SOURCE_DIR&rsquo; which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the &lsquo;data/&rsquo; subdirectory of ASPECT.

(parameters:Postprocess/Boundary_20strain_20rate_20residual_20statistics/Data_20file_20name)=
### __Parameter name:__ Data file name
**Default value:** box_3d_boundary_strain_rate.txt

**Pattern:** [Anything]

**Documentation:** The file name of the input surface strain rate an ascii data. The file has one column in addition to the coordinate columns corresponding to the second invariant of strain rate.

(parameters:Postprocess/Boundary_20strain_20rate_20residual_20statistics/Scale_20factor)=
### __Parameter name:__ Scale factor
**Default value:** 1.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Scalar factor, which is applied to the model data. You might want to use this to scale the input to a reference model.

(parameters:Postprocess/Boundary_20velocity_20residual_20statistics)=
## **Subsection:** Postprocess / Boundary velocity residual statistics
(parameters:Postprocess/Boundary_20velocity_20residual_20statistics/Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/boundary-velocity/gplates/

**Pattern:** [DirectoryName]

**Documentation:** The name of a directory that contains the GPlates model or the ascii data. This path may either be absolute (if starting with a &lsquo;/&rsquo;) or relative to the current directory. The path may also include the special text &lsquo;$ASPECT_SOURCE_DIR&rsquo; which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the &lsquo;data/&rsquo; subdirectory of ASPECT.

(parameters:Postprocess/Boundary_20velocity_20residual_20statistics/Data_20file_20name)=
### __Parameter name:__ Data file name
**Default value:** current_day.gpml

**Pattern:** [Anything]

**Documentation:** The file name of the input velocity as a GPlates model or an ascii data. For the GPlates model, provide file in the same format as described in the &rsquo;gplates&rsquo; boundary velocity plugin. For the ascii data, provide file in the same format as described in  &rsquo;ascii data&rsquo; initial composition plugin.

(parameters:Postprocess/Boundary_20velocity_20residual_20statistics/Scale_20factor)=
### __Parameter name:__ Scale factor
**Default value:** 1.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Scalar factor, which is applied to the model data. You might want to use this to scale the input to a reference model. Another way to use this factor is to convert units of the input files. For instance, if you provide velocities in cm/year set this factor to 0.01.

(parameters:Postprocess/Boundary_20velocity_20residual_20statistics/Use_20ascii_20data)=
### __Parameter name:__ Use ascii data
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Use ascii data files (e.g., GPS) for computing residual velocities instead of GPlates data.

(parameters:Postprocess/Boundary_20velocity_20residual_20statistics/Use_20spherical_20unit_20vectors)=
### __Parameter name:__ Use spherical unit vectors
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Specify velocity as r, phi, and theta components instead of x, y, and z. Positive velocities point up, east, and north (in 3D) or out and clockwise (in 2D). This setting only makes sense for spherical geometries.GPlates data is always interpreted to be in east and north directions and is not affected by this parameter.

(parameters:Postprocess/Command)=
## **Subsection:** Postprocess / Command
(parameters:Postprocess/Command/Command)=
### __Parameter name:__ Command
**Default value:**

**Pattern:** [Anything]

**Documentation:** Command to execute.

(parameters:Postprocess/Command/Run_20on_20all_20processes)=
### __Parameter name:__ Run on all processes
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Whether to run command from all processes (true), or only on process 0 (false).

(parameters:Postprocess/Command/Terminate_20on_20failure)=
### __Parameter name:__ Terminate on failure
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Select whether ASPECT should terminate if the command returns a non-zero exit status.

(parameters:Postprocess/Depth_20average)=
## **Subsection:** Postprocess / Depth average
(parameters:Postprocess/Depth_20average/Depth_20boundaries_20of_20zones)=
### __Parameter name:__ Depth boundaries of zones
**Default value:**

**Pattern:** [List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** The depth boundaries of zones within which we are to compute averages. By default this list is empty and we subdivide the entire domain into equidistant depth zones and compute averages within each of these zones. If this list is not empty it has to contain one more entry than the &rsquo;Number of zones&rsquo; parameter, representing the upper and lower depth boundary of each zone. It is not necessary to cover the whole depth-range (i.e. you can select to only average in a single layer by choosing 2 arbitrary depths as the boundaries of that layer).

(parameters:Postprocess/Depth_20average/List_20of_20output_20variables)=
### __Parameter name:__ List of output variables
**Default value:** all

**Pattern:** [MultipleSelection all|temperature|composition|adiabatic temperature|adiabatic pressure|adiabatic density|adiabatic density derivative|velocity magnitude|sinking velocity|rising velocity|Vs|Vp|log viscosity|viscosity|vertical heat flux|vertical mass flux|composition mass ]

**Documentation:** A comma separated list which specifies which quantities to average in each depth slice. It defaults to averaging all available quantities, but this can be an expensive operation, so you may want to select only a few.

Specifically, the sinking velocity is defined as the scalar product of the velocity and a unit vector in the direction of gravity, if positive (being zero if this product is negative, which would correspond to an upward velocity). The rising velocity is the opposite: the scalar product of the velocity and a unit vector in the direction opposite of gravity, if positive (being zero for downward velocities).

List of options:
all|temperature|composition|adiabatic temperature|adiabatic pressure|adiabatic density|adiabatic density derivative|velocity magnitude|sinking velocity|rising velocity|Vs|Vp|log viscosity|viscosity|vertical heat flux|vertical mass flux|composition mass

(parameters:Postprocess/Depth_20average/Number_20of_20zones)=
### __Parameter name:__ Number of zones
**Default value:** 10

**Pattern:** [Integer range 1...2147483647 (inclusive)]

**Documentation:** The number of zones in depth direction within which we are to compute averages. By default, we subdivide the entire domain into 10 depth zones and compute temperature and other averages within each of these zones. However, if you have a very coarse mesh, it may not make much sense to subdivide the domain into so many zones and you may wish to choose less than this default. It may also make computations slightly faster. On the other hand, if you have an extremely highly resolved mesh, choosing more zones might also make sense.

(parameters:Postprocess/Depth_20average/Output_20format)=
### __Parameter name:__ Output format
**Default value:** gnuplot, txt

**Pattern:** [MultipleSelection none|dx|ucd|gnuplot|povray|eps|gmv|tecplot|tecplot_binary|vtk|vtu|hdf5|svg|deal.II intermediate|txt ]

**Documentation:** A list of formats in which the output shall be produced. The format in which the output is generated also determines the extension of the file into which data is written. The list of possible output formats that can be given here is documented in the appendix of the manual where the current parameter is described. By default the output is written as gnuplot file (for plotting), and as a simple text file.

(parameters:Postprocess/Depth_20average/Time_20between_20graphical_20output)=
### __Parameter name:__ Time between graphical output
**Default value:** 1e8

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The time interval between each generation of graphical output files. A value of zero indicates that output should be generated in each time step. Units: years if the &rsquo;Use years in output instead of seconds&rsquo; parameter is set; seconds otherwise.

(parameters:Postprocess/Dynamic_20core_20statistics)=
## **Subsection:** Postprocess / Dynamic core statistics
(parameters:Postprocess/Dynamic_20core_20statistics/Excess_20entropy_20only)=
### __Parameter name:__ Excess entropy only
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Output the excess entropy only instead the each entropy terms.

(parameters:Postprocess/Dynamic_20topography)=
## **Subsection:** Postprocess / Dynamic topography
(parameters:Postprocess/Dynamic_20topography/Density_20above)=
### __Parameter name:__ Density above
**Default value:** 0.

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Dynamic topography is calculated as the excess or lack of mass that is supported by mantle flow. This value depends on the density of material that is moved up or down, i.e. crustal rock, and the density of the material that is displaced (generally water or air). While the density of crustal rock is part of the material model, this parameter &lsquo;Density above&rsquo; allows the user to specify the density value of material that is displaced above the solid surface. By default this material is assumed to be air, with a density of 0. Units: \si{\kilogram\per\meter\cubed}.

(parameters:Postprocess/Dynamic_20topography/Density_20below)=
### __Parameter name:__ Density below
**Default value:** 9900.

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Dynamic topography is calculated as the excess or lack of mass that is supported by mantle flow. This value depends on the density of material that is moved up or down, i.e. mantle above CMB, and the density of the material that is displaced (generally outer core material). While the density of mantle rock is part of the material model, this parameter &lsquo;Density below&rsquo; allows the user to specify the density value of material that is displaced below the solid surface. By default this material is assumed to be outer core material with a density of 9900. Units: \si{\kilogram\per\meter\cubed}.

(parameters:Postprocess/Dynamic_20topography/Output_20bottom)=
### __Parameter name:__ Output bottom
**Default value:** true

**Pattern:** [Bool]

**Documentation:** Whether to output a file containing the bottom (i.e., CMB) dynamic topography.

(parameters:Postprocess/Dynamic_20topography/Output_20surface)=
### __Parameter name:__ Output surface
**Default value:** true

**Pattern:** [Bool]

**Documentation:** Whether to output a file containing the surface dynamic topography.

(parameters:Postprocess/Geoid)=
## **Subsection:** Postprocess / Geoid
(parameters:Postprocess/Geoid/Also_20output_20the_20gravity_20anomaly)=
### __Parameter name__: Also output the gravity anomaly
**Alias:** [Output gravity anomaly](parameters:Postprocess/Geoid/Output_20gravity_20anomaly)

**Deprecation Status:** false

(parameters:Postprocess/Geoid/Also_20output_20the_20spherical_20harmonic_20coefficients_20of_20CMB_20dynamic_20topography_20contribution)=
### __Parameter name__: Also output the spherical harmonic coefficients of CMB dynamic topography contribution
**Alias:** [Output CMB topography contribution coefficients](parameters:Postprocess/Geoid/Output_20CMB_20topography_20contribution_20coefficients)

**Deprecation Status:** false

(parameters:Postprocess/Geoid/Also_20output_20the_20spherical_20harmonic_20coefficients_20of_20density_20anomaly_20contribution)=
### __Parameter name__: Also output the spherical harmonic coefficients of density anomaly contribution
**Alias:** [Output density anomaly contribution coefficients](parameters:Postprocess/Geoid/Output_20density_20anomaly_20contribution_20coefficients)

**Deprecation Status:** false

(parameters:Postprocess/Geoid/Also_20output_20the_20spherical_20harmonic_20coefficients_20of_20geoid_20anomaly)=
### __Parameter name__: Also output the spherical harmonic coefficients of geoid anomaly
**Alias:** [Output geoid anomaly coefficients](parameters:Postprocess/Geoid/Output_20geoid_20anomaly_20coefficients)

**Deprecation Status:** false

(parameters:Postprocess/Geoid/Also_20output_20the_20spherical_20harmonic_20coefficients_20of_20surface_20dynamic_20topography_20contribution)=
### __Parameter name__: Also output the spherical harmonic coefficients of surface dynamic topography contribution
**Alias:** [Output surface topography contribution coefficients](parameters:Postprocess/Geoid/Output_20surface_20topography_20contribution_20coefficients)

**Deprecation Status:** false

(parameters:Postprocess/Geoid/Density_20above)=
### __Parameter name:__ Density above
**Default value:** 0.

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The density value above the surface boundary.

(parameters:Postprocess/Geoid/Density_20below)=
### __Parameter name:__ Density below
**Default value:** 9900.

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The density value below the CMB boundary.

(parameters:Postprocess/Geoid/Include_20CMB_20topography_20contribution)=
### __Parameter name:__ Include CMB topography contribution
**Default value:** true

**Pattern:** [Bool]

**Documentation:** Option to include the contribution from CMB topography on geoid. The default is true.

(parameters:Postprocess/Geoid/Include_20surface_20topography_20contribution)=
### __Parameter name:__ Include surface topography contribution
**Default value:** true

**Pattern:** [Bool]

**Documentation:** Option to include the contribution from surface topography on geoid. The default is true.

(parameters:Postprocess/Geoid/Include_20the_20contributon_20from_20dynamic_20topography)=
### __Parameter name:__ Include the contributon from dynamic topography
**Default value:** true

**Pattern:** [Bool]

**Documentation:** Option to include the contribution from dynamic topography on geoid. The default is true.

(parameters:Postprocess/Geoid/Maximum_20degree)=
### __Parameter name:__ Maximum degree
**Default value:** 20

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** This parameter can be a random positive integer. However, the value normally should not exceed the maximum degree of the initial perturbed temperature field. For example, if the initial temperature uses S40RTS, the maximum degree should not be larger than 40.

(parameters:Postprocess/Geoid/Minimum_20degree)=
### __Parameter name:__ Minimum degree
**Default value:** 2

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** This parameter normally is set to 2 since the perturbed gravitational potential at degree 1 always vanishes in a reference frame with the planetary center of mass same as the center of figure.

(parameters:Postprocess/Geoid/Output_20CMB_20topography_20contribution_20coefficients)=
### __Parameter name:__ Output CMB topography contribution coefficients
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Option to output the spherical harmonic coefficients of the CMB topography contribution to the maximum degree. The default is false.

(parameters:Postprocess/Geoid/Output_20data_20in_20geographical_20coordinates)=
### __Parameter name:__ Output data in geographical coordinates
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Option to output the geoid anomaly in geographical coordinates (latitude and longitude). The default is false, so postprocess will output the data in geocentric coordinates (x,y,z) as normally.

(parameters:Postprocess/Geoid/Output_20density_20anomaly_20contribution_20coefficients)=
### __Parameter name:__ Output density anomaly contribution coefficients
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Option to output the spherical harmonic coefficients of the density anomaly contribution to the maximum degree. The default is false.

(parameters:Postprocess/Geoid/Output_20geoid_20anomaly_20coefficients)=
### __Parameter name:__ Output geoid anomaly coefficients
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Option to output the spherical harmonic coefficients of the geoid anomaly up to the maximum degree. The default is false, so postprocess will only output the geoid anomaly in grid format.

(parameters:Postprocess/Geoid/Output_20gravity_20anomaly)=
### __Parameter name:__ Output gravity anomaly
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Option to output the free-air gravity anomaly up to the maximum degree. The unit of the output is in SI, hence $m/s^2$ ($1mgal = 10^-5 m/s^2$). The default is false.

(parameters:Postprocess/Geoid/Output_20surface_20topography_20contribution_20coefficients)=
### __Parameter name:__ Output surface topography contribution coefficients
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Option to output the spherical harmonic coefficients of the surface topography contribution to the maximum degree. The default is false.

(parameters:Postprocess/Global_20statistics)=
## **Subsection:** Postprocess / Global statistics
(parameters:Postprocess/Global_20statistics/Write_20statistics_20for_20each_20nonlinear_20iteration)=
### __Parameter name:__ Write statistics for each nonlinear iteration
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Whether to put every nonlinear iteration into a separate line in the statistics file (if true), or to output only one line per time step that contains the total number of iterations of the Stokes and advection linear system solver.

(parameters:Postprocess/Gravity_20calculation)=
## **Subsection:** Postprocess / Gravity calculation
(parameters:Postprocess/Gravity_20calculation/List_20of_20latitude)=
### __Parameter name:__ List of latitude
**Default value:**

**Pattern:** [List of <[Double -90...90 (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** Parameter for the list of points sampling scheme: List of satellite latitude coordinates.

(parameters:Postprocess/Gravity_20calculation/List_20of_20longitude)=
### __Parameter name:__ List of longitude
**Default value:**

**Pattern:** [List of <[Double -180...180 (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** Parameter for the list of points sampling scheme: List of satellite longitude coordinates.

(parameters:Postprocess/Gravity_20calculation/List_20of_20radius)=
### __Parameter name:__ List of radius
**Default value:**

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** Parameter for the list of points sampling scheme: List of satellite radius coordinates. Just specify one radius if all points values have the same radius. If not, make sure there are as many radius as longitude and latitude

(parameters:Postprocess/Gravity_20calculation/Maximum_20latitude)=
### __Parameter name:__ Maximum latitude
**Default value:** 90

**Pattern:** [Double -90...90 (inclusive)]

**Documentation:** Parameter for the uniform distribution sampling scheme: Gravity may be calculated for a sets of points along the latitude between a minimum and maximum latitude.

(parameters:Postprocess/Gravity_20calculation/Maximum_20longitude)=
### __Parameter name:__ Maximum longitude
**Default value:** 180.

**Pattern:** [Double -180...180 (inclusive)]

**Documentation:** Parameter for the uniform distribution sampling scheme: Gravity may be calculated for a sets of points along the longitude between a minimum and maximum longitude.

(parameters:Postprocess/Gravity_20calculation/Maximum_20radius)=
### __Parameter name:__ Maximum radius
**Default value:** 0.

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Parameter for the map sampling scheme: Maximum radius can be defined in or outside the model.

(parameters:Postprocess/Gravity_20calculation/Minimum_20latitude)=
### __Parameter name:__ Minimum latitude
**Default value:** -90.

**Pattern:** [Double -90...90 (inclusive)]

**Documentation:** Parameter for the uniform distribution sampling scheme: Gravity may be calculated for a sets of points along the latitude between a minimum and maximum latitude.

(parameters:Postprocess/Gravity_20calculation/Minimum_20longitude)=
### __Parameter name:__ Minimum longitude
**Default value:** -180.

**Pattern:** [Double -180...180 (inclusive)]

**Documentation:** Parameter for the uniform distribution sampling scheme: Gravity may be calculated for a sets of points along the longitude between a minimum and maximum longitude.

(parameters:Postprocess/Gravity_20calculation/Minimum_20radius)=
### __Parameter name:__ Minimum radius
**Default value:** 0.

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Parameter for the map sampling scheme: Minimum radius may be defined in or outside the model. Prescribe a minimum radius for a sampling coverage at a specific height.

(parameters:Postprocess/Gravity_20calculation/Number_20points_20fibonacci_20spiral)=
### __Parameter name:__ Number points fibonacci spiral
**Default value:** 200

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** Parameter for the fibonacci spiral sampling scheme: This specifies the desired number of satellites per radius layer. The default value is 200. Note that sampling becomes more uniform with increasing number of satellites

(parameters:Postprocess/Gravity_20calculation/Number_20points_20latitude)=
### __Parameter name:__ Number points latitude
**Default value:** 1

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** Parameter for the map sampling scheme: This specifies the number of points along the latitude (e.g. gravity map) between a minimum and maximum latitude.

(parameters:Postprocess/Gravity_20calculation/Number_20points_20longitude)=
### __Parameter name:__ Number points longitude
**Default value:** 1

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** Parameter for the map sampling scheme: This specifies the number of points along the longitude (e.g. gravity map) between a minimum and maximum longitude.

(parameters:Postprocess/Gravity_20calculation/Number_20points_20radius)=
### __Parameter name:__ Number points radius
**Default value:** 1

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** Parameter for the map sampling scheme: This specifies the number of points along the radius (e.g. depth profile) between a minimum and maximum radius.

(parameters:Postprocess/Gravity_20calculation/Precision_20in_20gravity_20output)=
### __Parameter name:__ Precision in gravity output
**Default value:** 12

**Pattern:** [Integer range 1...2147483647 (inclusive)]

**Documentation:** Set the precision of gravity acceleration, potential and gradients in the gravity output and statistics file.

(parameters:Postprocess/Gravity_20calculation/Quadrature_20degree_20increase)=
### __Parameter name:__ Quadrature degree increase
**Default value:** 0

**Pattern:** [Integer range -1...2147483647 (inclusive)]

**Documentation:** Quadrature degree increase over the velocity element degree may be required when gravity is calculated near the surface or inside the model. An increase in the quadrature element adds accuracy to the gravity solution from noise due to the model grid.

(parameters:Postprocess/Gravity_20calculation/Reference_20density)=
### __Parameter name:__ Reference density
**Default value:** 3300.

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Gravity anomalies may be computed using density anomalies relative to a reference density.

(parameters:Postprocess/Gravity_20calculation/Sampling_20scheme)=
### __Parameter name:__ Sampling scheme
**Default value:** map

**Pattern:** [Selection map|list|list of points|fibonacci spiral ]

**Documentation:** Choose the sampling scheme. By default, the map produces a grid of equally angled points between a minimum and maximum radius, longitude, and latitude. A list of points contains the specific coordinates of the satellites. The fibonacci spiral sampling scheme produces a uniformly distributed map on the surface of sphere defined by a minimum and/or maximum radius.

(parameters:Postprocess/Gravity_20calculation/Time_20between_20gravity_20output)=
### __Parameter name:__ Time between gravity output
**Default value:** 1e8

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The time interval between each generation of gravity output files. A value of 0 indicates that output should be generated in each time step. Units: years if the &rsquo;Use years in output instead of seconds&rsquo; parameter is set; seconds otherwise.

(parameters:Postprocess/Gravity_20calculation/Time_20steps_20between_20gravity_20output)=
### __Parameter name:__ Time steps between gravity output
**Default value:** 2147483647

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** The maximum number of time steps between each generation of gravity output files.

(parameters:Postprocess/Memory_20statistics)=
## **Subsection:** Postprocess / Memory statistics
(parameters:Postprocess/Memory_20statistics/Output_20peak_20virtual_20memory_20_28VmPeak_29)=
### __Parameter name:__ Output peak virtual memory (VmPeak)
**Default value:** true

**Pattern:** [Bool]

**Documentation:** If set to &rsquo;true&rsquo;, also output the peak virtual memory usage (computed as the maximum over all processors).

(parameters:Postprocess/Particles)=
## **Subsection:** Postprocess / Particles
(parameters:Postprocess/Particles/Allow_20cells_20without_20particles)=
### __Parameter name:__ Allow cells without particles
**Default value:** false

**Pattern:** [Bool]

**Documentation:** By default, every cell needs to contain particles to use this interpolator plugin. If this parameter is set to true, cells are allowed to have no particles. In case both the current cell and its neighbors are empty, the interpolator will return 0 for the current cell&rsquo;s properties.

(parameters:Postprocess/Particles/Data_20output_20format)=
### __Parameter name:__ Data output format
**Default value:** vtu

**Pattern:** [MultipleSelection none|dx|ucd|gnuplot|povray|eps|gmv|tecplot|tecplot_binary|vtk|vtu|hdf5|svg|deal.II intermediate|ascii ]

**Documentation:** A comma separated list of file formats to be used for graphical output. The list of possible output formats that can be given here is documented in the appendix of the manual where the current parameter is described.

(parameters:Postprocess/Particles/Exclude_20output_20properties)=
### __Parameter name:__ Exclude output properties
**Default value:**

**Pattern:** [Anything]

**Documentation:** A comma separated list of particle properties that should *not* be output. If this list contains the entry &lsquo;all&rsquo;, only the id of particles will be provided in graphical output files.

(parameters:Postprocess/Particles/Integration_20scheme)=
### __Parameter name:__ Integration scheme
**Default value:** rk2

**Pattern:** [Selection euler|rk2|rk4 ]

**Documentation:** This parameter is used to decide which method to use to solve the equation that describes the position of particles, i.e., $\frac{d}{dt}\mathbf x_k(t) = \mathbf u(\mathbf x_k(t),t)$, where $k$ is an index that runs over all particles, and $\mathbf u(\mathbf x,t)$ is the velocity field that results from the Stokes equations.

In practice, the exact velocity $\mathbf u(\mathbf x,t)$ is of course not available, but only a numerical approximation $\mathbf u_h(\mathbf x,t)$. Furthermore, this approximation is only available at discrete time steps, $\mathbf u^n(\mathbf x)=\mathbf u(\mathbf x,t^n)$, and these need to be interpolated between time steps if the integrator for the equation above requires an evaluation at time points between the discrete time steps. If we denote this interpolation in time by $\tilde{\mathbf u}_h(\mathbf x,t)$ where $\tilde{\mathbf u}_h(\mathbf x,t^n)=\mathbf u^n(\mathbf x)$, then the equation the differential equation solver really tries to solve is $\frac{d}{dt}\tilde{\mathbf x}_k(t) =  \tilde{\mathbf u}_h(\mathbf x_k(t),t)$.

As a consequence of these considerations, if you try to assess convergence properties of an ODE integrator -- for example to verify that the RK4 integrator converges with fourth order --, it is important to recall that the integrator may not solve the equation you think it solves. If, for example, we call the numerical solution of the ODE $\tilde{\mathbf x}_{k,h}(t)$, then the error will typically satisfy a relationship like \[  \| \tilde{\mathbf x}_k(T) - \tilde{\mathbf x}_{k,h}(T) \|  \le  C(T) \Delta t^p\] where $\Delta t$ is the time step and $p$ the convergence order of the method, and $C(T)$ is a (generally unknown) constant that depends on the end time $T$ at which one compares the solutions. On the other hand, an analytically computed trajectory would likely use the *exact* velocity, and one may be tempted to compute $\| \mathbf x_k(T) - \tilde{\mathbf x}_{k,h}(T) \|$, but this quantity will, in the best case, only satisfy an estimate of the form \[  \| \mathbf x_k(T) - \tilde{\mathbf x}_{k,h}(T) \|  \le  C_1(T) \Delta t^p  + C_2(T) \| \mathbf u-\mathbf u_h \|  + C_3(T) \| \mathbf u_h-\tilde{\mathbf u}_h \|\] with appropriately chosen norms for the second and third term. These second and third terms typically converge to zero at relatively low rates (compared to the order $p$ of the integrator, which can often be chosen relatively high) in the mesh size $h$ and the time step size $\\Delta t$, limiting the overall accuracy of the ODE integrator.

Select one of the following models:

&lsquo;euler&rsquo;: Explicit Euler scheme integrator, where $y_{n+1} = y_n + \Delta t \, v(y_n)$. This requires only one integration substep per timestep.

&lsquo;rk2&rsquo;: Second Order Runge Kutta integrator $y_{n+1} = y_n + \Delta t\, v(t_{n+1/2}, y_{n} + \frac{1}{2} k_1)$ where $k_1 = \Delta t\, v(t_{n}, y_{n})$

&lsquo;rk4&rsquo;: Runge Kutta fourth order integrator, where $y_{n+1} = y_n + \frac{1}{6} k_1 + \frac{1}{3} k_2 + \frac{1}{3} k_3 + \frac{1}{6} k_4$ and $k_1$, $k_2$, $k_3$, $k_4$ are defined as usual.

(parameters:Postprocess/Particles/Interpolation_20scheme)=
### __Parameter name:__ Interpolation scheme
**Default value:** cell average

**Pattern:** [Selection bilinear least squares|cell average|harmonic average|nearest neighbor|quadratic least squares ]

**Documentation:** Select one of the following models:

&lsquo;bilinear least squares&rsquo;: Uses linear least squares to obtain the slopes and center of a 2D or 3D plane from the particle positions and a particular property value on those particles. Interpolate this property onto a vector of points. If the limiter is enabled then it will ensure the interpolated properties do not exceed the range of the minimum and maximum of the values of the property on the particles. Note that deal.II must be configured with BLAS and LAPACK to support this operation.

&lsquo;cell average&rsquo;: Return the arithmetic average of all particle properties in the given cell, or in the neighboring cells if the given cell is empty. In case the neighboring cells are also empty, and &rsquo;Allow cells without particles&rsquo; is set to true, the interpolator returns 0. Otherwise, an exception is thrown.

&lsquo;harmonic average&rsquo;: Return the harmonic average of all particle properties in the given cell. If the cell contains no particles, return the harmonic average of the properties in the neighboring cells. In case the neighboring cells are also empty, and &rsquo;Allow cells without particles&rsquo; is set to true, the interpolator returns 0. Otherwise, an exception is thrown.

&lsquo;nearest neighbor&rsquo;: Return the properties of the nearest neighboring particle in the current cell, or nearest particle in nearest neighboring cell if current cell is empty. In case the neighboring cells are also empty, and &rsquo;Allow cells without particles&rsquo; is set to true, the interpolator returns 0. Otherwise, an exception is thrown.

&lsquo;quadratic least squares&rsquo;: Interpolates particle properties onto a vector of points using a quadratic least squares method. Note that deal.II must be configured with BLAS/LAPACK.

(parameters:Postprocess/Particles/List_20of_20particle_20properties)=
### __Parameter name:__ List of particle properties
**Default value:**

**Pattern:** [MultipleSelection composition|elastic stress|function|initial composition|initial position|integrated strain|integrated strain invariant|melt particle|pT path|position|reference position|velocity|viscoplastic strain invariants ]

**Documentation:** A comma separated list of particle properties that should be tracked. By default none is selected, which means only position, velocity and id of the particles are output.

The following properties are available:

&lsquo;composition&rsquo;: Implementation of a plugin in which the particle property is defined by the compositional fields in the model. This can be used to track solid compositionevolution over time.

&lsquo;elastic stress&rsquo;: A plugin in which the particle property tensor is defined as the total elastic stress a particle has accumulated. See the viscoelastic material model documentation for more detailed information.

&lsquo;function&rsquo;: Implementation of a model in which the particle property is set by evaluating an explicit function at the initial position of each particle. The function is defined in the parameters in section &ldquo;Particles|Function&rdquo;. The format of these functions follows the syntax understood by the muparser library, see Section~\ref{sec:muparser-format}.

&lsquo;initial composition&rsquo;: Implementation of a plugin in which the particle property is given as the initial composition at the particle&rsquo;s initial position. The particle gets as many properties as there are compositional fields.

&lsquo;initial position&rsquo;: Implementation of a plugin in which the particle property is given as the initial position of the particle. This property is vector-valued with as many components as there are space dimensions. In practice, it is often most useful to only visualize one of the components of this vector, or the magnitude of the vector. For example, in a spherical mantle simulation, the magnitude of this property equals the starting radius of a particle, and is thereby indicative of which part of the mantle a particle comes from.

&lsquo;integrated strain&rsquo;: A plugin in which the particle property tensor is defined as the deformation gradient tensor $\mathbf F$ this particle has experienced. $\mathbf F$ can be polar-decomposed into the left stretching tensor $\mathbf L$ (the finite strain we are interested in), and the rotation tensor $\mathbf Q$. See the corresponding cookbook in the manual for more detailed information.

&lsquo;integrated strain invariant&rsquo;: A plugin in which the particle property is defined as the finite strain invariant ($\varepsilon_{ii}$). This property is calculated with the timestep ($dt$) and the second invariant of the deviatoric strain rate tensor ($\dot{\varepsilon}_{ii}$), where the value at time step $n$ is $\varepsilon_{ii}^{n} = \varepsilon_{ii}^{n-1} + dt\dot{\varepsilon}_{ii}$.

&lsquo;melt particle&rsquo;: Implementation of a plugin in which the particle property is defined as presence of melt above a threshold, which can be set as an input parameter. This property is set to 0 if melt is not present and set to 1 if melt is present.

&lsquo;pT path&rsquo;: Implementation of a plugin in which the particle property is defined as the current pressure and temperature at this position. This can be used to generate pressure-temperature paths of material points over time.

&lsquo;position&rsquo;: Implementation of a plugin in which the particle property is defined as the current position.

&lsquo;reference position&rsquo;: Implementation of a plugin in which the particle property is defined as the current reference position.

&lsquo;velocity&rsquo;: Implementation of a plugin in which the particle property is defined as the recent velocity at this position.

&lsquo;viscoplastic strain invariants&rsquo;: A plugin that calculates the finite strain invariant a particle has experienced and assigns it to either the plastic and/or viscous strain field based on whether the material is plastically yielding, or the total strain field used in the visco plastic material model. The implementation of this property is equivalent to the implementation for compositional fields that is located in the plugin in `benchmarks/buiter\_et\_al\_2008\_jgr/plugin/`,and is effectively the same as what the visco plastic material model uses for compositional fields.

(parameters:Postprocess/Particles/Load_20balancing_20strategy)=
### __Parameter name:__ Load balancing strategy
**Default value:** repartition

**Pattern:** [MultipleSelection none|remove particles|add particles|remove and add particles|repartition ]

**Documentation:** Strategy that is used to balance the computational load across processors for adaptive meshes.

(parameters:Postprocess/Particles/Maximum_20particles_20per_20cell)=
### __Parameter name:__ Maximum particles per cell
**Default value:** 100

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** Upper limit for particle number per cell. This limit is useful for adaptive meshes to prevent coarse cells from slowing down the whole model. It will be checked and enforced after mesh refinement, after MPI transfer of particles and after particle movement. If there are `n\_number\_of\_particles` $>$ `max\_particles\_per\_cell` particles in one cell then `n\_number\_of\_particles` - `max\_particles\_per\_cell` particles in this cell are randomly chosen and destroyed.

(parameters:Postprocess/Particles/Minimum_20particles_20per_20cell)=
### __Parameter name:__ Minimum particles per cell
**Default value:** 0

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** Lower limit for particle number per cell. This limit is useful for adaptive meshes to prevent fine cells from being empty of particles. It will be checked and enforced after mesh refinement and after particle movement. If there are `n\_number\_of\_particles` $<$ `min\_particles\_per\_cell` particles in one cell then `min\_particles\_per\_cell` - `n\_number\_of\_particles` particles are generated and randomly placed in this cell. If the particles carry properties the individual property plugins control how the properties of the new particles are initialized.

(parameters:Postprocess/Particles/Number_20of_20grouped_20files)=
### __Parameter name:__ Number of grouped files
**Default value:** 16

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** VTU file output supports grouping files from several CPUs into a given number of files using MPI I/O when writing on a parallel filesystem. Select 0 for no grouping. This will disable parallel file output and instead write one file per processor. A value of 1 will generate one big file containing the whole solution, while a larger value will create that many files (at most as many as there are MPI ranks).

(parameters:Postprocess/Particles/Number_20of_20particles)=
### __Parameter name:__ Number of particles
**Default value:** 1000

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Total number of particles to create (not per processor or per element). The number is parsed as a floating point number (so that one can specify, for example, &rsquo;1e4&rsquo; particles) but it is interpreted as an integer, of course.

(parameters:Postprocess/Particles/Particle_20generator_20name)=
### __Parameter name:__ Particle generator name
**Default value:** random uniform

**Pattern:** [Selection ascii file|probability density function|quadrature points|random uniform|reference cell|uniform box|uniform radial ]

**Documentation:** Select one of the following models:

&lsquo;ascii file&rsquo;: Generates a distribution of particles from coordinates specified in an Ascii data file. The file format is a simple text file, with as many columns as spatial dimensions and as many lines as particles to be generated. Initial comment lines starting with &lsquo;#&rsquo; will be discarded. Note that this plugin always generates as many particles as there are coordinates in the data file, the &ldquo;Postprocess/Particles/Number of particles&rdquo; parameter has no effect on this plugin. All of the values that define this generator are read from a section &ldquo;Postprocess/Particles/Generator/Ascii file&rdquo; in the input file, see Section~\ref{parameters:Postprocess/Particles/Generator/Ascii_20file}.

&lsquo;probability density function&rsquo;: Generate a random distribution of particles over the entire simulation domain. The probability density is prescribed in the form of a user-prescribed function. The format of this function follows the syntax understood by the muparser library, see Section~\ref{sec:muparser-format}. The return value of the function is always checked to be a non-negative probability density but it can be zero in parts of the domain.

&lsquo;quadrature points&rsquo;: Generates particles at the quadrature points of each active cell of the triangulation. Here, Gauss quadrature of degree (velocity\_degree + 1), is used similarly to the assembly of Stokes matrix.

&lsquo;random uniform&rsquo;: Generates a random uniform distribution of particles over the entire simulation domain.

&lsquo;reference cell&rsquo;: Generates a uniform distribution of particles per cell and spatial direction in the unit cell and transforms each of the particles back to real region in the model domain. Uniform here means the particles will be generated with an equal spacing in each spatial dimension.

&lsquo;uniform box&rsquo;: Generate a uniform distribution of particles over a rectangular domain in 2D or 3D. Uniform here means the particles will be generated with an equal spacing in each spatial dimension. Note that in order to produce a regular distribution the number of generated particles might not exactly match the one specified in the input file.

&lsquo;uniform radial&rsquo;: Generate a uniform distribution of particles over a spherical domain in 2D or 3D. Uniform here means the particles will be generated with an equal spacing in each spherical spatial dimension, i.e., the particles are created at positions that increase linearly with equal spacing in radius, colatitude and longitude around a certain center point. Note that in order to produce a regular distribution the number of generated particles might not exactly match the one specified in the input file.

(parameters:Postprocess/Particles/Particle_20weight)=
### __Parameter name:__ Particle weight
**Default value:** 10

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** Weight that is associated with the computational load of a single particle. The sum of particle weights will be added to the sum of cell weights to determine the partitioning of the mesh if the &lsquo;repartition&rsquo; particle load balancing strategy is selected. The optimal weight depends on the used integrator and particle properties. In general for a more expensive integrator and more expensive properties a larger particle weight is recommended. Before adding the weights of particles, each cell already carries a weight of 1000 to account for the cost of field-based computations.

(parameters:Postprocess/Particles/Temporary_20output_20location)=
### __Parameter name:__ Temporary output location
**Default value:**

**Pattern:** [Anything]

**Documentation:** On large clusters it can be advantageous to first write the output to a temporary file on a local file system and later move this file to a network file system. If this variable is set to a non-empty string it will be interpreted as a temporary storage location.

(parameters:Postprocess/Particles/Time_20between_20data_20output)=
### __Parameter name:__ Time between data output
**Default value:** 1e8

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The time interval between each generation of output files. A value of zero indicates that output should be generated every time step.

Units: years if the &rsquo;Use years in output instead of seconds&rsquo; parameter is set; seconds otherwise.

(parameters:Postprocess/Particles/Update_20ghost_20particles)=
### __Parameter name:__ Update ghost particles
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Some particle interpolation algorithms require knowledge about particles in neighboring cells. To allow this, particles in ghost cells need to be exchanged between the processes neighboring this cell. This parameter determines whether this transport is happening.

(parameters:Postprocess/Particles/Write_20in_20background_20thread)=
### __Parameter name:__ Write in background thread
**Default value:** false

**Pattern:** [Bool]

**Documentation:** File operations can potentially take a long time, blocking the progress of the rest of the model run. Setting this variable to &lsquo;true&rsquo; moves this process into a background thread, while the rest of the model continues.

(parameters:Postprocess/Particles/Function)=
## **Subsection:** Postprocess / Particles / Function
(parameters:Postprocess/Particles/Function/Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**

**Pattern:** [Anything]

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form &lsquo;var1=value1, var2=value2, ...&rsquo;.

A typical example would be to set this runtime parameter to &lsquo;pi=3.1415926536&rsquo; and then use &lsquo;pi&rsquo; in the expression of the actual formula. (That said, for convenience this class actually defines both &lsquo;pi&rsquo; and &lsquo;Pi&rsquo; by default, but you get the idea.)

(parameters:Postprocess/Particles/Function/Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0

**Pattern:** [Anything]

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as &lsquo;sin&rsquo; or &lsquo;cos&rsquo;. In addition, it may contain expressions like &lsquo;if(x>0, 1, -1)&rsquo; where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon.

(parameters:Postprocess/Particles/Function/Number_20of_20components)=
### __Parameter name:__ Number of components
**Default value:** 1

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** The number of function components where each component is described by a function expression delimited by a &rsquo;;&rsquo;.

(parameters:Postprocess/Particles/Function/Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,y,t

**Pattern:** [Anything]

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are &lsquo;x&rsquo; (in 1d), &lsquo;x,y&rsquo; (in 2d) or &lsquo;x,y,z&rsquo; (in 3d) for spatial coordinates and &lsquo;t&rsquo; for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to &lsquo;r,phi,theta,t&rsquo; and then use these variable names in your function expression.

(parameters:Postprocess/Particles/Generator)=
## **Subsection:** Postprocess / Particles / Generator
(parameters:Postprocess/Particles/Generator/Ascii_20file)=
## **Subsection:** Postprocess / Particles / Generator / Ascii file
(parameters:Postprocess/Particles/Generator/Ascii_20file/Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/particle/generator/ascii/

**Pattern:** [DirectoryName]

**Documentation:** The name of a directory that contains the particle data. This path may either be absolute (if starting with a &rsquo;/&rsquo;) or relative to the current directory. The path may also include the special text &rsquo;$ASPECT_SOURCE_DIR&rsquo; which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the &lsquo;data/&rsquo; subdirectory of ASPECT.

(parameters:Postprocess/Particles/Generator/Ascii_20file/Data_20file_20name)=
### __Parameter name:__ Data file name
**Default value:** particle.dat

**Pattern:** [Anything]

**Documentation:** The name of the particle file.

(parameters:Postprocess/Particles/Generator/Probability_20density_20function)=
## **Subsection:** Postprocess / Particles / Generator / Probability density function
(parameters:Postprocess/Particles/Generator/Probability_20density_20function/Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**

**Pattern:** [Anything]

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form &lsquo;var1=value1, var2=value2, ...&rsquo;.

A typical example would be to set this runtime parameter to &lsquo;pi=3.1415926536&rsquo; and then use &lsquo;pi&rsquo; in the expression of the actual formula. (That said, for convenience this class actually defines both &lsquo;pi&rsquo; and &lsquo;Pi&rsquo; by default, but you get the idea.)

(parameters:Postprocess/Particles/Generator/Probability_20density_20function/Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0

**Pattern:** [Anything]

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as &lsquo;sin&rsquo; or &lsquo;cos&rsquo;. In addition, it may contain expressions like &lsquo;if(x>0, 1, -1)&rsquo; where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon.

(parameters:Postprocess/Particles/Generator/Probability_20density_20function/Random_20cell_20selection)=
### __Parameter name:__ Random cell selection
**Default value:** true

**Pattern:** [Bool]

**Documentation:** If true, particle numbers per cell are calculated randomly according to their respective probability density. This means particle numbers per cell can deviate statistically from the integral of the probability density. If false, first determine how many particles each cell should have based on the integral of the density over each of the cells, and then once we know how many particles we want on each cell, choose their locations randomly within each cell.

(parameters:Postprocess/Particles/Generator/Probability_20density_20function/Random_20number_20seed)=
### __Parameter name:__ Random number seed
**Default value:** 5432

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** The seed for the random number generator that controls the particle generation. Keep constant to generate identical particle distributions in subsequent model runs. Change to get a different distribution. In parallel computations the seed is further modified on each process to ensure different particle patterns on different processes. Note that the number of particles per processor is not affected by the seed.

(parameters:Postprocess/Particles/Generator/Probability_20density_20function/Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,y,t

**Pattern:** [Anything]

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are &lsquo;x&rsquo; (in 1d), &lsquo;x,y&rsquo; (in 2d) or &lsquo;x,y,z&rsquo; (in 3d) for spatial coordinates and &lsquo;t&rsquo; for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to &lsquo;r,phi,theta,t&rsquo; and then use these variable names in your function expression.

(parameters:Postprocess/Particles/Generator/Reference_20cell)=
## **Subsection:** Postprocess / Particles / Generator / Reference cell
(parameters:Postprocess/Particles/Generator/Reference_20cell/Number_20of_20particles_20per_20cell_20per_20direction)=
### __Parameter name:__ Number of particles per cell per direction
**Default value:** 2

**Pattern:** [List of <[Integer range 1...2147483647 (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** List of number of particles to create per cell and spatial dimension. The size of the list is the number of spatial dimensions. If only one value is given, then each spatial dimension is set to the same value. The list of numbers are parsed as a floating point number (so that one can specify, for example, &rsquo;1e4&rsquo; particles) but it is interpreted as an integer, of course.

(parameters:Postprocess/Particles/Generator/Uniform_20box)=
## **Subsection:** Postprocess / Particles / Generator / Uniform box
(parameters:Postprocess/Particles/Generator/Uniform_20box/Maximum_20x)=
### __Parameter name:__ Maximum x
**Default value:** 1.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Maximum x coordinate for the region of particles.

(parameters:Postprocess/Particles/Generator/Uniform_20box/Maximum_20y)=
### __Parameter name:__ Maximum y
**Default value:** 1.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Maximum y coordinate for the region of particles.

(parameters:Postprocess/Particles/Generator/Uniform_20box/Maximum_20z)=
### __Parameter name:__ Maximum z
**Default value:** 1.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Maximum z coordinate for the region of particles.

(parameters:Postprocess/Particles/Generator/Uniform_20box/Minimum_20x)=
### __Parameter name:__ Minimum x
**Default value:** 0.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Minimum x coordinate for the region of particles.

(parameters:Postprocess/Particles/Generator/Uniform_20box/Minimum_20y)=
### __Parameter name:__ Minimum y
**Default value:** 0.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Minimum y coordinate for the region of particles.

(parameters:Postprocess/Particles/Generator/Uniform_20box/Minimum_20z)=
### __Parameter name:__ Minimum z
**Default value:** 0.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Minimum z coordinate for the region of particles.

(parameters:Postprocess/Particles/Generator/Uniform_20radial)=
## **Subsection:** Postprocess / Particles / Generator / Uniform radial
(parameters:Postprocess/Particles/Generator/Uniform_20radial/Center_20x)=
### __Parameter name:__ Center x
**Default value:** 0.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** x coordinate for the center of the spherical region, where particles are generated.

(parameters:Postprocess/Particles/Generator/Uniform_20radial/Center_20y)=
### __Parameter name:__ Center y
**Default value:** 0.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** y coordinate for the center of the spherical region, where particles are generated.

(parameters:Postprocess/Particles/Generator/Uniform_20radial/Center_20z)=
### __Parameter name:__ Center z
**Default value:** 0.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** z coordinate for the center of the spherical region, where particles are generated.

(parameters:Postprocess/Particles/Generator/Uniform_20radial/Maximum_20latitude)=
### __Parameter name:__ Maximum latitude
**Default value:** 180.

**Pattern:** [Double 0...180 (inclusive)]

**Documentation:** Maximum latitude coordinate for the region of particles in degrees. Measured from the center position, and from the north pole.

(parameters:Postprocess/Particles/Generator/Uniform_20radial/Maximum_20longitude)=
### __Parameter name:__ Maximum longitude
**Default value:** 360.

**Pattern:** [Double -180...360 (inclusive)]

**Documentation:** Maximum longitude coordinate for the region of particles in degrees. Measured from the center position.

(parameters:Postprocess/Particles/Generator/Uniform_20radial/Maximum_20radius)=
### __Parameter name:__ Maximum radius
**Default value:** 1.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Maximum radial coordinate for the region of particles. Measured from the center position.

(parameters:Postprocess/Particles/Generator/Uniform_20radial/Minimum_20latitude)=
### __Parameter name:__ Minimum latitude
**Default value:** 0.

**Pattern:** [Double 0...180 (inclusive)]

**Documentation:** Minimum latitude coordinate for the region of particles in degrees. Measured from the center position, and from the north pole.

(parameters:Postprocess/Particles/Generator/Uniform_20radial/Minimum_20longitude)=
### __Parameter name:__ Minimum longitude
**Default value:** 0.

**Pattern:** [Double -180...360 (inclusive)]

**Documentation:** Minimum longitude coordinate for the region of particles in degrees. Measured from the center position.

(parameters:Postprocess/Particles/Generator/Uniform_20radial/Minimum_20radius)=
### __Parameter name:__ Minimum radius
**Default value:** 0.

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Minimum radial coordinate for the region of particles. Measured from the center position.

(parameters:Postprocess/Particles/Generator/Uniform_20radial/Radial_20layers)=
### __Parameter name:__ Radial layers
**Default value:** 1

**Pattern:** [Integer range 1...2147483647 (inclusive)]

**Documentation:** The number of radial shells of particles that will be generated around the central point.

(parameters:Postprocess/Particles/Interpolator)=
## **Subsection:** Postprocess / Particles / Interpolator
(parameters:Postprocess/Particles/Interpolator/Bilinear_20least_20squares)=
## **Subsection:** Postprocess / Particles / Interpolator / Bilinear least squares
(parameters:Postprocess/Particles/Interpolator/Bilinear_20least_20squares/Use_20boundary_20extrapolation)=
### __Parameter name:__ Use boundary extrapolation
**Default value:** false

**Pattern:** [List of <[Bool]> of length 0...4294967295 (inclusive)]

**Documentation:** Extends the range used by &rsquo;Use linear least squares limiter&rsquo; by linearly interpolating values at cell boundaries from neighboring cells. If more than one value is given, it will be treated as a list with one component per particle property. Enabling &rsquo;Use boundary extrapolation&rsquo; requires enabling &rsquo;Use linear least squares limiter&rsquo;.

(parameters:Postprocess/Particles/Interpolator/Bilinear_20least_20squares/Use_20linear_20least_20squares_20limiter)=
### __Parameter name:__ Use linear least squares limiter
**Default value:** false

**Pattern:** [List of <[Bool]> of length 0...4294967295 (inclusive)]

**Documentation:** Limit the interpolation of particle properties onto the cell, so that the value of each property is no smaller than its minimum and no larger than its maximum on the particles of each cell, and the average of neighboring cells. If more than one value is given, it will be treated as a list with one component per particle property.

(parameters:Postprocess/Particles/Interpolator/Quadratic_20least_20squares)=
## **Subsection:** Postprocess / Particles / Interpolator / Quadratic least squares
(parameters:Postprocess/Particles/Interpolator/Quadratic_20least_20squares/Use_20boundary_20extrapolation)=
### __Parameter name:__ Use boundary extrapolation
**Default value:** false

**Pattern:** [List of <[Bool]> of length 0...4294967295 (inclusive)]

**Documentation:** Extends the range used by &rsquo;Use quadratic least squares limiter&rsquo; by linearly interpolating values at cell boundaries from neighboring cells. If more than one value is given, it will be treated as a list with one component per particle property. Enabling &rsquo;Use boundary extrapolation&rsquo; requires enabling &rsquo;Use quadratic least squares limiter&rsquo;.

(parameters:Postprocess/Particles/Interpolator/Quadratic_20least_20squares/Use_20quadratic_20least_20squares_20limiter)=
### __Parameter name:__ Use quadratic least squares limiter
**Default value:** true

**Pattern:** [List of <[Bool]> of length 0...4294967295 (inclusive)]

**Documentation:** Limit the interpolation of particle properties onto the cell, so that the value of each property is no smaller than its minimum and no larger than its maximum on the particles of each cell, and the average of neighboring cells. If more than one value is given, it will be treated as a list with one component per particle property.

(parameters:Postprocess/Particles/Melt_20particle)=
## **Subsection:** Postprocess / Particles / Melt particle
(parameters:Postprocess/Particles/Melt_20particle/Threshold_20for_20melt_20presence)=
### __Parameter name:__ Threshold for melt presence
**Default value:** 1e-3

**Pattern:** [Double 0...1 (inclusive)]

**Documentation:** The minimum porosity that has to be present at the position of a particle for it to be considered a melt particle (in the sense that the melt presence property is set to 1).

(parameters:Postprocess/Point_20values)=
## **Subsection:** Postprocess / Point values
(parameters:Postprocess/Point_20values/Evaluation_20points)=
### __Parameter name:__ Evaluation points
**Default value:**

**Pattern:** [List of <[List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 2...2 (inclusive)]> of length 0...4294967295 (inclusive) separated by <;>]

**Documentation:** The list of points at which the solution should be evaluated. Points need to be separated by semicolons, and coordinates of each point need to be separated by commas.

(parameters:Postprocess/Point_20values/Time_20between_20point_20values_20output)=
### __Parameter name:__ Time between point values output
**Default value:** 0.

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The time interval between each generation of point values output. A value of zero indicates that output should be generated in each time step. Units: years if the &rsquo;Use years in output instead of seconds&rsquo; parameter is set; seconds otherwise.

(parameters:Postprocess/Point_20values/Use_20natural_20coordinates)=
### __Parameter name:__ Use natural coordinates
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Whether or not the Evaluation points are specified in the natural coordinates of the geometry model, e.g. radius, lon, lat for the chunk model. Currently, natural coordinates for the spherical shell and sphere geometries are not supported.

(parameters:Postprocess/Rotation_20statistics)=
## **Subsection:** Postprocess / Rotation statistics
(parameters:Postprocess/Rotation_20statistics/Output_20full_20moment_20of_20inertia_20tensor)=
### __Parameter name:__ Output full moment of inertia tensor
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Whether to write the full moment of inertia tensor into the statistics output instead of its norm for the current rotation axis. This is a second-order symmetric tensor with 6 components in 3D. In 2D this option has no effect, because the rotation axis is fixed and thus the moment of inertia is always a scalar.

(parameters:Postprocess/Rotation_20statistics/Use_20constant_20density_20of_20one)=
### __Parameter name:__ Use constant density of one
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Whether to use a constant density of one for the computation of the angular momentum and moment of inertia. This is an approximation that assumes that the &rsquo;volumetric&rsquo; rotation is equal to the &rsquo;mass&rsquo; rotation. If this parameter is true this postprocessor computes &rsquo;net rotation&rsquo; instead of &rsquo;angular momentum&rsquo;.

(parameters:Postprocess/Topography)=
## **Subsection:** Postprocess / Topography
(parameters:Postprocess/Topography/Output_20to_20file)=
### __Parameter name:__ Output to file
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Whether or not to write topography to a text file named named &rsquo;topography.NNNNN&rsquo; in the output directory

(parameters:Postprocess/Topography/Time_20between_20text_20output)=
### __Parameter name:__ Time between text output
**Default value:** 0.

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The time interval between each generation of text output files. A value of zero indicates that output should be generated in each time step. Units: years if the &rsquo;Use years in output instead of seconds&rsquo; parameter is set; seconds otherwise.

(parameters:Postprocess/Visualization)=
## **Subsection:** Postprocess / Visualization
(parameters:Postprocess/Visualization/Filter_20output)=
### __Parameter name:__ Filter output
**Default value:** false

**Pattern:** [Bool]

**Documentation:** deal.II offers the possibility to filter duplicate vertices for HDF5 output files. This merges the vertices of adjacent cells and therefore saves disk space, but misrepresents discontinuous output properties. Activating this function reduces the disk space by about a factor of $2^{dim}$ for HDF5 output, and currently has no effect on other output formats. \note{**Warning:** Setting this flag to true will result in visualization output that does not accurately represent discontinuous fields. This may be because you are using a discontinuous finite element for the pressure, temperature, or compositional variables, or because you use a visualization postprocessor that outputs quantities as discontinuous fields (e.g., the strain rate, viscosity, etc.). These will then all be visualized as *continuous* quantities even though, internally, ASPECT considers them as discontinuous fields.}

(parameters:Postprocess/Visualization/Interpolate_20output)=
### __Parameter name:__ Interpolate output
**Default value:** true

**Pattern:** [Bool]

**Documentation:** deal.II offers the possibility to linearly interpolate output fields of higher order elements to a finer resolution. This somewhat compensates the fact that most visualization software only offers linear interpolation between grid points and therefore the output file is a very coarse representation of the actual solution field. Activating this option increases the spatial resolution in each dimension by a factor equal to the polynomial degree used for the velocity finite element (usually 2). In other words, instead of showing one quadrilateral or hexahedron in the visualization per cell on which ASPECT computes, it shows multiple (for quadratic elements, it will describe each cell of the mesh on which we compute as $2\times 2$ or $2\times 2\times 2$ cells in 2d and 3d, respectively; correspondingly more subdivisions are used if you use cubic, quartic, or even higher order elements for the velocity).

The effect of using this option can be seen in the following picture showing a variation of the output produced with the input files from Section~\ref{sec:shell-simple-2d}:

\begin{center}  \includegraphics[width=0.5\textwidth]{viz/parameters/build-patches}\end{center}Here, the left picture shows one visualization cell per computational cell (i.e., the option is switched off), and the right picture shows the same simulation with the option switched on (which is the default). The images show the same data, demonstrating that interpolating the solution onto bilinear shape functions as is commonly done in visualizing data loses information.

Of course, activating this option also greatly increases the amount of data ASPECT will write to disk: approximately by a factor of 4 in 2d, and a factor of 8 in 3d, when using quadratic elements for the velocity, and correspondingly more for even higher order elements.

(parameters:Postprocess/Visualization/List_20of_20output_20variables)=
### __Parameter name:__ List of output variables
**Default value:**

**Pattern:** [MultipleSelection ISA rotation timescale|Vp anomaly|Vs anomaly|adiabat|artificial viscosity|artificial viscosity composition|boundary indicators|boundary strain rate residual|boundary velocity residual|compositional vector|depth|dynamic topography|error indicator|geoid|grain lag angle|gravity|heat flux map|heating|material properties|maximum horizontal compressive stress|melt fraction|melt material properties|named additional outputs|nonadiabatic pressure|nonadiabatic temperature|particle count|partition|principal stress|shear stress|spd factor|spherical velocity components|strain rate|strain rate tensor|stress|stress second invariant|surface dynamic topography|surface stress|temperature anomaly|vertical heat flux|volume of fluid values|volumetric strain rate|density|specific heat|thermal conductivity|thermal diffusivity|thermal expansivity|viscosity ]

**Documentation:** A comma separated list of visualization objects that should be run whenever writing graphical output. By default, the graphical output files will always contain the primary variables velocity, pressure, and temperature. However, one frequently wants to also visualize derived quantities, such as the thermodynamic phase that corresponds to a given temperature-pressure value, or the corresponding seismic wave speeds. The visualization objects do exactly this: they compute such derived quantities and place them into the output file. The current parameter is the place where you decide which of these additional output variables you want to have in your output file.

The following postprocessors are available:

&lsquo;ISA rotation timescale&rsquo;: A visualization output object that generates output showing the timescale for the rotation of grains toward the infinite strain axis. Kaminski and Ribe (see \cite{Kaminski2002}) call this quantity $\tau_\text{ISA}$ and define it as $\tau_\text{ISA} \approx \frac{1}{\dot{\epsilon}}$ where $\dot{\epsilon}$ is the largest eigenvalue of the strain rate tensor. It can be used, along with the grain lag angle $\Theta$, to calculate the grain orientation lag parameter.

Physical units: \si{\second}.

&lsquo;Vp anomaly&rsquo;: A visualization output object that generates output showing the percentage anomaly in the seismic compressional wave speed $V_p$ as a spatially variable function with one value per cell. This anomaly is either shown as a percentage anomaly relative to the reference profile given by adiabatic conditions (with the compositions given by the current composition, such that the reference could potentially change through time), or as a percentage change relative to the laterally averaged velocity at the depth of the cell. This velocity is calculated by linear interpolation between average values calculated within equally thick depth slices. The number of depth slices in the domain is user-defined. Typically, the best results will be obtained if the number of depth slices is balanced between being large enough to capture step changes in velocities, but small enough to maintain a reasonable number of evaluation points per slice. Bear in mind that lateral averaging subsamples the finite element mesh. Note that this plugin requires a material model that provides seismic velocities.

Physical units: None, the quantity being output is a fractional change provided as a percentage.

&lsquo;Vs anomaly&rsquo;: A visualization output object that generates output showing the percentage anomaly in the seismic shear wave speed $V_s$ as a spatially variable function with one value per cell. This anomaly is either shown as a percentage anomaly relative to the reference profile given by adiabatic conditions (with the compositions given by the current composition, such that the reference could potentially change through time), or as a percentage change relative to the laterally averaged velocity at the depth of the cell. This velocity is calculated by linear interpolation between average values calculated within equally thick depth slices. The number of depth slices in the domain is user-defined. Typically, the best results will be obtained if the number of depth slices is balanced between being large enough to capture step changes in velocities, but small enough to maintain a reasonable number of evaluation points per slice. Bear in mind that lateral averaging subsamples the finite element mesh. Note that this plugin requires a material model that provides seismic velocities.

Physical units: None, the quantity being output is a fractional change provided as a percentage.

&lsquo;adiabat&rsquo;: A visualization output object that generates adiabatic temperature, pressure, density, and density derivative (with regard to depth)as produced by the `AdiabaticConditions` class.

Physical units: \si{\kelvin}, \si{\pascal}, \si{\kilo\gram\per\meter\cubed\per\meter}, respectively, for the four components.

&lsquo;artificial viscosity&rsquo;: A visualization output object that generates output showing the value of the artificial viscosity on each cell.

Physical units: \si{\watt\per\meter\per\kelvin}.

&lsquo;artificial viscosity composition&rsquo;: A visualization output object that generates output showing the value of the artificial viscosity for a compositional field on each cell.

Physical units: \si{\meter\squared\per\second}.

&lsquo;boundary indicators&rsquo;: A visualization output object that generates output about the used boundary indicators. In a loop over the active cells, if a cell lies at a domain boundary, the boundary indicator of the face along the boundary is requested. In case the cell does not lie along any domain boundary, the cell is assigned the value of the largest used boundary indicator plus one. When a cell is situated in one of the corners of the domain, multiple faces will have a boundary indicator. This postprocessor returns the value of the first face along a boundary that is encountered in a loop over all the faces.

Physical units: None.

&lsquo;boundary strain rate residual&rsquo;: A visualization output object that generates output for the strain rate residual at the top surface. The residual is computed at each point at the surface as the difference between the strain rate invariant in the model and the input data, where the invariant is computed like in the &rsquo;strain rate&rsquo; postprocessor. The user chooses the input data as ascii data files with coordinate columns and column corresponding to the surface strain rate norm.

Physical units: $\frac{1}{\text{s}}$ or $\frac{1}{\text{year}}$, depending on settings in the input file.

&lsquo;boundary velocity residual&rsquo;: A visualization output object that generates output for the velocity residual at the top surface. The residual is computed at each point at the surface as the difference between the modeled velocities and the input data velocities for each vector component. The user has an option to choose the input data as ascii data files (e.g. GPS velocities) with columns in the same format as described for the &rsquo;ascii data&rsquo; initial temperature plugin or a velocity field computed from the GPlates program as described in the gplates boundary velocity plugin.

Physical units: $\frac{\text{m}}{\text{s}}$ or $\frac{\text{m}}{\text{year}}$, depending on settings in the input file.

&lsquo;compositional vector&rsquo;: A visualization output object that outputs vectors whose components are derived from compositional fields. Input parameters for this postprocessor are defined in section Postprocess/Visualization/Compositional fields as vectors.

Physical units: None.

&lsquo;depth&rsquo;: A visualization output postprocessor that outputs the depth for all points inside the domain, as determined by the geometry model.

Physical units: \si{\meter}.

&lsquo;dynamic topography&rsquo;: A visualization output object that generates output for the dynamic topography at the top and bottom of the model space. The approach to determine the dynamic topography requires us to compute the stress tensor and evaluate the component of it in the direction in which gravity acts. In other words, we compute $\sigma_{rr}={\hat g}^T(2 \eta \varepsilon(\mathbf u)-\frac 13 (\textrm{div}\;\mathbf u)I)\hat g - p_d$ where $\hat g = \mathbf g/\|\mathbf g\|$ is the direction of the gravity vector $\mathbf g$ and $p_d=p-p_a$ is the dynamic pressure computed by subtracting the adiabatic pressure $p_a$ from the total pressure $p$ computed as part of the Stokes solve. From this, the dynamic topography is computed using the formula $h=\frac{\sigma_{rr}}{(\mathbf g \cdot \mathbf n)  \rho}$ where $\rho$ is the density at the cell center. For the bottom surface we chose the convection that positive values are up (out) and negative values are in (down), analogous to the deformation of the upper surface. Note that this implementation takes the direction of gravity into account, which means that reversing the flow in backward advection calculations will not reverse the instantaneous topography because the reverse flow will be divided by the reverse surface gravity.

Strictly speaking, the dynamic topography is of course a quantity that is only of interest at the surface. However, we compute it everywhere to make things fit into the framework within which we produce data for visualization. You probably only want to visualize whatever data this postprocessor generates at the surface of your domain and simply ignore the rest of the data generated.

Alternatively, consider using the "surface dynamic topography" visualization postprocessor to only output the dynamic topography at the boundary of the domain.

Physical units: \si{\meter}.

&lsquo;error indicator&rsquo;: A visualization output object that generates output showing the estimated error or other mesh refinement indicator as a spatially variable function with one value per cell.

Physical units: None. (Strictly speaking, errors have physical units of course, but because error *indicators* can be computed from different solution components and other input, we consider error indicators unitless.)

&lsquo;geoid&rsquo;: Visualization for the geoid solution. The geoid is given by the equivalent water column height due to a gravity perturbation.

Physical units: \si{\meter}.

&lsquo;grain lag angle&rsquo;: A visualization output object that generates output showing the angle between the ~infinite strain axis and the flow velocity. Kaminski and Ribe (see \cite{Kaminski2002}) call this quantity $\Theta$ and define it as $\Theta = \cos^{-1}(\hat{u}\cdot\hat{e})$  where $\hat{u}=\vec{u}/|{u}|$, $\vec{u}$ is the local flow velocity, and $\hat{e}$ is the local infinite strain axis, which we calculate as the first eigenvector of the &rsquo;left stretch&rsquo; tensor. $\Theta$ can be used to calculate the grain orientation lag parameter.

Physical units: \si{\radian}.

&lsquo;gravity&rsquo;: A visualization output object that outputs the gravity vector.

Physical units: \si {\meter\per\second\squared} .

&lsquo;heat flux map&rsquo;: A visualization output object that generates output for the heat flux density across the top and bottom boundary in outward direction. The heat flux is computed as sum of advective heat flux and conductive heat flux through Neumann boundaries, both computed as integral over the boundary area, and conductive heat flux through Dirichlet boundaries, which is computed using the consistent boundary flux method as described in &ldquo;Gresho, Lee, Sani, Maslanik, Eaton (1987). The consistent Galerkin FEM for computing derived boundary quantities in thermal and or fluids problems. International Journal for Numerical Methods in Fluids, 7(4), 371-394.&rdquo; If only conductive heat flux through Dirichlet boundaries is of interest, the postprocessor can produce output of higher resolution by evaluating the CBF solution vector point-wise instead of computing cell-wise averaged values.

Physical units: \si{\watt\per\meter\squared}.

&lsquo;heating&rsquo;: A visualization output object that generates output for all the heating terms used in the energy equation.

Physical units: \si{\watt\per\cubic\meter}.

&lsquo;material properties&rsquo;: A visualization output object that generates output for the material properties given by the material model. The current postprocessor allows to output a (potentially large) subset of all of the information provided by material models at once, with just a single material model evaluation per output point. Although individual properties can still be listed in the &ldquo;List of output variables&rdquo;, this visualization plugin is called internally to avoid duplicated evaluations of the material model.

In almost all places inside ASPECT, the program can use &ldquo;averaged&rdquo; material properties, for example for the assembly of matrices and right hand side vectors. To accurately reflect the material parameters used internally, this visualization postprocessor averages in the same way as is used to do the assembly, and consequently the graphical output will reflect not pointwise properties, but averaged properties.

Physical units: Various.

&lsquo;maximum horizontal compressive stress&rsquo;: A plugin that computes the direction and magnitude of the maximum horizontal component of the compressive stress as a vector field. The direction of this vector can often be used to visualize the principal mode of deformation (e.g., at normal faults or extensional margins) and can be correlated with seismic anisotropy. Recall that the *compressive* stress is simply the negative stress, $\sigma_c=-\sigma=-\left[     2\eta (\varepsilon(\mathbf u)             - \frac 13 (\nabla \cdot \mathbf u) I)     + pI\right]$.

Following \cite{LundTownend07}, we define the maximum horizontal stress direction as that *horizontal* direction $\mathbf n$ that maximizes $\mathbf n^T \sigma_c \mathbf n$. We call a vector *horizontal* if it is perpendicular to the gravity vector $\mathbf g$.

In two space dimensions, $\mathbf n$ is simply a vector that is horizontal (we choose one of the two possible choices). This direction is then scaled by the size of the horizontal stress in this direction, i.e., the plugin outputs the vector $\mathbf w = (\mathbf n^T \sigma_c \mathbf n) \; \mathbf n$.

In three space dimensions, given two horizontal, perpendicular, unit length, but otherwise arbitrarily chosen vectors $\mathbf u,\mathbf v$, we can express $\mathbf n = (\cos \alpha)\mathbf u + (\sin\alpha)\mathbf v$ where $\alpha$ maximizes the expression \begin{align*}  f(\alpha) = \mathbf n^T \sigma_c \mathbf n  = (\mathbf u^T \sigma_c \mathbf u)(\cos\alpha)^2    +2(\mathbf u^T \sigma_c \mathbf v)(\cos\alpha)(\sin\alpha)    +(\mathbf v^T \sigma_c \mathbf v)(\sin\alpha)^2.\end{align*}

The maximum of $f(\alpha)$ is attained where $f&rsquo;(\alpha)=0$. Evaluating the derivative and using trigonometric identities, one finds that $\alpha$ has to satisfy the equation \begin{align*}  \tan(2\alpha) = \frac{2.0\mathbf u^T \sigma_c \mathbf v}                          {\mathbf u^T \sigma_c \mathbf u                            - \mathbf v^T \sigma_c \mathbf v}.\end{align*}Since the transform $\alpha\mapsto\alpha+\pi$ flips the direction of $\mathbf n$, we only need to seek a solution to this equation in the interval $\alpha\in[0,\pi)$. These are given by $\alpha_1=\frac 12 \arctan \frac{\mathbf u^T \sigma_c \mathbf v}{\mathbf u^T \sigma_c \mathbf u - \mathbf v^T \sigma_c \mathbf v}$ and $\alpha_2=\alpha_1+\frac{\pi}{2}$, one of which will correspond to a minimum and the other to a maximum of $f(\alpha)$. One checks the sign of $f&rdquo;(\alpha)=-2(\mathbf u^T \sigma_c \mathbf u - \mathbf v^T \sigma_c \mathbf v)\cos(2\alpha) - 2 (\mathbf u^T \sigma_c \mathbf v) \sin(2\alpha)$ for each of these to determine the $\alpha$ that maximizes $f(\alpha)$, and from this immediately arrives at the correct form for the maximum horizontal stress $\mathbf n$.

The description above computes a 3d *direction* vector $\mathbf n$. If one were to scale this vector the same way as done in 2d, i.e., with the magnitude of the stress in this direction, one will typically get vectors whose length is principally determined by the hydrostatic pressure at a given location simply because the hydrostatic pressure is the largest component of the overall stress. On the other hand, the hydrostatic pressure does not determine any principal direction because it is an isotropic, anti-compressive force. As a consequence, there are often points in simulations (e.g., at the center of convection rolls) where the stress has no dominant horizontal direction, and the algorithm above will then in essence choose a random direction because the stress is approximately equal in all horizontal directions. If one scaled the output by the magnitude of the stress in this direction (i.e., approximately equal to the hydrostatic pressure at this point), one would get randomly oriented vectors at these locations with significant lengths.

To avoid this problem, we scale the maximal horizontal compressive stress direction $\mathbf n$ by the *difference* between the stress in the maximal and minimal horizontal stress directions. In other words, let $\mathbf n_\perp=(\sin \alpha)\mathbf u - (\cos\alpha)\mathbf v$ be the horizontal direction perpendicular to $\mathbf n$, then this plugin outputs the vector quantity $\mathbf w = (\mathbf n^T \sigma_c \mathbf n                -\mathbf n^T_\perp \sigma_c \mathbf n_\perp)               \; \mathbf n$. In other words, the length of the vector produced indicates *how dominant* the direction of maximal horizontal compressive strength is.

Fig.~\ref{fig:max-horizontal-compressive-stress} shows a simple example for this kind of visualization in 3d.

\begin{figure}  \includegraphics[width=0.3\textwidth]    {viz/plugins/maximum_horizontal_compressive_stress/temperature.png}  \hfill  \includegraphics[width=0.3\textwidth]    {viz/plugins/maximum_horizontal_compressive_stress/velocity.png}  \hfill  \includegraphics[width=0.3\textwidth]    {viz/plugins/maximum_horizontal_compressive_stress/horizontal-stress.png}  \caption{\it Illustration of the &lsquo;maximum horizontal     compressive stress&rsquo; visualization plugin. The left     figure shows a ridge-like temperature anomaly. Together     with no-slip boundary along all six boundaries, this     results in two convection rolls (center). The maximal     horizontal compressive strength at the bottom center     of the domain is perpendicular to the ridge because     the flow comes together there from the left and right,     yielding a compressive force in left-right direction.     At the top of the model, the flow separates outward,     leading to a *negative* compressive stress     in left-right direction; because there is no flow     in front-back direction, the compressive strength     in front-back direction is zero, making the along-ridge     direction the dominant one. At the center of the     convection rolls, both horizontal directions yield     the same stress; the plugin therefore chooses an     essentially arbitrary horizontal vector, but then     uses a zero magnitude given that the difference     between the maximal and minimal horizontal stress     is zero at these points.}  \label{fig:max-horizontal-compressive-stress}\end{figure}

Physical units: \si{\pascal}.

&lsquo;melt fraction&rsquo;: A visualization output object that generates output for the melt fraction at the temperature and pressure of the current point. If the material model computes a melt fraction, this is the quantity that will be visualized. Otherwise, a specific parametrization for batch melting (as described in the following) will be used. It does not take into account latent heat. If there are no compositional fields, or no fields called &rsquo;pyroxenite&rsquo;,  this postprocessor will visualize the melt fraction of peridotite (calculated using the anhydrous model of Katz, 2003). If there is a compositional field called &rsquo;pyroxenite&rsquo;, the postprocessor assumes that this compositional field is the content of pyroxenite, and will visualize the melt fraction for a mixture of peridotite and pyroxenite (using the melting model of Sobolev, 2011 for pyroxenite). All the parameters that were used in these calculations can be changed in the input file, the most relevant maybe being the mass fraction of Cpx in peridotite in the Katz melting model (Mass fraction cpx), which right now has a default of 15\%. The corresponding $p$-$T$-diagrams can be generated by running the tests melt\_postprocessor\_peridotite and melt\_postprocessor\_pyroxenite.

Physical units: None.

&lsquo;melt material properties&rsquo;: A visualization output object that generates output for melt related properties of the material model. Note that this postprocessor always outputs the compaction pressure, but can output a large range of additional properties, as selected in the &ldquo;List of properties&rdquo; parameter.

Physical units: Various, depending on what is being output.

&lsquo;named additional outputs&rsquo;: Some material models can compute quantities other than those that typically appear in the equations that ASPECT solves (such as the viscosity, density, etc). Examples of quantities material models may be able to compute are seismic velocities, or other quantities that can be derived from the state variables and the material coefficients such as the stress or stress anisotropies. These quantities are generically referred to as &lsquo;named outputs&rsquo; because they are given an explicit name different from the usual outputs of material models.

This visualization postprocessor outputs whatever quantities the material model can compute. What quantities these are is specific to the material model in use for a simulation, and for many models in fact does not contain any named outputs at all.

Physical units: Various, depending on what is being output.

&lsquo;nonadiabatic pressure&rsquo;: A visualization output object that generates output for the non-adiabatic component of the pressure.

The variable that is outputted this way is computed by taking the pressure at each point and subtracting from it the adiabatic pressure computed at the beginning of the simulation. Because the adiabatic pressure is one way of defining a static pressure background field, what this visualization postprocessor therefore produces is *one* way to compute a *dynamic pressure*. There are, however, other ways as well, depending on the choice of the &ldquo;background pressure&rdquo;.

Physical units: \si{\pascal}.

&lsquo;nonadiabatic temperature&rsquo;: A visualization output object that generates output for the non-adiabatic component of the temperature.

Physical units: \si{\kelvin}.

&lsquo;particle count&rsquo;: A visualization output object that generates output about the number of particles per cell.

Physical units: None.

&lsquo;partition&rsquo;: A visualization output object that generates output for the parallel partition that every cell of the mesh is associated with.

Physical units: None.

&lsquo;principal stress&rsquo;: A visualization output object that outputs the principal stress values and directions, i.e., the eigenvalues and eigenvectors of the stress tensor. The postprocessor can either operate on the full stress tensor or only on the deviatoric stress tensor, depending on what run-time parameters are set.

Physical units: \si{\pascal}.

&lsquo;shear stress&rsquo;: A visualization output object that generates output for the 3 (in 2d) or 6 (in 3d) components of the shear stress tensor, i.e., for the components of the tensor $-2\eta\varepsilon(\mathbf u)$ in the incompressible case and $-2\eta\left[\varepsilon(\mathbf u)-\tfrac 13(\textrm{tr}\;\varepsilon(\mathbf u))\mathbf I\right]$ in the compressible case. If elasticity is used, the elastic contribution is being accounted for. The shear stress differs from the full stress tensor by the absence of the pressure. Note that the convention of positive compressive stress is followed.

Physical units: \si{\pascal}.

&lsquo;spd factor&rsquo;: A visualization output object that generates output for the spd factor. The spd factor is a factor which scales a part of the Jacobian used for the Newton solver to make sure that the Jacobian remains positive definite.

Physical units: None.

&lsquo;spherical velocity components&rsquo;: A visualization output object that outputs the polar coordinates components $v_r$ and $v_\phi$ of the velocity field in 2D and the spherical coordinates components $v_r$, $v_{\phi}$ and $v_{\theta}$ of the velocity field in 3D.

Physical units: $\frac{\text{m}}{\text{s}}$ or $\frac{\text{m}}{\text{year}}$, depending on settings in the input file.

&lsquo;strain rate&rsquo;: A visualization output object that generates output for the norm of the strain rate, i.e., for the quantity $\sqrt{\varepsilon(\mathbf u):\varepsilon(\mathbf u)}$ in the incompressible case and $\sqrt{[\varepsilon(\mathbf u)-\tfrac 13(\textrm{tr}\;\varepsilon(\mathbf u))\mathbf I]:[\varepsilon(\mathbf u)-\tfrac 13(\textrm{tr}\;\varepsilon(\mathbf u))\mathbf I]}$ in the compressible case.

Physical units: \si{\per\second}.

&lsquo;strain rate tensor&rsquo;: A visualization output object that generates output for the 4 (in 2d) or 9 (in 3d) components of the strain rate tensor, i.e., for the components of the tensor $\varepsilon(\mathbf u)$ in the incompressible case and $\varepsilon(\mathbf u)-\tfrac 13(\textrm{tr}\;\varepsilon(\mathbf u))\mathbf I$ in the compressible case.

Physical units: \si{\per\second}.

&lsquo;stress&rsquo;: A visualization output object that generates output for the 3 (in 2d) or 6 (in 3d) components of the stress tensor, i.e., for the components of the tensor $-2\eta\varepsilon(\mathbf u)+pI$ in the incompressible case and $-2\eta\left[\varepsilon(\mathbf u)-\tfrac 13(\textrm{tr}\;\varepsilon(\mathbf u))\mathbf I\right]+pI$ in the compressible case. If elasticity is used, the elastic contribution is being accounted for. Note that the convention of positive compressive stress is followed.

Physical units: \si{\pascal}.

&lsquo;stress second invariant&rsquo;: A visualization output object that outputs the second moment invariant of the deviatoric stress tensor.

Physical units: \si{\pascal}.

&lsquo;surface dynamic topography&rsquo;: A visualization output object that generates output for the dynamic topography at the top and bottom of the model space. The approach to determine the dynamic topography requires us to compute the stress tensor and evaluate the component of it in the direction in which gravity acts. In other words, we compute $\sigma_{rr}={\hat g}^T(2 \eta \varepsilon(\mathbf u)-\frac 13 (\textrm{div}\;\mathbf u)I)\hat g - p_d$ where $\hat g = \mathbf g/\|\mathbf g\|$ is the direction of the gravity vector $\mathbf g$ and $p_d=p-p_a$ is the dynamic pressure computed by subtracting the adiabatic pressure $p_a$ from the total pressure $p$ computed as part of the Stokes solve. From this, the dynamic topography is computed using the formula $h=\frac{\sigma_{rr}}{(\mathbf g \cdot \mathbf n)  \rho}$ where $\rho$ is the density at the cell center. For the bottom surface we chose the convection that positive values are up (out) and negative values are in (down), analogous to the deformation of the upper surface. Note that this implementation takes the direction of gravity into account, which means that reversing the flow in backward advection calculations will not reverse the instantaneous topography because the reverse flow will be divided by the reverse surface gravity.

In contrast to the &lsquo;dynamic topography&rsquo; visualization postprocessor, this plugin really only evaluates the dynamic topography at faces of cells that are adjacent to &lsquo;bottom&rsquo; and &lsquo;top&rsquo; boundaries, and only outputs information on the surface of the domain, rather than padding the information with zeros in the interior of the domain.

Physical units: \si{\meter}.

&lsquo;surface stress&rsquo;: A visualization output object that generates output on the surface of the domain for the 3 (in 2d) or 6 (in 3d) components of the stress tensor, i.e., for the components of the tensor $-2\eta\varepsilon(\mathbf u)+pI$ in the incompressible case and $-2\eta\left[\varepsilon(\mathbf u)-\tfrac 13(\textrm{tr}\;\varepsilon(\mathbf u))\mathbf I\right]+pI$ in the compressible case. If elasticity is included, its contribution is accounted for. Note that the convention of positive compressive stress is followed.

Physical units: \si{\pascal}.

&lsquo;temperature anomaly&rsquo;: A visualization output postprocessor that outputs the temperature minus the depth-average of the temperature.The average temperature is calculated using the lateral averaging function from the &ldquo;depth average&rdquo; postprocessor and interpolated linearly between the layers specified through &ldquo;Number of depth slices&rdquo;.

Physical units: \si{\kelvin}.

&lsquo;vertical heat flux&rsquo;: A visualization output object that generates output for the heat flux in the vertical direction, which is the sum of the advective and the conductive heat flux, with the sign convention of positive flux upwards.

Physical units: \si{\watt\per\square\meter}.

&lsquo;volume of fluid values&rsquo;: A visualization output object that outputs the volume fraction and optionally a level set field and the interface normal vectors of volume of fluid fields.

Physical units: None.

&lsquo;volumetric strain rate&rsquo;: A visualization output object that generates output for the volumetric strain rate, i.e., for the quantity $\nabla\cdot\mathbf u = \textrm{div}\; \mathbf u = \textrm{trace}\; \varepsilon(\mathbf u)$. This should be zero (in some average sense) in incompressible convection models, but can be non-zero in compressible models and models with melt transport.

Physical units: \si{\per\second}.

(parameters:Postprocess/Visualization/Number_20of_20grouped_20files)=
### __Parameter name:__ Number of grouped files
**Default value:** 16

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** VTU file output supports grouping files from several CPUs into a given number of files using MPI I/O when writing on a parallel filesystem. Select 0 for no grouping. This will disable parallel file output and instead write one file per processor. A value of 1 will generate one big file containing the whole solution, while a larger value will create that many files (at most as many as there are MPI ranks).

(parameters:Postprocess/Visualization/Output_20format)=
### __Parameter name:__ Output format
**Default value:** vtu

**Pattern:** [Selection none|dx|ucd|gnuplot|povray|eps|gmv|tecplot|tecplot_binary|vtk|vtu|hdf5|svg|deal.II intermediate ]

**Documentation:** The file format to be used for graphical output. The list of possible output formats that can be given here is documented in the appendix of the manual where the current parameter is described.

(parameters:Postprocess/Visualization/Output_20mesh_20velocity)=
### __Parameter name:__ Output mesh velocity
**Default value:** false

**Pattern:** [Bool]

**Documentation:** For computations with deforming meshes, ASPECT uses an Arbitrary-Lagrangian-Eulerian formulation to handle deforming the domain, so the mesh has its own velocity field.  This may be written as an output field by setting this parameter to true.

(parameters:Postprocess/Visualization/Point_2dwise_20stress_20and_20strain)=
### __Parameter name:__ Point-wise stress and strain
**Default value:** false

**Pattern:** [Bool]

**Documentation:** If set to true, quantities related to stress and strain are computed in each vertex. Otherwise, an average per cell is computed.

(parameters:Postprocess/Visualization/Temporary_20output_20location)=
### __Parameter name:__ Temporary output location
**Default value:**

**Pattern:** [Anything]

**Documentation:** On large clusters it can be advantageous to first write the output to a temporary file on a local file system and later move this file to a network file system. If this variable is set to a non-empty string it will be interpreted as a temporary storage location.

(parameters:Postprocess/Visualization/Time_20between_20graphical_20output)=
### __Parameter name:__ Time between graphical output
**Default value:** 1e8

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The time interval between each generation of graphical output files. A value of zero indicates that output should be generated in each time step. Units: years if the &rsquo;Use years in output instead of seconds&rsquo; parameter is set; seconds otherwise.

(parameters:Postprocess/Visualization/Time_20steps_20between_20graphical_20output)=
### __Parameter name:__ Time steps between graphical output
**Default value:** 2147483647

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** The maximum number of time steps between each generation of graphical output files.

(parameters:Postprocess/Visualization/Write_20higher_20order_20output)=
### __Parameter name:__ Write higher order output
**Default value:** false

**Pattern:** [Bool]

**Documentation:** deal.II offers the possibility to write vtu files with higher order representations of the output data. This means each cell will correctly show the higher order representation of the output data instead of the linear interpolation between vertices that ParaView and VisIt usually show. Note that activating this option is safe and recommended, but requires that (i) &ldquo;Output format&rdquo; is set to &ldquo;vtu&rdquo;, (ii) &ldquo;Interpolate output&rdquo; is set to true, (iii) you use a sufficiently new version of Paraview or VisIt to read the files (Paraview version 5.5 or newer, and VisIt version to be determined), and (iv) you use deal.II version 9.1.0 or newer.
The effect of using this option can be seen in the following picture:

\begin{center}  \includegraphics[width=0.5\textwidth]{viz/parameters/higher-order-output}\end{center}The top figure shows the plain output without interpolation or higher order output. The middle figure shows output that was interpolated as discussed for the &ldquo;Interpolate output&rdquo; option. The bottom panel shows higher order output that achieves better accuracy than the interpolated output at a lower memory cost.

(parameters:Postprocess/Visualization/Write_20in_20background_20thread)=
### __Parameter name:__ Write in background thread
**Default value:** false

**Pattern:** [Bool]

**Documentation:** File operations can potentially take a long time, blocking the progress of the rest of the model run. Setting this variable to &lsquo;true&rsquo; moves this process into a background thread, while the rest of the model continues.

(parameters:Postprocess/Visualization/Artificial_20viscosity_20composition)=
## **Subsection:** Postprocess / Visualization / Artificial viscosity composition
(parameters:Postprocess/Visualization/Artificial_20viscosity_20composition/Name_20of_20compositional_20field)=
### __Parameter name:__ Name of compositional field
**Default value:**

**Pattern:** [Anything]

**Documentation:** The name of the compositional field whose output should be visualized.

(parameters:Postprocess/Visualization/Compositional_20fields_20as_20vectors)=
## **Subsection:** Postprocess / Visualization / Compositional fields as vectors
(parameters:Postprocess/Visualization/Compositional_20fields_20as_20vectors/Names_20of_20fields)=
### __Parameter name:__ Names of fields
**Default value:**

**Pattern:** [Anything]

**Documentation:** A list of sets of compositional fields which should be output as vectors. Sets are separated from each other by semicolons and vector components within each set are separated by commas (e.g. $vec1_x$, $vec1_y$ ; $vec2_x$, $vec2_y$) where each name must be a defined named compositional field. If only one name is given in a set, it is interpreted as the first in a sequence of dim consecutive compositional fields.

(parameters:Postprocess/Visualization/Compositional_20fields_20as_20vectors/Names_20of_20vectors)=
### __Parameter name:__ Names of vectors
**Default value:**

**Pattern:** [List of <[Anything]> of length 0...4294967295 (inclusive)]

**Documentation:** Names of vectors as they will appear in the output.

(parameters:Postprocess/Visualization/Heat_20flux_20map)=
## **Subsection:** Postprocess / Visualization / Heat flux map
(parameters:Postprocess/Visualization/Heat_20flux_20map/Output_20point_20wise_20heat_20flux)=
### __Parameter name:__ Output point wise heat flux
**Default value:** false

**Pattern:** [Bool]

**Documentation:** A boolean flag that controls whether to output the heat flux map as a point wise value, or as a cell-wise averaged value. The point wise output is more accurate, but it currently omits prescribed heat flux values at boundaries and advective heat flux that is caused by velocities non-tangential to boundaries. If you do not use these two features it is recommended to switch this setting on to benefit from the increased output resolution.

(parameters:Postprocess/Visualization/Material_20properties)=
## **Subsection:** Postprocess / Visualization / Material properties
(parameters:Postprocess/Visualization/Material_20properties/List_20of_20material_20properties)=
### __Parameter name:__ List of material properties
**Default value:** density,thermal expansivity,specific heat,viscosity

**Pattern:** [MultipleSelection viscosity|density|thermal expansivity|specific heat|thermal conductivity|thermal diffusivity|compressibility|entropy derivative temperature|entropy derivative pressure|reaction terms|melt fraction ]

**Documentation:** A comma separated list of material properties that should be written whenever writing graphical output. By default, the material properties will always contain the density, thermal expansivity, specific heat and viscosity. The following material properties are available:

viscosity|density|thermal expansivity|specific heat|thermal conductivity|thermal diffusivity|compressibility|entropy derivative temperature|entropy derivative pressure|reaction terms|melt fraction

(parameters:Postprocess/Visualization/Melt_20fraction)=
## **Subsection:** Postprocess / Visualization / Melt fraction
(parameters:Postprocess/Visualization/Melt_20fraction/A1)=
### __Parameter name:__ A1
**Default value:** 1085.7

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Constant parameter in the quadratic function that approximates the solidus of peridotite. Units: \si{\degreeCelsius}.

(parameters:Postprocess/Visualization/Melt_20fraction/A2)=
### __Parameter name:__ A2
**Default value:** 1.329e-7

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Prefactor of the linear pressure term in the quadratic function that approximates the solidus of peridotite. \si{\degreeCelsius\per\pascal}.

(parameters:Postprocess/Visualization/Melt_20fraction/A3)=
### __Parameter name:__ A3
**Default value:** -5.1e-18

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Prefactor of the quadratic pressure term in the quadratic function that approximates the solidus of peridotite. \si{\degreeCelsius\per\pascal\squared}.

(parameters:Postprocess/Visualization/Melt_20fraction/B1)=
### __Parameter name:__ B1
**Default value:** 1475.0

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Constant parameter in the quadratic function that approximates the lherzolite liquidus used for calculating the fraction of peridotite-derived melt. Units: \si{\degreeCelsius}.

(parameters:Postprocess/Visualization/Melt_20fraction/B2)=
### __Parameter name:__ B2
**Default value:** 8.0e-8

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Prefactor of the linear pressure term in the quadratic function that approximates the  lherzolite liquidus used for calculating the fraction of peridotite-derived melt. \si{\degreeCelsius\per\pascal}.

(parameters:Postprocess/Visualization/Melt_20fraction/B3)=
### __Parameter name:__ B3
**Default value:** -3.2e-18

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Prefactor of the quadratic pressure term in the quadratic function that approximates the  lherzolite liquidus used for calculating the fraction of peridotite-derived melt. \si{\degreeCelsius\per\pascal\squared}.

(parameters:Postprocess/Visualization/Melt_20fraction/C1)=
### __Parameter name:__ C1
**Default value:** 1780.0

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Constant parameter in the quadratic function that approximates the liquidus of peridotite. Units: \si{\degreeCelsius}.

(parameters:Postprocess/Visualization/Melt_20fraction/C2)=
### __Parameter name:__ C2
**Default value:** 4.50e-8

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Prefactor of the linear pressure term in the quadratic function that approximates the liquidus of peridotite. \si{\degreeCelsius\per\pascal}.

(parameters:Postprocess/Visualization/Melt_20fraction/C3)=
### __Parameter name:__ C3
**Default value:** -2.0e-18

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Prefactor of the quadratic pressure term in the quadratic function that approximates the liquidus of peridotite. \si{\degreeCelsius\per\pascal\squared}.

(parameters:Postprocess/Visualization/Melt_20fraction/D1)=
### __Parameter name:__ D1
**Default value:** 976.0

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Constant parameter in the quadratic function that approximates the solidus of pyroxenite. Units: \si{\degreeCelsius}.

(parameters:Postprocess/Visualization/Melt_20fraction/D2)=
### __Parameter name:__ D2
**Default value:** 1.329e-7

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Prefactor of the linear pressure term in the quadratic function that approximates the solidus of pyroxenite. Note that this factor is different from the value given in Sobolev, 2011, because they use the potential temperature whereas we use the absolute temperature. \si{\degreeCelsius\per\pascal}.

(parameters:Postprocess/Visualization/Melt_20fraction/D3)=
### __Parameter name:__ D3
**Default value:** -5.1e-18

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Prefactor of the quadratic pressure term in the quadratic function that approximates the solidus of pyroxenite. \si{\degreeCelsius\per\pascal\squared}.

(parameters:Postprocess/Visualization/Melt_20fraction/E1)=
### __Parameter name:__ E1
**Default value:** 663.8

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Prefactor of the linear depletion term in the quadratic function that approximates the melt fraction of pyroxenite. \si{\degreeCelsius\per\pascal}.

(parameters:Postprocess/Visualization/Melt_20fraction/E2)=
### __Parameter name:__ E2
**Default value:** -611.4

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Prefactor of the quadratic depletion term in the quadratic function that approximates the melt fraction of pyroxenite. \si{\degreeCelsius\per\pascal\squared}.

(parameters:Postprocess/Visualization/Melt_20fraction/Mass_20fraction_20cpx)=
### __Parameter name:__ Mass fraction cpx
**Default value:** 0.15

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Mass fraction of clinopyroxene in the peridotite to be molten. Units: non-dimensional.

(parameters:Postprocess/Visualization/Melt_20fraction/beta)=
### __Parameter name:__ beta
**Default value:** 1.5

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Exponent of the melting temperature in the melt fraction calculation. Units: non-dimensional.

(parameters:Postprocess/Visualization/Melt_20fraction/r1)=
### __Parameter name:__ r1
**Default value:** 0.5

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Constant in the linear function that approximates the clinopyroxene reaction coefficient. Units: non-dimensional.

(parameters:Postprocess/Visualization/Melt_20fraction/r2)=
### __Parameter name:__ r2
**Default value:** 8e-11

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Prefactor of the linear pressure term in the linear function that approximates the clinopyroxene reaction coefficient. Units: \si{\per\pascal}.

(parameters:Postprocess/Visualization/Melt_20material_20properties)=
## **Subsection:** Postprocess / Visualization / Melt material properties
(parameters:Postprocess/Visualization/Melt_20material_20properties/List_20of_20properties)=
### __Parameter name:__ List of properties
**Default value:** compaction viscosity,permeability

**Pattern:** [MultipleSelection compaction viscosity|fluid viscosity|permeability|fluid density|fluid density gradient|is melt cell|darcy coefficient|darcy coefficient no cutoff|compaction length ]

**Documentation:** A comma separated list of melt properties that should be written whenever writing graphical output. The following material properties are available:

compaction viscosity|fluid viscosity|permeability|fluid density|fluid density gradient|is melt cell|darcy coefficient|darcy coefficient no cutoff|compaction length

(parameters:Postprocess/Visualization/Principal_20stress)=
## **Subsection:** Postprocess / Visualization / Principal stress
(parameters:Postprocess/Visualization/Principal_20stress/Use_20deviatoric_20stress)=
### __Parameter name:__ Use deviatoric stress
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Whether to use the deviatoric stress tensor instead of the full stress tensor to compute principal stress directions and values.

(parameters:Postprocess/Visualization/Temperature_20anomaly)=
## **Subsection:** Postprocess / Visualization / Temperature anomaly
(parameters:Postprocess/Visualization/Temperature_20anomaly/Number_20of_20depth_20slices)=
### __Parameter name:__ Number of depth slices
**Default value:** 20

**Pattern:** [Integer range 1...2147483647 (inclusive)]

**Documentation:** Number of depth slices used to define average temperature.

(parameters:Postprocess/Visualization/Temperature_20anomaly/Use_20maximal_20temperature_20for_20bottom)=
### __Parameter name:__ Use maximal temperature for bottom
**Default value:** true

**Pattern:** [Bool]

**Documentation:** If true, use the specified boundary temperatures as average temperatures at the surface. If false, extrapolate the temperature gradient between the first and second cells to the surface. This option will only work for models with a fixed surface temperature.

(parameters:Postprocess/Visualization/Temperature_20anomaly/Use_20minimal_20temperature_20for_20surface)=
### __Parameter name:__ Use minimal temperature for surface
**Default value:** true

**Pattern:** [Bool]

**Documentation:** Whether to use the minimal specified boundary temperature as the bottom boundary temperature. This option will only work for models with a fixed bottom boundary temperature.

(parameters:Postprocess/Visualization/Volume_20of_20Fluid)=
## **Subsection:** Postprocess / Visualization / Volume of Fluid
(parameters:Postprocess/Visualization/Volume_20of_20Fluid/Output_20interface_20normals)=
### __Parameter name:__ Output interface normals
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Include the internal data for the interface normal on the unit cells.

(parameters:Postprocess/Visualization/Volume_20of_20Fluid/Output_20interface_20reconstruction_20contour)=
### __Parameter name:__ Output interface reconstruction contour
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Include fields defined such that the 0 contour is the fluid interface.

(parameters:Postprocess/Visualization/Vp_20anomaly)=
## **Subsection:** Postprocess / Visualization / Vp anomaly
(parameters:Postprocess/Visualization/Vp_20anomaly/Average_20velocity_20scheme)=
### __Parameter name:__ Average velocity scheme
**Default value:** reference profile

**Pattern:** [Selection reference profile|lateral average ]

**Documentation:** Scheme to compute the average velocity-depth profile. The reference profile option evaluates the conditions along the reference adiabat according to the material model. The lateral average option instead calculates a lateral average from subdivision of the mesh. The lateral average option may produce spurious results where there are sharp velocity changes.

(parameters:Postprocess/Visualization/Vp_20anomaly/Number_20of_20depth_20slices)=
### __Parameter name:__ Number of depth slices
**Default value:** 50

**Pattern:** [Integer range 1...2147483647 (inclusive)]

**Documentation:** Number of depth slices used to define average seismic compressional wave velocities from which anomalies are calculated. Units: non-dimensional.

(parameters:Postprocess/Visualization/Vs_20anomaly)=
## **Subsection:** Postprocess / Visualization / Vs anomaly
(parameters:Postprocess/Visualization/Vs_20anomaly/Average_20velocity_20scheme)=
### __Parameter name:__ Average velocity scheme
**Default value:** reference profile

**Pattern:** [Selection reference profile|lateral average ]

**Documentation:** Scheme to compute the average velocity-depth profile. The reference profile option evaluates the conditions along the reference adiabat according to the material model. The lateral average option instead calculates a lateral average from subdivision of the mesh. The lateral average option may produce spurious results where there are sharp velocity changes.

(parameters:Postprocess/Visualization/Vs_20anomaly/Number_20of_20depth_20slices)=
### __Parameter name:__ Number of depth slices
**Default value:** 50

**Pattern:** [Integer range 1...2147483647 (inclusive)]

**Documentation:** Number of depth slices used to define average seismic shear wave velocities from which anomalies are calculated. Units: non-dimensional.
