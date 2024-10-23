(parameters:global)=
# Global parameters


## **Subsection:** No subsection


(parameters:Additional_20shared_20libraries)=
### __Parameter name:__ Additional shared libraries
**Default value:**

**Pattern:** [List of <[FileName (Type: input)]> of length 0...4294967295 (inclusive)]

**Documentation:** A list of names of additional shared libraries that should be loaded upon starting up the program. The names of these files can contain absolute or relative paths (relative to the directory in which you call ASPECT). In fact, file names that do not contain any directory information (i.e., only the name of a file such as <libmyplugin.so> will not be found if they are not located in one of the directories listed in the `LD_LIBRARY_PATH` environment variable. In order to load a library in the current directory, use <./libmyplugin.so> instead.

If you specify <./libmyplugin.so>, ASPECT will open either <./libmyplugin.debug.so> or <./libmyplugin.release.so> depending on the current ASPECT build type.

The typical use of this parameter is so that you can implement additional plugins in your own directories, rather than in the ASPECT source directories. You can then simply compile these plugins into a shared library without having to re-compile all of ASPECT. See the section of the manual discussing writing extensions for more information on how to compile additional files into a shared library.

(parameters:Adiabatic_20surface_20temperature)=
### __Parameter name:__ Adiabatic surface temperature
**Default value:** 0.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** In order to make the problem in the first time step easier to solve, we need a reasonable guess for the temperature and pressure. To obtain it, we use an adiabatic pressure and temperature field. This parameter describes what the &lsquo;adiabatic&rsquo; temperature would be at the surface of the domain (i.e. at depth zero). Note that this value need not coincide with the boundary condition posed at this point. Rather, the boundary condition may differ significantly from the adiabatic value, and then typically induce a thermal boundary layer.

For more information, see the section in the manual that discusses the general mathematical model.

(parameters:CFL_20number)=
### __Parameter name:__ CFL number
**Default value:** 1.0

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** In computations, the time step $k$ is chosen according to $k = c \min_K \frac {h_K} {\|u\|_{\infty,K} p_T}$ where $h_K$ is the diameter of cell $K$, and the denominator is the maximal magnitude of the velocity on cell $K$ times the polynomial degree $p_T$ of the temperature discretization. The dimensionless constant $c$ is called the CFL number in this program. For time discretizations that have explicit components, $c$ must be less than a constant that depends on the details of the time discretization and that is no larger than one. On the other hand, for implicit discretizations such as the one chosen here, one can choose the time step as large as one wants (in particular, one can choose $c>1$) though a CFL number significantly larger than one will yield rather diffusive solutions. Units: None.

(parameters:Dimension)=
### __Parameter name:__ Dimension
**Default value:** 2

**Pattern:** [Integer range 2...3 (inclusive)]

**Documentation:** The number of space dimensions you want to run this program in. ASPECT can run in 2 and 3 space dimensions.

(parameters:End_20time)=
### __Parameter name:__ End time
**Default value:** 5.69e+300

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** The end time of the simulation. The default value is a number so that when converted from years to seconds it is approximately equal to the largest number representable in floating point arithmetic. For all practical purposes, this equals infinity. Units: Years if the &rsquo;Use years in output instead of seconds&rsquo; parameter is set; seconds otherwise.

(parameters:Max_20nonlinear_20iterations)=
### __Parameter name:__ Max nonlinear iterations
**Default value:** 10

**Pattern:** [Integer range 1...2147483647 (inclusive)]

**Documentation:** The maximal number of nonlinear iterations to be performed.

(parameters:Max_20nonlinear_20iterations_20in_20pre_2drefinement)=
### __Parameter name:__ Max nonlinear iterations in pre-refinement
**Default value:** 2147483647

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** The maximal number of nonlinear iterations to be performed in the pre-refinement steps. This does not include the last refinement step before moving to timestep 1. When this parameter has a larger value than max nonlinear iterations, the latter is used.

(parameters:Maximum_20first_20time_20step)=
### __Parameter name:__ Maximum first time step
**Default value:** 5.69e+300

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Set a maximum time step size for only the first timestep. Generally the time step based on the CFL number should be sufficient, but for complicated models or benchmarking it may be useful to limit the first time step to some value, especially when using the free surface, which needs to settle to prevent instabilities. This should in that case be combined with a value set for &ldquo;Maximum relative increase in time step&rdquo;. The default value is a value so that when converted from years into seconds it equals the largest number representable by a floating point number, implying an unlimited time step. Units: Years or seconds, depending on the &ldquo;Use years in output instead of seconds&rdquo; parameter.

(parameters:Maximum_20relative_20increase_20in_20time_20step)=
### __Parameter name:__ Maximum relative increase in time step
**Default value:** 91.0

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Set a percentage with which the length of the time step is limited to increase. Generally the time step based on the CFL number should be sufficient, but for complicated models which may suddenly drastically change behavior, it may be useful to limit the increase in the time step, without limiting the time step size of the whole simulation to a particular number. For example, if this parameter is set to $50$, then that means that the length of a time step can at most increase by 50\% from one time step to the next, or by a factor of 1.5.

Here, the default value is set to be 91\% because the best available step-size ratio bound guaranteeing stability in the PDE context seems to be 1.91, see {cite}`Denner:2014`. In that thesis, the bound was proved in the context of semilinear parabolic problem, but it appears reasonable to also use this value as an upper bound in the current context.

Units: \%.

(parameters:Maximum_20time_20step)=
### __Parameter name:__ Maximum time step
**Default value:** 5.69e+300

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Set a maximum time step size for the solver to use. Generally the time step based on the CFL number should be sufficient, but for complicated models or benchmarking it may be useful to limit the time step to some value. The default value is a value so that when converted from years into seconds it equals the largest number representable by a floating point number, implying an unlimited time step.Units: Years or seconds, depending on the &ldquo;Use years in output instead of seconds&rdquo; parameter.

(parameters:Nonlinear_20solver_20failure_20strategy)=
### __Parameter name:__ Nonlinear solver failure strategy
**Default value:** continue with next timestep

**Pattern:** [Selection continue with next timestep|cut timestep size|abort program ]

**Documentation:** Select the strategy on what to do if the nonlinear solver scheme fails to converge. The options are:
&lsquo;continue with next timestep&lsquo;: ignore error and continue to the next timestep
&lsquo;cut timestep size&lsquo;: reduce the current timestep size by a specified factor and redo the timestep
&lsquo;abort program&lsquo;: abort the program with an error message.

(parameters:Nonlinear_20solver_20scheme)=
### __Parameter name:__ Nonlinear solver scheme
**Default value:** single Advection, single Stokes

**Pattern:** [Selection single Advection, single Stokes|iterated Advection and Stokes|single Advection, iterated Stokes|no Advection, iterated Stokes|no Advection, single Stokes|no Advection, iterated defect correction Stokes|single Advection, iterated defect correction Stokes|iterated Advection and defect correction Stokes|iterated Advection and Newton Stokes|single Advection, iterated Newton Stokes|single Advection, no Stokes|first timestep only, single Stokes|no Advection, no Stokes ]

**Documentation:** The kind of scheme used to resolve the nonlinearity in the system. &lsquo;single Advection, single Stokes&rsquo; means that no nonlinear iterations are done, and the temperature, compositional fields and Stokes equations are solved exactly once per time step, one after the other. The &lsquo;iterated Advection and Stokes&rsquo; scheme iterates this decoupled approach by alternating the solution of the temperature, composition and Stokes systems. The &lsquo;single Advection, iterated Stokes&rsquo; scheme solves the temperature and composition equation once at the beginning of each time step and then iterates out the solution of the Stokes equation. The &lsquo;no Advection, iterated Stokes&rsquo; scheme only solves the Stokes system, iterating out the solution, and ignores compositions and the temperature equation (careful, the material model must not depend on the temperature or composition; this is mostly useful for Stokes benchmarks).  The &lsquo;no Advection, single Stokes&rsquo; scheme only solves the Stokes system once per timestep. This is also mostly useful for Stokes benchmarks. The &lsquo;single Advection, no Stokes&rsquo; scheme only solves the temperature and other advection systems once, and instead of solving for the Stokes system, a prescribed velocity and pressure is used. The &lsquo;iterated Advection and Newton Stokes&rsquo; scheme iterates by alternating the solution of the temperature, composition and Stokes equations, using Picard iterations for the temperature and composition, and Newton iterations for the Stokes system. The &lsquo;single Advection, iterated Newton Stokes&rsquo; scheme solves the temperature and composition equations once at the beginning of each time step and then iterates out the solution of the Stokes equation, using Newton iterations for the Stokes system. The &lsquo;iterated Advection and defect correction Stokes&rsquo; scheme iterates by alternating the solution of the temperature, composition and Stokes equations, using Picard iterations for the temperature and composition, and defect correction Picard iterations for the Stokes system. The &lsquo;single Advection, iterated defect correction Stokes&rsquo; scheme solves the temperature and composition equations once at the beginning of each time step and then iterates out the solution of the Stokes equation, using defect correction Picard iterations for the Stokes system. The &lsquo;no Advection, iterated defect correction Stokes&rsquo; scheme solves the temperature and composition equations once at the beginning of each time step and then iterates out the solution of the Stokes equation, using defect correction Picard iterations for the Stokes system. The &lsquo;first timestep only, single Stokes&rsquo; scheme solves the Stokes equations exactly once, at the first time step. No nonlinear iterations are done, and the temperature and composition systems are not solved.

(parameters:Nonlinear_20solver_20tolerance)=
### __Parameter name:__ Nonlinear solver tolerance
**Default value:** 1e-5

**Pattern:** [Double 0...1 (inclusive)]

**Documentation:** A relative tolerance up to which the nonlinear solver will iterate. This parameter is only relevant if the &lsquo;Nonlinear solver scheme&rsquo; does nonlinear iterations, in other words, if it is set to something other than &lsquo;single Advection, single Stokes&rsquo; or &lsquo;single Advection, no Stokes&rsquo;.

(parameters:Output_20directory)=
### __Parameter name:__ Output directory
**Default value:** output

**Pattern:** [DirectoryName]

**Documentation:** The name of the directory into which all output files should be placed. This may be an absolute or a relative path. ASPECT will write output such as statistics files or visualization files into this directory or into directories further nested within.

(parameters:Output_20directory_20LFS_20stripe_20count)=
### __Parameter name:__ Output directory LFS stripe count
**Default value:** 0

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** Many large clusters use the Lustre file system (LFS) that allows to &rsquo;stripe&rsquo; files, i.e., to use multiple file servers to store a single file. This is useful when writing very large files from multiple MPI processes, such as when creating graphical output or creating checkpoints. In those cases, if all MPI processes try to route their data to a single file server, that file server and the disks it manages may be saturated by data and everything slows down. File striping instead ensures that the data is sent to several file servers, improving performance. A description of how Lustre manages file striping can be found at https://doc.lustre.org/lustre_manual.xhtml#managingstripingfreespace . How file striping can be configured is discussed at https://wiki.lustre.org/Configuring_Lustre_File_Striping .

When this parameter is set to anything other than zero, ASPECT will call the Lustre support tool, &lsquo;lst&lsquo;, as follows: &lsquo;lst setstripe -c N OUTPUT_DIR&lsquo;, where &lsquo;N&lsquo; is the value of the input parameter discussed here, and &lsquo;OUTPUT_DIR&lsquo; is the directory into which ASPECT writes its output. The file striping so set on the output directory are also inherited by the sub-directories ASPECT creates within it.

In order to use this parameter, your cluster must obviously be using the Lustre file system. What the correct value for the stripe count is is something you will have to find out from your cluster&rsquo;s local documentation, or your cluster administrator. It depends on the physical details and configuration of the file servers attached to a cluster.

(parameters:Pressure_20normalization)=
### __Parameter name:__ Pressure normalization
**Default value:** surface

**Pattern:** [Selection surface|volume|no ]

**Documentation:** If and how to normalize the pressure after the solution step. This is necessary because depending on boundary conditions, in many cases the pressure is only determined by the model up to a constant. On the other hand, we often would like to have a well-determined pressure, for example for table lookups of material properties in models or for comparing solutions. If the given value is &lsquo;surface&rsquo;, then normalization at the end of each time steps adds a constant value to the pressure in such a way that the average pressure at the surface of the domain is what is set in the &lsquo;Surface pressure&rsquo; parameter; the surface of the domain is determined by asking the geometry model whether a particular face of the geometry has a zero or small &lsquo;depth&rsquo;. If the value of this parameter is &lsquo;volume&rsquo; then the pressure is normalized so that the domain average is zero. If &lsquo;no&rsquo; is given, the no pressure normalization is performed.

(parameters:Resume_20computation)=
### __Parameter name:__ Resume computation
**Default value:** false

**Pattern:** [Selection true|false|auto ]

**Documentation:** A flag indicating whether the computation should be resumed from a previously saved state (if true) or start from scratch (if false). If auto is selected, models will be resumed if there is an existing checkpoint file, otherwise started from scratch.

(parameters:Start_20time)=
### __Parameter name:__ Start time
**Default value:** 0.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** The start time of the simulation. Units: Years if the &rsquo;Use years in output instead of seconds&rsquo; parameter is set; seconds otherwise.

(parameters:Surface_20pressure)=
### __Parameter name:__ Surface pressure
**Default value:** 0.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** The value the pressure is normalized to in each time step when &lsquo;Pressure normalization&rsquo; is set to &lsquo;surface&rsquo; with default value 0. This setting is ignored in all other cases.

The mathematical equations that describe thermal convection only determine the pressure up to an arbitrary constant. On the other hand, for comparison and for looking up material parameters it is important that the pressure be normalized somehow. We do this by enforcing a particular average pressure value at the surface of the domain, where the geometry model determines where the surface is. This parameter describes what this average surface pressure value is supposed to be. By default, it is set to zero, but one may want to choose a different value for example for simulating only the volume of the mantle below the lithosphere, in which case the surface pressure should be the lithostatic pressure at the bottom of the lithosphere.

For more information, see the section in the manual that discusses the general mathematical model.

(parameters:Timing_20output_20frequency)=
### __Parameter name:__ Timing output frequency
**Default value:** 100

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** How frequently in timesteps to output timing information. This is generally adjusted only for debugging and timing purposes. If the value is set to zero it will also output timing information at the initiation timesteps.

(parameters:Use_20conduction_20timestep)=
### __Parameter name:__ Use conduction timestep
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Mantle convection simulations are often focused on convection dominated systems. However, these codes can also be used to investigate systems where heat conduction plays a dominant role. This parameter indicates whether the simulator should also use heat conduction in determining the length of each time step.

(parameters:Use_20operator_20splitting)=
### __Parameter name:__ Use operator splitting
**Default value:** false

**Pattern:** [Bool]

**Documentation:** If set to true, the advection and reactions of compositional fields and temperature are solved separately, and can use different time steps. Note that this will only work if the material/heating model fills the reaction\_rates/heating\_reaction\_rates structures. Operator splitting can be used with any existing solver schemes that solve the temperature/composition equations.

(parameters:Use_20years_20in_20output_20instead_20of_20seconds)=
### __Parameter name:__ Use years in output instead of seconds
**Default value:** true

**Pattern:** [Bool]

**Documentation:** When computing results for mantle convection simulations, it is often difficult to judge the order of magnitude of results when they are stated in MKS units involving seconds. Rather, some kinds of results such as velocities are often stated in terms of meters per year (or, sometimes, centimeters per year). On the other hand, for non-dimensional computations, one wants results in their natural unit system as used inside the code. If this flag is set to &lsquo;true&rsquo; conversion to years happens; if it is &lsquo;false&rsquo;, no such conversion happens.

Contrary to the word &ldquo;output&rdquo; in the name of this parameter, a number of plugins also use this parameter to determine how to interpret their *inputs*. For example, when &lsquo;true&rsquo;, several of the boundary velocity models described in {ref}`parameters:Boundary_20velocity_20model` interpret both specific times in years instead of seconds, and velocities in meters per year instead of meters per second.

For the purposes of this parameter, a year consists of 60*60*24*365.2425 seconds. In other words, a year is taken to have 365.2425 days.

(parameters:World_20builder_20file)=
### __Parameter name:__ World builder file
**Default value:**

**Pattern:** [FileName (Type: input)]

**Documentation:** Name of the world builder file. If empty, the world builder is not initialized.
