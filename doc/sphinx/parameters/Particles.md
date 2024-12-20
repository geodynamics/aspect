(parameters:Particles)=
# Particles


## **Subsection:** Particles


(parameters:Particles/Allow_20cells_20without_20particles)=
### __Parameter name:__ Allow cells without particles
**Default value:** false

**Pattern:** [Bool]

**Documentation:** By default, every cell needs to contain particles to use this interpolator plugin. If this parameter is set to true, cells are allowed to have no particles. In case both the current cell and its neighbors are empty, the interpolator will return 0 for the current cell&rsquo;s properties.

(parameters:Particles/Integration_20scheme)=
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

(parameters:Particles/Interpolation_20scheme)=
### __Parameter name:__ Interpolation scheme
**Default value:** cell average

**Pattern:** [Selection bilinear least squares|cell average|distance weighted average|harmonic average|nearest neighbor|quadratic least squares ]

**Documentation:** Select one of the following models:

&lsquo;bilinear least squares&rsquo;: Uses linear least squares to obtain the slopes and center of a 2d or 3d plane from the particle positions and a particular property value on those particles. Interpolate this property onto a vector of points. If the limiter is enabled then it will ensure the interpolated properties do not exceed the range of the minimum and maximum of the values of the property on the particles. Note that deal.II must be configured with BLAS and LAPACK to support this operation.

&lsquo;cell average&rsquo;: Return the arithmetic average of all particle properties in the given cell, or in the neighboring cells if the given cell is empty. In case the neighboring cells are also empty, and &rsquo;Allow cells without particles&rsquo; is set to true, the interpolator returns 0. Otherwise, an exception is thrown.

&lsquo;distance weighted average&rsquo;: Interpolates particle properties onto a vector of points using a distance weighed averaging method.

&lsquo;harmonic average&rsquo;: Return the harmonic average of all particle properties in the given cell. If the cell contains no particles, return the harmonic average of the properties in the neighboring cells. In case the neighboring cells are also empty, and &rsquo;Allow cells without particles&rsquo; is set to true, the interpolator returns 0. Otherwise, an exception is thrown.

&lsquo;nearest neighbor&rsquo;: Return the properties of the nearest neighboring particle in the current cell, or nearest particle in nearest neighboring cell if current cell is empty. In case the neighboring cells are also empty, and &rsquo;Allow cells without particles&rsquo; is set to true, the interpolator returns 0. Otherwise, an exception is thrown.

&lsquo;quadratic least squares&rsquo;: Interpolates particle properties onto a vector of points using a quadratic least squares method. Note that deal.II must be configured with BLAS/LAPACK.

(parameters:Particles/List_20of_20particle_20properties)=
### __Parameter name:__ List of particle properties
**Default value:**

**Pattern:** [MultipleSelection composition|cpo bingham average|cpo elastic tensor|crystal preferred orientation|elastic stress|elastic tensor decomposition|function|grain size|initial composition|initial position|integrated strain|integrated strain invariant|melt particle|pT path|position|reference position|strain rate|velocity|viscoplastic strain invariants ]

**Documentation:** A comma separated list of particle properties that should be tracked. By default none is selected, which means only position, velocity and id of the particles are output.

The following properties are available:

&lsquo;composition&rsquo;: Implementation of a plugin in which the particle property is defined by the compositional fields in the model. This can be used to track solid compositionevolution over time.

&lsquo;cpo bingham average&rsquo;: This is a particle property plugin which computes the Bingham average for the Crystal Preferred Orientation particle property plugin so that it can be visualized.

&lsquo;cpo elastic tensor&rsquo;: A plugin in which the particle property tensor is defined as the Voigt average of the elastic tensors of the minerals in the textured rock.Currently only Olivine and Enstatite are supported.

&lsquo;crystal preferred orientation&rsquo;: WARNING: all the CPO plugins are a work in progress and not ready for production use yet. See https://github.com/geodynamics/aspect/issues/3885 for current status and alternatives. The plugin manages and computes the evolution of Lattice/Crystal Preferred Orientations (LPO/CPO) on particles. Each ASPECT particle can be assigned many grains. Each grain is assigned a size and a orientation matrix. This allows for CPO evolution tracking with polycrystalline kinematic CrystalPreferredOrientation evolution models such as D-Rex (Kaminski and Ribe, 2001; Kaminski et al., 2004).

&lsquo;elastic stress&rsquo;: A plugin in which the particle property tensor is defined as the total elastic stress a particle has accumulated. See the viscoelastic material model documentation for more detailed information.

&lsquo;elastic tensor decomposition&rsquo;: A plugin which decomposes the elastic tensor into different approximations (Isotropic, Hexagonal, Tetragonal, Orthorhombic, Monoclinic and Triclinic) and provides the eigenvectors of the tensor.

&lsquo;function&rsquo;: Implementation of a model in which the particle property is set by evaluating an explicit function at the initial position of each particle. The function is defined in the parameters in section &ldquo;Particles|Function&rdquo;. The format of these functions follows the syntax understood by the muparser library, see {ref}`sec:run-aspect:parameters-overview:muparser-format`.

&lsquo;grain size&rsquo;: A plugin in which the particle property is defined as the evolving grain size of a particle. See the grain_size material model documentation for more detailed information.

&lsquo;initial composition&rsquo;: Implementation of a plugin in which the particle property is given as the initial composition at the particle&rsquo;s initial position. The particle gets as many properties as there are compositional fields.

&lsquo;initial position&rsquo;: Implementation of a plugin in which the particle property is given as the initial position of the particle. This property is vector-valued with as many components as there are space dimensions. In practice, it is often most useful to only visualize one of the components of this vector, or the magnitude of the vector. For example, in a spherical mantle simulation, the magnitude of this property equals the starting radius of a particle, and is thereby indicative of which part of the mantle a particle comes from.

&lsquo;integrated strain&rsquo;: A plugin in which the particle property tensor is defined as the deformation gradient tensor $\mathbf F$ this particle has experienced. $\mathbf F$ can be polar-decomposed into the left stretching tensor $\mathbf L$ (the finite strain we are interested in), and the rotation tensor $\mathbf Q$. See the corresponding cookbook in the manual for more detailed information.

&lsquo;integrated strain invariant&rsquo;: A plugin in which the particle property is defined as the finite strain invariant ($\varepsilon_{ii}$). This property is calculated with the timestep ($dt$) and the second invariant of the deviatoric strain rate tensor ($\dot{\varepsilon}_{ii}$), where the value at time step $n$ is $\varepsilon_{ii}^{n} = \varepsilon_{ii}^{n-1} + dt\dot{\varepsilon}_{ii}$.

&lsquo;melt particle&rsquo;: Implementation of a plugin in which the particle property is defined as presence of melt above a threshold, which can be set as an input parameter. This property is set to 0 if melt is not present and set to 1 if melt is present.

&lsquo;pT path&rsquo;: Implementation of a plugin in which the particle property is defined as the current pressure and temperature at this position. This can be used to generate pressure-temperature paths of material points over time.

&lsquo;position&rsquo;: Implementation of a plugin in which the particle property is defined as the current position.

&lsquo;reference position&rsquo;: Implementation of a plugin in which the particle property is defined as the current reference position.

&lsquo;strain rate&rsquo;: Implementation of a plugin in which the time evolution of strain rate is saved and stored on the particles.

&lsquo;velocity&rsquo;: Implementation of a plugin in which the particle property is defined as the recent velocity at this position.

&lsquo;viscoplastic strain invariants&rsquo;: A plugin that calculates the finite strain invariant a particle has experienced and assigns it to either the plastic and/or viscous strain field based on whether the material is plastically yielding, or the total strain field used in the visco plastic material model. The implementation of this property is equivalent to the implementation for compositional fields that is located in the plugin in `benchmarks/buiter\_et\_al\_2008\_jgr/plugin/`,and is effectively the same as what the visco plastic material model uses for compositional fields.

(parameters:Particles/Load_20balancing_20strategy)=
### __Parameter name:__ Load balancing strategy
**Default value:** repartition

**Pattern:** [MultipleSelection none|remove particles|add particles|remove and add particles|repartition ]

**Documentation:** Strategy that is used to balance the computational load across processors for adaptive meshes.

(parameters:Particles/Maximum_20particles_20per_20cell)=
### __Parameter name:__ Maximum particles per cell
**Default value:** 100

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** Upper limit for particle number per cell. This limit is useful for adaptive meshes to prevent coarse cells from slowing down the whole model. It will be checked and enforced after mesh refinement, after MPI transfer of particles and after particle movement. If there are `n\_number\_of\_particles` $>$ `max\_particles\_per\_cell` particles in one cell then `n\_number\_of\_particles` - `max\_particles\_per\_cell` particles in this cell are randomly chosen and destroyed.

(parameters:Particles/Minimum_20particles_20per_20cell)=
### __Parameter name:__ Minimum particles per cell
**Default value:** 0

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** Lower limit for particle number per cell. This limit is useful for adaptive meshes to prevent fine cells from being empty of particles. It will be checked and enforced after mesh refinement and after particle movement. If there are `n\_number\_of\_particles` $<$ `min\_particles\_per\_cell` particles in one cell then `min\_particles\_per\_cell` - `n\_number\_of\_particles` particles are generated and randomly placed in this cell. If the particles carry properties the individual property plugins control how the properties of the new particles are initialized.

(parameters:Particles/Number_20of_20particle_20systems)=
### __Parameter name:__ Number of particle systems
**Default value:** 1

**Pattern:** [Integer range 0...2 (inclusive)]

**Documentation:** The number of particle systems to be created. The maximum number of particle systems is set by the CMake variable &lsquo;ASPECT_MAX_NUM_PARTICLE_SYSTEMS&lsquo; and is by default 2.

(parameters:Particles/Particle_20generator_20name)=
### __Parameter name:__ Particle generator name
**Default value:** random uniform

**Pattern:** [Selection ascii file|probability density function|quadrature points|random uniform|reference cell|uniform box|uniform radial ]

**Documentation:** Select one of the following models:

&lsquo;ascii file&rsquo;: Generates a distribution of particles from coordinates specified in an Ascii data file. The file format is a simple text file, with as many columns as spatial dimensions and as many lines as particles to be generated. Initial comment lines starting with &lsquo;#&rsquo; will be discarded. Note that this plugin always generates as many particles as there are coordinates in the data file, the &ldquo;Particles/Number of particles&rdquo; parameter has no effect on this plugin. All of the values that define this generator are read from a section &ldquo;Particles/Generator/Ascii file&rdquo; in the input file, see {ref}`parameters:Particles/Generator/Ascii_20file`.

&lsquo;probability density function&rsquo;: Generate a random distribution of particles over the entire simulation domain. The probability density is prescribed in the form of a user-prescribed function. The format of this function follows the syntax understood by the muparser library, see {ref}`sec:run-aspect:parameters-overview:muparser-format`. The return value of the function is always checked to be a non-negative probability density but it can be zero in parts of the domain.

&lsquo;quadrature points&rsquo;: Generates particles at the quadrature points of each active cell of the triangulation. Here, Gauss quadrature of degree (velocity\_degree + 1), is used similarly to the assembly of Stokes matrix.

&lsquo;random uniform&rsquo;: Generates a random uniform distribution of particles over the entire simulation domain. This generator can be understood as the special case of the &rsquo;probability density function&rsquo; generator where the probability density is constant over the domain.

&lsquo;reference cell&rsquo;: Generates a uniform distribution of particles per cell and spatial direction in the unit cell and transforms each of the particles back to real region in the model domain. Uniform here means the particles will be generated with an equal spacing in each spatial dimension.

&lsquo;uniform box&rsquo;: Generate a uniform distribution of particles over a rectangular domain in 2d or 3d. Uniform here means the particles will be generated with an equal spacing in each spatial dimension. Note that in order to produce a regular distribution the number of generated particles might not exactly match the one specified in the input file.

&lsquo;uniform radial&rsquo;: Generate a uniform distribution of particles over a spherical domain in 2d or 3d. Uniform here means the particles will be generated with an equal spacing in each spherical spatial dimension, i.e., the particles are created at positions that increase linearly with equal spacing in radius, colatitude and longitude around a certain center point. Note that in order to produce a regular distribution the number of generated particles might not exactly match the one specified in the input file.

(parameters:Particles/Particle_20weight)=
### __Parameter name:__ Particle weight
**Default value:** 10

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** Weight that is associated with the computational load of a single particle. The sum of particle weights will be added to the sum of cell weights to determine the partitioning of the mesh if the &lsquo;repartition&rsquo; particle load balancing strategy is selected. The optimal weight depends on the used integrator and particle properties. In general for a more expensive integrator and more expensive properties a larger particle weight is recommended. Before adding the weights of particles, each cell already carries a weight of 1000 to account for the cost of field-based computations.

(parameters:Particles/Update_20ghost_20particles)=
### __Parameter name:__ Update ghost particles
**Default value:** true

**Pattern:** [Bool]

**Documentation:** Some particle interpolation algorithms require knowledge about particles in neighboring cells. To allow this, particles in ghost cells need to be exchanged between the processes neighboring this cell. This parameter determines whether this transport is happening. This parameter is deprecated and will be removed in the future. Ghost particle updates are always performed. Please set the parameter to &lsquo;true&rsquo;.

(parameters:Particles/CPO_20Bingham_20Average)=
## **Subsection:** Particles / CPO Bingham Average
(parameters:Particles/CPO_20Bingham_20Average/Number_20of_20samples)=
### __Parameter name:__ Number of samples
**Default value:** 0

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** This determines how many samples are taken when using the random draw volume averaging. Setting it to zero means that the number of samples is set to be equal to the number of grains.

(parameters:Particles/CPO_20Bingham_20Average/Random_20number_20seed)=
### __Parameter name:__ Random number seed
**Default value:** 1

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** The seed used to generate random numbers. This will make sure that results are reproducible as long as the problem is run with the same amount of MPI processes. It is implemented as final seed = Random number seed + MPI Rank.

(parameters:Particles/Crystal_20Preferred_20Orientation)=
## **Subsection:** Particles / Crystal Preferred Orientation
(parameters:Particles/Crystal_20Preferred_20Orientation/CPO_20derivatives_20algorithm)=
### __Parameter name:__ CPO derivatives algorithm
**Default value:** Spin tensor

**Pattern:** [List of <[Anything]> of length 0...4294967295 (inclusive)]

**Documentation:** Options: Spin tensor

(parameters:Particles/Crystal_20Preferred_20Orientation/Number_20of_20grains_20per_20particle)=
### __Parameter name:__ Number of grains per particle
**Default value:** 50

**Pattern:** [Integer range 1...2147483647 (inclusive)]

**Documentation:** The number of grains of each different mineral each particle contains.

(parameters:Particles/Crystal_20Preferred_20Orientation/Property_20advection_20max_20iterations)=
### __Parameter name:__ Property advection max iterations
**Default value:** 100

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** The Backward Euler property advection method involve internal iterations. This option allows for setting the maximum number of iterations. Note that when the iteration is ended by the max iteration amount an assert is thrown.

(parameters:Particles/Crystal_20Preferred_20Orientation/Property_20advection_20method)=
### __Parameter name:__ Property advection method
**Default value:** Backward Euler

**Pattern:** [Anything]

**Documentation:** Options: Forward Euler, Backward Euler

(parameters:Particles/Crystal_20Preferred_20Orientation/Property_20advection_20tolerance)=
### __Parameter name:__ Property advection tolerance
**Default value:** 1e-10

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The Backward Euler property advection method involve internal iterations. This option allows for setting a tolerance. When the norm of tensor new - tensor old is smaller than this tolerance, the iteration is stopped.

(parameters:Particles/Crystal_20Preferred_20Orientation/Random_20number_20seed)=
### __Parameter name:__ Random number seed
**Default value:** 1

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** The seed used to generate random numbers. This will make sure that results are reproducible as long as the problem is run with the same number of MPI processes. It is implemented as final seed = user seed + MPI Rank.

(parameters:Particles/Crystal_20Preferred_20Orientation/D_2dRex_202004)=
## **Subsection:** Particles / Crystal Preferred Orientation / D-Rex 2004
(parameters:Particles/Crystal_20Preferred_20Orientation/D_2dRex_202004/Exponents_20p)=
### __Parameter name:__ Exponents p
**Default value:** 1.5

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** This is exponent p as defined in equation 11 of Kaminski et al., 2004.

(parameters:Particles/Crystal_20Preferred_20Orientation/D_2dRex_202004/Mobility)=
### __Parameter name:__ Mobility
**Default value:** 50

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The dimensionless intrinsic grain boundary mobility for both olivine and enstatite.

(parameters:Particles/Crystal_20Preferred_20Orientation/D_2dRex_202004/Nucleation_20efficiency)=
### __Parameter name:__ Nucleation efficiency
**Default value:** 5

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** This is the dimensionless nucleation rate as defined in equation 8 of Kaminski et al., 2004.

(parameters:Particles/Crystal_20Preferred_20Orientation/D_2dRex_202004/Stress_20exponents)=
### __Parameter name:__ Stress exponents
**Default value:** 3.5

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** This is the power law exponent that characterizes the rheology of the slip systems. It is used in equation 11 of Kaminski et al., 2004.

(parameters:Particles/Crystal_20Preferred_20Orientation/D_2dRex_202004/Threshold_20GBS)=
### __Parameter name:__ Threshold GBS
**Default value:** 0.3

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The Dimensionless Grain Boundary Sliding (GBS) threshold. This is a grain size threshold below which grain deform by GBS and become strain-free grains.

(parameters:Particles/Crystal_20Preferred_20Orientation/D_2dRex_202004/Volume_20fractions_20minerals)=
### __Parameter name:__ Volume fractions minerals
**Default value:** 0.5, 0.5

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** The volume fraction for the different minerals. There need to be the same amount of values as there are minerals

(parameters:Particles/Crystal_20Preferred_20Orientation/Initial_20grains)=
## **Subsection:** Particles / Crystal Preferred Orientation / Initial grains
(parameters:Particles/Crystal_20Preferred_20Orientation/Initial_20grains/Minerals)=
### __Parameter name:__ Minerals
**Default value:** Olivine: Karato 2008, Enstatite

**Pattern:** [List of <[Anything]> of length 0...4294967295 (inclusive)]

**Documentation:** This determines what minerals and fabrics or fabric selectors are used used for the LPO/CPO calculation. The options are Olivine: Passive, A-fabric, Olivine: B-fabric, Olivine: C-fabric, Olivine: D-fabric, Olivine: E-fabric, Olivine: Karato 2008 or Enstatite. Passive sets all RRSS entries to the maximum. The Karato 2008 selector selects a fabric based on stress and water content as defined in figure 4 of the Karato 2008 review paper (doi: 10.1146/annurev.earth.36.031207.124120).

(parameters:Particles/Crystal_20Preferred_20Orientation/Initial_20grains/Model_20name)=
### __Parameter name:__ Model name
**Default value:** Uniform grains and random uniform rotations

**Pattern:** [Anything]

**Documentation:** The model used to initialize the CPO for all particles. Currently &rsquo;Uniform grains and random uniform rotations&rsquo; and &rsquo;World Builder&rsquo; are the only valid option.

(parameters:Particles/Crystal_20Preferred_20Orientation/Initial_20grains/Volume_20fractions_20minerals)=
### __Parameter name:__ Volume fractions minerals
**Default value:** 0.7, 0.3

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** The volume fractions for the different minerals. There need to be the same number of values as there are minerals.Note that the currently implemented scheme is incompressible and does not allow chemical interaction or the formation of new phases

(parameters:Particles/Function)=
## **Subsection:** Particles / Function
(parameters:Particles/Function/Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**

**Pattern:** [Anything]

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form &lsquo;var1=value1, var2=value2, ...&rsquo;.

A typical example would be to set this runtime parameter to &lsquo;pi=3.1415926536&rsquo; and then use &lsquo;pi&rsquo; in the expression of the actual formula. (That said, for convenience this class actually defines both &lsquo;pi&rsquo; and &lsquo;Pi&rsquo; by default, but you get the idea.)

(parameters:Particles/Function/Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0

**Pattern:** [Anything]

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as &lsquo;sin&rsquo; or &lsquo;cos&rsquo;. In addition, it may contain expressions like &lsquo;if(x>0, 1, -1)&rsquo; where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon.

(parameters:Particles/Function/Number_20of_20components)=
### __Parameter name:__ Number of components
**Default value:** 1

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** The number of function components where each component is described by a function expression delimited by a &rsquo;;&rsquo;.

(parameters:Particles/Function/Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,y,t

**Pattern:** [Anything]

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are &lsquo;x&rsquo; (in 1d), &lsquo;x,y&rsquo; (in 2d) or &lsquo;x,y,z&rsquo; (in 3d) for spatial coordinates and &lsquo;t&rsquo; for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to &lsquo;r,phi,theta,t&rsquo; and then use these variable names in your function expression.

(parameters:Particles/Generator)=
## **Subsection:** Particles / Generator
(parameters:Particles/Generator/Ascii_20file)=
## **Subsection:** Particles / Generator / Ascii file
(parameters:Particles/Generator/Ascii_20file/Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/particle/generator/ascii/

**Pattern:** [DirectoryName]

**Documentation:** The name of a directory that contains the particle data. This path may either be absolute (if starting with a &rsquo;/&rsquo;) or relative to the current directory. The path may also include the special text &rsquo;$ASPECT_SOURCE_DIR&rsquo; which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the &lsquo;data/&rsquo; subdirectory of ASPECT.

(parameters:Particles/Generator/Ascii_20file/Data_20file_20name)=
### __Parameter name:__ Data file name
**Default value:** particle.dat

**Pattern:** [Anything]

**Documentation:** The name of the particle file.

(parameters:Particles/Generator/Probability_20density_20function)=
## **Subsection:** Particles / Generator / Probability density function
(parameters:Particles/Generator/Probability_20density_20function/Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**

**Pattern:** [Anything]

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form &lsquo;var1=value1, var2=value2, ...&rsquo;.

A typical example would be to set this runtime parameter to &lsquo;pi=3.1415926536&rsquo; and then use &lsquo;pi&rsquo; in the expression of the actual formula. (That said, for convenience this class actually defines both &lsquo;pi&rsquo; and &lsquo;Pi&rsquo; by default, but you get the idea.)

(parameters:Particles/Generator/Probability_20density_20function/Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 1.0

**Pattern:** [Anything]

**Documentation:** The formula that denotes the spatially variable probability density function. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as &lsquo;sin&rsquo; or &lsquo;cos&rsquo;. In addition, it may contain expressions like &lsquo;if(x>0, 1, 0)&rsquo; where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise; this example would result in no particles at all in that part of the domain where $x==0$, and a constant particle density in the rest of the domain. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/. Note that the function has to be non-negative everywhere in the domain, and needs to be positive in at least some parts of the domain.

(parameters:Particles/Generator/Probability_20density_20function/Number_20of_20particles)=
### __Parameter name:__ Number of particles
**Default value:** 1000

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Total number of particles to create (not per processor or per element). The number is parsed as a floating point number (so that one can specify, for example, &rsquo;1e4&rsquo; particles) but it is interpreted as an integer, of course.

(parameters:Particles/Generator/Probability_20density_20function/Random_20cell_20selection)=
### __Parameter name:__ Random cell selection
**Default value:** true

**Pattern:** [Bool]

**Documentation:** If true, particle numbers per cell are calculated randomly according to their respective probability density. This means particle numbers per cell can deviate statistically from the integral of the probability density. If false, first determine how many particles each cell should have based on the integral of the density over each of the cells, and then once we know how many particles we want on each cell, choose their locations randomly within each cell.

(parameters:Particles/Generator/Probability_20density_20function/Random_20number_20seed)=
### __Parameter name:__ Random number seed
**Default value:** 5432

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** The seed for the random number generator that controls the particle generation. Keep constant to generate identical particle distributions in subsequent model runs. Change to get a different distribution. In parallel computations the seed is further modified on each process to ensure different particle patterns on different processes. Note that the number of particles per processor is not affected by the seed.

(parameters:Particles/Generator/Probability_20density_20function/Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,y,t

**Pattern:** [Anything]

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are &lsquo;x&rsquo; (in 1d), &lsquo;x,y&rsquo; (in 2d) or &lsquo;x,y,z&rsquo; (in 3d) for spatial coordinates and &lsquo;t&rsquo; for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to &lsquo;r,phi,theta,t&rsquo; and then use these variable names in your function expression.

(parameters:Particles/Generator/Random_20uniform)=
## **Subsection:** Particles / Generator / Random uniform
(parameters:Particles/Generator/Random_20uniform/Number_20of_20particles)=
### __Parameter name:__ Number of particles
**Default value:** 1000

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Total number of particles to create (not per processor or per element). The number is parsed as a floating point number (so that one can specify, for example, &rsquo;1e4&rsquo; particles) but it is interpreted as an integer, of course.

(parameters:Particles/Generator/Random_20uniform/Random_20cell_20selection)=
### __Parameter name:__ Random cell selection
**Default value:** true

**Pattern:** [Bool]

**Documentation:** If true, particle numbers per cell are calculated randomly according to their respective probability density. This means particle numbers per cell can deviate statistically from the integral of the probability density. If false, first determine how many particles each cell should have based on the integral of the density over each of the cells, and then once we know how many particles we want on each cell, choose their locations randomly within each cell.

(parameters:Particles/Generator/Random_20uniform/Random_20number_20seed)=
### __Parameter name:__ Random number seed
**Default value:** 5432

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** The seed for the random number generator that controls the particle generation. Keep constant to generate identical particle distributions in subsequent model runs. Change to get a different distribution. In parallel computations the seed is further modified on each process to ensure different particle patterns on different processes. Note that the number of particles per processor is not affected by the seed.

(parameters:Particles/Generator/Reference_20cell)=
## **Subsection:** Particles / Generator / Reference cell
(parameters:Particles/Generator/Reference_20cell/Number_20of_20particles_20per_20cell_20per_20direction)=
### __Parameter name:__ Number of particles per cell per direction
**Default value:** 2

**Pattern:** [List of <[Integer range 1...2147483647 (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** List of number of particles to create per cell and spatial dimension. The size of the list is the number of spatial dimensions. If only one value is given, then each spatial dimension is set to the same value. The list of numbers are parsed as a floating point number (so that one can specify, for example, &rsquo;1e4&rsquo; particles) but it is interpreted as an integer, of course.

(parameters:Particles/Generator/Uniform_20box)=
## **Subsection:** Particles / Generator / Uniform box
(parameters:Particles/Generator/Uniform_20box/Maximum_20x)=
### __Parameter name:__ Maximum x
**Default value:** 1.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Maximum x coordinate for the region of particles.

(parameters:Particles/Generator/Uniform_20box/Maximum_20y)=
### __Parameter name:__ Maximum y
**Default value:** 1.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Maximum y coordinate for the region of particles.

(parameters:Particles/Generator/Uniform_20box/Maximum_20z)=
### __Parameter name:__ Maximum z
**Default value:** 1.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Maximum z coordinate for the region of particles.

(parameters:Particles/Generator/Uniform_20box/Minimum_20x)=
### __Parameter name:__ Minimum x
**Default value:** 0.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Minimum x coordinate for the region of particles.

(parameters:Particles/Generator/Uniform_20box/Minimum_20y)=
### __Parameter name:__ Minimum y
**Default value:** 0.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Minimum y coordinate for the region of particles.

(parameters:Particles/Generator/Uniform_20box/Minimum_20z)=
### __Parameter name:__ Minimum z
**Default value:** 0.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Minimum z coordinate for the region of particles.

(parameters:Particles/Generator/Uniform_20box/Number_20of_20particles)=
### __Parameter name:__ Number of particles
**Default value:** 1000

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Total number of particles to create (not per processor or per element). The number is parsed as a floating point number (so that one can specify, for example, &rsquo;1e4&rsquo; particles) but it is interpreted as an integer, of course.

(parameters:Particles/Generator/Uniform_20radial)=
## **Subsection:** Particles / Generator / Uniform radial
(parameters:Particles/Generator/Uniform_20radial/Center_20x)=
### __Parameter name:__ Center x
**Default value:** 0.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** x coordinate for the center of the spherical region, where particles are generated.

(parameters:Particles/Generator/Uniform_20radial/Center_20y)=
### __Parameter name:__ Center y
**Default value:** 0.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** y coordinate for the center of the spherical region, where particles are generated.

(parameters:Particles/Generator/Uniform_20radial/Center_20z)=
### __Parameter name:__ Center z
**Default value:** 0.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** z coordinate for the center of the spherical region, where particles are generated.

(parameters:Particles/Generator/Uniform_20radial/Maximum_20latitude)=
### __Parameter name:__ Maximum latitude
**Default value:** 180.

**Pattern:** [Double 0...180 (inclusive)]

**Documentation:** Maximum latitude coordinate for the region of particles in degrees. Measured from the center position, and from the north pole.

(parameters:Particles/Generator/Uniform_20radial/Maximum_20longitude)=
### __Parameter name:__ Maximum longitude
**Default value:** 360.

**Pattern:** [Double -180...360 (inclusive)]

**Documentation:** Maximum longitude coordinate for the region of particles in degrees. Measured from the center position.

(parameters:Particles/Generator/Uniform_20radial/Maximum_20radius)=
### __Parameter name:__ Maximum radius
**Default value:** 1.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Maximum radial coordinate for the region of particles. Measured from the center position.

(parameters:Particles/Generator/Uniform_20radial/Minimum_20latitude)=
### __Parameter name:__ Minimum latitude
**Default value:** 0.

**Pattern:** [Double 0...180 (inclusive)]

**Documentation:** Minimum latitude coordinate for the region of particles in degrees. Measured from the center position, and from the north pole.

(parameters:Particles/Generator/Uniform_20radial/Minimum_20longitude)=
### __Parameter name:__ Minimum longitude
**Default value:** 0.

**Pattern:** [Double -180...360 (inclusive)]

**Documentation:** Minimum longitude coordinate for the region of particles in degrees. Measured from the center position.

(parameters:Particles/Generator/Uniform_20radial/Minimum_20radius)=
### __Parameter name:__ Minimum radius
**Default value:** 0.

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Minimum radial coordinate for the region of particles. Measured from the center position.

(parameters:Particles/Generator/Uniform_20radial/Number_20of_20particles)=
### __Parameter name:__ Number of particles
**Default value:** 1000

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Total number of particles to create (not per processor or per element). The number is parsed as a floating point number (so that one can specify, for example, &rsquo;1e4&rsquo; particles) but it is interpreted as an integer, of course.

(parameters:Particles/Generator/Uniform_20radial/Radial_20layers)=
### __Parameter name:__ Radial layers
**Default value:** 1

**Pattern:** [Integer range 1...2147483647 (inclusive)]

**Documentation:** The number of radial shells of particles that will be generated around the central point.

(parameters:Particles/Integrator)=
## **Subsection:** Particles / Integrator
(parameters:Particles/Integrator/RK2)=
## **Subsection:** Particles / Integrator / RK2
(parameters:Particles/Integrator/RK2/Higher_20order_20accurate_20in_20time)=
### __Parameter name:__ Higher order accurate in time
**Default value:** true

**Pattern:** [Bool]

**Documentation:** Whether to correctly evaluate old and current velocity solution to reach higher-order accuracy in time. If set to &rsquo;false&rsquo; only the old velocity solution is evaluated to simulate a first order method in time. This is only recommended for benchmark purposes.

(parameters:Particles/Interpolator)=
## **Subsection:** Particles / Interpolator
(parameters:Particles/Interpolator/Bilinear_20least_20squares)=
## **Subsection:** Particles / Interpolator / Bilinear least squares
(parameters:Particles/Interpolator/Bilinear_20least_20squares/Use_20boundary_20extrapolation)=
### __Parameter name:__ Use boundary extrapolation
**Default value:** false

**Pattern:** [List of <[Bool]> of length 0...4294967295 (inclusive)]

**Documentation:** Extends the range used by &rsquo;Use linear least squares limiter&rsquo; by linearly interpolating values at cell boundaries from neighboring cells. If more than one value is given, it will be treated as a list with one component per particle property. Enabling &rsquo;Use boundary extrapolation&rsquo; requires enabling &rsquo;Use linear least squares limiter&rsquo;.

(parameters:Particles/Interpolator/Bilinear_20least_20squares/Use_20linear_20least_20squares_20limiter)=
### __Parameter name:__ Use linear least squares limiter
**Default value:** true

**Pattern:** [List of <[Bool]> of length 0...4294967295 (inclusive)]

**Documentation:** Limit the interpolation of particle properties onto the cell, so that the value of each property is no smaller than its minimum and no larger than its maximum on the particles of each cell, and the average of neighboring cells. If more than one value is given, it will be treated as a list with one component per particle property.

(parameters:Particles/Interpolator/Quadratic_20least_20squares)=
## **Subsection:** Particles / Interpolator / Quadratic least squares
(parameters:Particles/Interpolator/Quadratic_20least_20squares/Use_20boundary_20extrapolation)=
### __Parameter name:__ Use boundary extrapolation
**Default value:** false

**Pattern:** [List of <[Bool]> of length 0...4294967295 (inclusive)]

**Documentation:** Extends the range used by &rsquo;Use quadratic least squares limiter&rsquo; by linearly interpolating values at cell boundaries from neighboring cells. If more than one value is given, it will be treated as a list with one component per particle property. Enabling &rsquo;Use boundary extrapolation&rsquo; requires enabling &rsquo;Use quadratic least squares limiter&rsquo;.

(parameters:Particles/Interpolator/Quadratic_20least_20squares/Use_20quadratic_20least_20squares_20limiter)=
### __Parameter name:__ Use quadratic least squares limiter
**Default value:** true

**Pattern:** [List of <[Bool]> of length 0...4294967295 (inclusive)]

**Documentation:** Limit the interpolation of particle properties onto the cell, so that the value of each property is no smaller than its minimum and no larger than its maximum on the particles of each cell, and the average of neighboring cells. If more than one value is given, it will be treated as a list with one component per particle property.

(parameters:Particles/Melt_20particle)=
## **Subsection:** Particles / Melt particle
(parameters:Particles/Melt_20particle/Threshold_20for_20melt_20presence)=
### __Parameter name:__ Threshold for melt presence
**Default value:** 1e-3

**Pattern:** [Double 0...1 (inclusive)]

**Documentation:** The minimum porosity that has to be present at the position of a particle for it to be considered a melt particle (in the sense that the melt presence property is set to 1).
