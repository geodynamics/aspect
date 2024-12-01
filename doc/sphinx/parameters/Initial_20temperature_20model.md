(parameters:Initial_20temperature_20model)=
# Initial temperature model


## **Subsection:** Initial temperature model


(parameters:Initial_20temperature_20model/List_20of_20model_20names)=
### __Parameter name:__ List of model names
**Default value:**

**Pattern:** [MultipleSelection S40RTS perturbation|SAVANI perturbation|adiabatic|adiabatic boundary|ascii data|ascii data layered|ascii profile|continental geotherm|function|harmonic perturbation|inclusion shape perturbation|lithosphere mask|mandelbox|patch on S40RTS|perturbed box|polar box|prescribed temperature|random Gaussian perturbation|spherical gaussian perturbation|spherical hexagonal perturbation|world builder ]

**Documentation:** A comma-separated list of initial temperature models that will be used to initialize the temperature. These plugins are loaded in the order given, and modify the existing temperature field via the operators listed in &rsquo;List of model operators&rsquo;.

The following initial temperature models are available:

&lsquo;S40RTS perturbation&rsquo;: An initial temperature field in which the temperature is perturbed following the S20RTS or S40RTS shear wave velocity model by Ritsema and others, which can be downloaded here \url{http://www.earth.lsa.umich.edu/~jritsema/research.html}. Information on the vs model can be found in Ritsema, J., Deuss, A., van Heijst, H.J. \& Woodhouse, J.H., 2011. S40RTS: a degree-40 shear-velocity model for the mantle from new Rayleigh wave dispersion, teleseismic traveltime and normal-mode splitting function measurements, Geophys. J. Int. 184, 1223-1236. The scaling between the shear wave perturbation and the density perturbation can be constant and set by the user with the &rsquo;Vs to density scaling&rsquo; parameter or depth-dependent and read in from a file. To convert density the user can specify the &rsquo;Thermal expansion coefficient in initial temperature scaling&rsquo; parameter. The scaling is as follows: $\delta \ln \rho (r,\theta,\phi) = \xi \cdot \delta \ln v_s(r,\theta, \phi)$ and $\delta T(r,\theta,\phi) = - \frac{1}{\alpha} \delta \ln \rho(r,\theta,\phi)$. $\xi$ is the &lsquo;vs to density scaling&rsquo; parameter and $\alpha$ is the &rsquo;Thermal expansion coefficient in initial temperature scaling&rsquo; parameter. The temperature perturbation is added to an otherwise constant temperature (incompressible model) or adiabatic reference profile (compressible model). If a depth is specified in &rsquo;Remove temperature heterogeneity down to specified depth&rsquo;, there is no temperature perturbation prescribed down to that depth.
Note the required file format if the vs to density scaling is read in from a file: The first lines may contain any number of comments if they begin with &rsquo;#&rsquo;, but one of these lines needs to contain the number of points in the reference state as for example &rsquo;# POINTS: 3&rsquo;. Following the comment lines there has to be a single line containing the names of all data columns, separated by arbitrarily many spaces. Column names are not allowed to contain spaces. The file can contain unnecessary columns, but for this plugin it needs to at least provide the columns named &lsquo;depth&rsquo; and &lsquo;vs\_to\_density&rsquo;. Note that the data lines in the file need to be sorted in order of increasing depth from 0 to the maximal depth in the model domain. Points in the model that are outside of the provided depth range will be assigned the maximum or minimum depth values, respectively. Points do not need to be equidistant, but the computation of properties is optimized in speed if they are.
If the plugin is used in 2d it will use an equatorial slice of the seismic tomography model.

&lsquo;SAVANI perturbation&rsquo;: An initial temperature field in which the temperature is perturbed following the SAVANI shear wave velocity model by Auer and others, which can be downloaded here \url{http://n.ethz.ch/~auerl/savani.tar.bz2}. Information on the vs model can be found in Auer, L., Boschi, L., Becker, T.W., Nissen-Meyer, T. \& Giardini, D., 2014. Savani: A variable resolution whole-mantle model of anisotropic shear velocity variations based on multiple data sets. Journal of Geophysical Research: Solid Earth 119.4 (2014): 3006-3034. The scaling between the shear wave perturbation and the density perturbation can be constant and set by the user with the &rsquo;Vs to density scaling&rsquo; parameter or depth-dependent and read in from a file. To convert density the user can specify the &rsquo;Thermal expansion coefficient in initial temperature scaling&rsquo; parameter. The scaling is as follows: $\delta \ln \rho (r,\theta,\phi) = \xi \cdot \delta \ln v_s(r,\theta, \phi)$ and $\delta T(r,\theta,\phi) = - \frac{1}{\alpha} \delta \ln \rho(r,\theta,\phi)$. $\xi$ is the &lsquo;vs to density scaling&rsquo; parameter and $\alpha$ is the &rsquo;Thermal expansion coefficient in initial temperature scaling&rsquo; parameter. The temperature perturbation is added to an otherwise constant temperature (incompressible model) or adiabatic reference profile (compressible model).If a depth is specified in &rsquo;Remove temperature heterogeneity down to specified depth&rsquo;, there is no temperature perturbation prescribed down to that depth.
Note the required file format if the vs to density scaling is read in from a file: The first lines may contain any number of comments if they begin with &rsquo;#&rsquo;, but one of these lines needs to contain the number of points in the reference state as for example &rsquo;# POINTS: 3&rsquo;. Following the comment lines there has to be a single line containing the names of all data columns, separated by arbitrarily many spaces. Column names are not allowed to contain spaces. The file can contain unnecessary columns, but for this plugin it needs to at least provide the columns named &lsquo;depth&rsquo; and &lsquo;vs\_to\_density&rsquo;. Note that the data lines in the file need to be sorted in order of increasing depth from 0 to the maximal depth in the model domain. Points in the model that are outside of the provided depth range will be assigned the maximum or minimum depth values, respectively. Points do not need to be equidistant, but the computation of properties is optimized in speed if they are.

&lsquo;adiabatic&rsquo;: Temperature is prescribed as an adiabatic profile with upper and lower thermal boundary layers, whose ages are given as input parameters. Note that this plugin uses the &rsquo;Adiabatic conditions model&rsquo; to compute the adiabat. Thus, the results depend on variables defined outside of this specific subsection; e.g. the globally defined &rsquo;Adiabatic surface temperature&rsquo;, and the variables defined in the &rsquo;Material model&rsquo; section including densities, heat capacities and thermal expansivities.

&lsquo;adiabatic boundary&rsquo;: An initial temperature condition that allows for discretizing the model domain into two layers separated by a user-defined isothermal boundary. The user includes an input ascii data file that is formatted as 3 columns of &lsquo;longitude(radians)&rsquo;, &lsquo;colatitude(radians)&rsquo;, and &lsquo;isotherm depth(meters)&rsquo;, where &lsquo;isotherm depth&rsquo; represents the depth of an initial temperature of 1673.15 K (by default). The first lines in the data file may contain any number of comments if they begin with &lsquo;#&rsquo;, but one of these lines needs to contain the number of grid points in each dimension as for example &lsquo;# POINTS: 69 121&rsquo;. Note that the coordinates need to be sorted in a specific order: the &lsquo;longitude&rsquo; coordinate needs to ascend first, followed by the &lsquo;colatitude&rsquo; coordinate in order to assign the correct data (isotherm depth) to the prescribed coordinates. The temperature is defined from the surface (273.15 K) to the isotherm depth (1673.15 K) as a linear gradient. Below the isotherm depth the temperature increases approximately adiabatically (0.0005 K per meter). This plugin should work for all geometry models, but is currently only tested for spherical models.

&lsquo;ascii data&rsquo;: Implementation of a model in which the initial temperature is derived from files containing data in ascii format. Note the required format of the input data: The first lines may contain any number of comments if they begin with &lsquo;#&rsquo;, but one of these lines needs to contain the number of grid points in each dimension as for example &lsquo;# POINTS: 3 3&rsquo;. The order of the data columns has to be &lsquo;x&rsquo;, &lsquo;y&rsquo;, &lsquo;Temperature [K]&rsquo; in a 2d model and  &lsquo;x&rsquo;, &lsquo;y&rsquo;, &lsquo;z&rsquo;, &lsquo;Temperature [K]&rsquo; in a 3d model, which means that there has to be a single column containing the temperature. Note that the data in the input files need to be sorted in a specific order: the first coordinate needs to ascend first, followed by the second and the third at last in order to assign the correct data to the prescribed coordinates. If you use a spherical model, then the assumed grid changes. &lsquo;x&rsquo; will be replaced by the radial distance of the point to the bottom of the model, &lsquo;y&rsquo; by the azimuth angle and &lsquo;z&rsquo; by the polar angle measured positive from the north pole. The grid will be assumed to be a latitude-longitude grid. Note that the order of spherical coordinates is &lsquo;r&rsquo;, &lsquo;phi&rsquo;, &lsquo;theta&rsquo; and not &lsquo;r&rsquo;, &lsquo;theta&rsquo;, &lsquo;phi&rsquo;, since this allows for dimension independent expressions.

&lsquo;ascii data layered&rsquo;: Implementation of a model in which the initial temperature is derived from files containing data in ascii format. Each file defines a surface on which temperature is defined. Between the surfaces, the temperatures can be chosen to be constant (with a value defined by the nearest shallower surface), or linearly interpolated between surfaces. Note the required format of the input ascii data file: The first lines may contain any number of comments if they begin with &lsquo;#&rsquo;, but one of these lines needs to contain the number of grid points in each dimension as for example &lsquo;# POINTS: 3 3&rsquo;. The order of the data columns has to be &lsquo;x&rsquo;, &lsquo;y&rsquo;, &lsquo;Temperature [K]&rsquo; in a 2d model and &lsquo;x&rsquo;, &lsquo;y&rsquo;, &lsquo;z&rsquo;, &lsquo;Temperature [K]&rsquo; in a 3d model; i.e. the last two columns always contain the position of the isotherm along the vertical direction, and the temperature at that point. The first column needs to ascend first, followed by the second in order to assign the correct data to the prescribed coordinates. If you use a spherical model, then the assumed grid changes. &lsquo;x&rsquo; will be replaced by the azimuth angle and &lsquo;y&rsquo; (if 3d) by the polar angle measured positive from the north pole. The last column will be the distance of the point from the origin (i.e. radial position). The grid in this case will be a latitude-longitude grid. Note that the order of spherical coordinates in 3d is &lsquo;phi&rsquo;, &lsquo;theta&rsquo;, &lsquo;r&rsquo;, &lsquo;T&rsquo;and not &lsquo;theta&rsquo;, &lsquo;phi&rsquo;, &lsquo;r&rsquo;, &lsquo;T&rsquo; as this is more consistent with other ASPECT plugins. Outside of the region defined by the grid, the plugin will use the value at the edge of the region.

&lsquo;ascii profile&rsquo;: Implementation of a model in which the initial temperature is read from a file that provides these values as a function of depth. Note the required format of the input data: The first lines may contain any number of comments if they begin with &lsquo;#&rsquo;, but one of these lines needs to contain the number of points in the temperature profile, for example &lsquo;# POINTS: 10&rsquo;. Following the comment lines, there has to be a single line containing the names of all data columns, separated by arbitrarily many spaces. Column names are not allowed to contain spaces. The file can contain unnecessary columns, but for this plugin it needs to at least provide columns named &lsquo;depth&rsquo; and&lsquo;temperature&rsquo;.Note that the data lines in the file need to be sorted in order of increasing depth from 0 to the maximal depth in the model domain. Points in the model that are outside of the provided depth range will be assigned the maximum or minimum depth values, respectively. Points do not need to be equidistant, but the computation of properties is optimized in speed if they are.

&lsquo;continental geotherm&rsquo;: This is a temperature initial condition that computes a continental geotherm based on the solution of the steady-state conductive equation $k\frac{d^2 T}{dy^2}+\rho H = 0$ as described in e.g. Turcotte and Schubert, Ch. 4.6, or Chapman (1986). As boundary conditions, we take the surface temperature and the temperature of the Lithosphere-Asthenosphere Boundary (LAB).
The geotherm is computed for a homogeneous lithosphere composed of an upper crust, lower crust and mantle layer. The crustal layers are assumed to have a constant radioactive heating, and all layers are assumed to have a constant thermal conductivity. Layer thicknesses, surface temperature and LAB temperature should be specified by the user. For consistency, the density, heat production and thermal conductivity of each layer are read from the visco plastic material model and the compositional heating model.
For any depths below the depth of the LAB, a unrealistically high temperature is returned, such that this plugin can be combined with another temperature plugin through the &rsquo;minimum&rsquo; operator.
Note that the current implementation only works for a 3-layer lithosphere, even though in principle the heat conduction equation can be solved for any number of layers. The naming of the compositional fields that represent the layers is also very specific, namely &lsquo;upper\_crust&rsquo;, &lsquo;lower\_crust&rsquo;, and &lsquo;lithospheric\_mantle&rsquo;.
Make sure the top and bottom temperatures of the lithosphere agree with temperatures set in for example the temperature boundary conditions.

&lsquo;function&rsquo;: Specify the initial temperature in terms of an explicit formula. The format of these functions follows the syntax understood by the muparser library, see {ref}`sec:run-aspect:parameters-overview:muparser-format`.

&lsquo;harmonic perturbation&rsquo;: An initial temperature field in which the temperature is perturbed following a harmonic function (spherical harmonic or sine depending on geometry and dimension) in lateral and radial direction from an otherwise constant temperature (incompressible model) or adiabatic reference profile (compressible model).

&lsquo;inclusion shape perturbation&rsquo;: An initial temperature field in which there is an inclusion in a constant-temperature box field. The size, shape, gradient, position, and temperature of the inclusion are defined by parameters.

&lsquo;lithosphere mask&rsquo;: Implementation of a model in which the initial temperature is set to a specified lithosphere temperature above the lithosphere-asthenosphere boundary (specified by an ascii file or maximum lithosphere depth value). Below this the initial temperature is set as NaN.  Note the required format of the input data file: The first lines may contain any number of comments if they begin with &rsquo;#&rsquo;, but one of these lines needs to contain the number of grid points in each dimension as for example &rsquo;# POINTS: 3 3&rsquo;. For a spherical model, the order of the data columns has to be &rsquo;phi&rsquo;, &rsquo;theta&rsquo;, &rsquo;depth (m)&rsquo;, where phi is the azimuth angle and theta is the polar angle measured positive from the north pole. This plug-in can be combined with another using the &rsquo;replace if valid&rsquo; operator.

&lsquo;mandelbox&rsquo;: Fractal-shaped temperature field.

&lsquo;patch on S40RTS&rsquo;: Implementation of a model in which the initial temperature is derived from a file containing shear wave velocity perturbations in ascii format (e.g. a high resolution upper mantle tomography) combined with S40RTS. Note the required format of the input ascii input data: The first lines may contain any number of comments if they begin with &rsquo;#&rsquo;, but one of these lines needs to contain the number of grid points in each dimension as for example &rsquo;# POINTS: 3 3 3&rsquo;. The order of the data columns has to be  &lsquo;x&rsquo;, &lsquo;y&rsquo;, &lsquo;z&rsquo;, &rsquo;Vs Perturbation&rsquo; in a 3d model, which means that there has to be a single column containing the temperature. Note that the data in the input files need to be sorted in a specific order: the first coordinate needs to ascend first, followed by the second and the third at last in order to assign the correct data to the prescribed coordinates. In the spherical model data will be handled as Cartesian, however, &lsquo;x&rsquo; will be replaced by the radial distance of the point to the bottom of the model, &lsquo;y&rsquo; by the azimuth angle and &lsquo;z&rsquo; by the polar angle measured positive from the north pole. The grid will be assumed to be a latitude-longitude grid. Note that the order of spherical coordinates is &lsquo;r&rsquo;, &lsquo;phi&rsquo;, &lsquo;theta&rsquo; and not &lsquo;r&rsquo;, &lsquo;theta&rsquo;, &lsquo;phi&rsquo;, since this allows for dimension independent expressions. See S40RTS documentation for details on input parameters in the S40RTS perturbation subsection. The boundary between the two tomography models is smoothed using a depth weighted combination of Vs values within the region of smoothing.

&lsquo;perturbed box&rsquo;: An initial temperature field in which the temperature is perturbed slightly from an otherwise constant value equal to one. The perturbation is chosen in such a way that the initial temperature is constant to one along the entire boundary.

&lsquo;polar box&rsquo;: An initial temperature field in which the temperature is perturbed slightly from an otherwise constant value equal to one. The perturbation is such that there are two poles on opposing corners of the box.

&lsquo;prescribed temperature&rsquo;: This model fixes the initial temperature to the prescribed temperature outputs computed by the material model. This only works if the material model implements prescribed temperature outputs.

&lsquo;random Gaussian perturbation&rsquo;: An initial temperature field in which the temperature is perturbed from a temperature of zero following a given number of Gaussian perturbations placed randomly throughout the model domain. The number, width, and maximum magnitude of the perturbations can be chosen as model parameters. This plugin is meant to be used in combination with another initial temperature model that determines the background temperature (such as the &rsquo;function&rsquo; or the &rsquo;adiabatic&rsquo; plugin) using the &rsquo;add&rsquo; operator to combine them.

&lsquo;spherical gaussian perturbation&rsquo;: An initial temperature field in which the temperature is perturbed by a single Gaussian added to an otherwise spherically symmetric state. Additional parameters are read from the parameter file in subsection &rsquo;Spherical gaussian perturbation&rsquo;.

&lsquo;spherical hexagonal perturbation&rsquo;: An initial temperature field in which the temperature is perturbed following an $N$-fold pattern in a specified direction from an otherwise spherically symmetric state. The class&rsquo;s name comes from previous versions when the only option was $N=6$.

&lsquo;world builder&rsquo;: Specify the initial temperature through the World Builder. More information on the World Builder can be found at \url{https://geodynamicworldbuilder.github.io}. Make sure to specify the location of the World Builder file in the parameter &rsquo;World builder file&rsquo;.

(parameters:Initial_20temperature_20model/List_20of_20model_20operators)=
### __Parameter name:__ List of model operators
**Default value:** add

**Pattern:** [MultipleSelection add|subtract|minimum|maximum|replace if valid ]

**Documentation:** A comma-separated list of operators that will be used to append the listed temperature models onto the previous models. If only one operator is given, the same operator is applied to all models.

(parameters:Initial_20temperature_20model/Model_20name)=
### __Parameter name:__ Model name
**Default value:** unspecified

**Pattern:** [Selection S40RTS perturbation|SAVANI perturbation|adiabatic|adiabatic boundary|ascii data|ascii data layered|ascii profile|continental geotherm|function|harmonic perturbation|inclusion shape perturbation|lithosphere mask|mandelbox|patch on S40RTS|perturbed box|polar box|prescribed temperature|random Gaussian perturbation|spherical gaussian perturbation|spherical hexagonal perturbation|world builder|unspecified ]

**Documentation:** Select one of the following models:

&lsquo;S40RTS perturbation&rsquo;: An initial temperature field in which the temperature is perturbed following the S20RTS or S40RTS shear wave velocity model by Ritsema and others, which can be downloaded here \url{http://www.earth.lsa.umich.edu/~jritsema/research.html}. Information on the vs model can be found in Ritsema, J., Deuss, A., van Heijst, H.J. \& Woodhouse, J.H., 2011. S40RTS: a degree-40 shear-velocity model for the mantle from new Rayleigh wave dispersion, teleseismic traveltime and normal-mode splitting function measurements, Geophys. J. Int. 184, 1223-1236. The scaling between the shear wave perturbation and the density perturbation can be constant and set by the user with the &rsquo;Vs to density scaling&rsquo; parameter or depth-dependent and read in from a file. To convert density the user can specify the &rsquo;Thermal expansion coefficient in initial temperature scaling&rsquo; parameter. The scaling is as follows: $\delta \ln \rho (r,\theta,\phi) = \xi \cdot \delta \ln v_s(r,\theta, \phi)$ and $\delta T(r,\theta,\phi) = - \frac{1}{\alpha} \delta \ln \rho(r,\theta,\phi)$. $\xi$ is the &lsquo;vs to density scaling&rsquo; parameter and $\alpha$ is the &rsquo;Thermal expansion coefficient in initial temperature scaling&rsquo; parameter. The temperature perturbation is added to an otherwise constant temperature (incompressible model) or adiabatic reference profile (compressible model). If a depth is specified in &rsquo;Remove temperature heterogeneity down to specified depth&rsquo;, there is no temperature perturbation prescribed down to that depth.
Note the required file format if the vs to density scaling is read in from a file: The first lines may contain any number of comments if they begin with &rsquo;#&rsquo;, but one of these lines needs to contain the number of points in the reference state as for example &rsquo;# POINTS: 3&rsquo;. Following the comment lines there has to be a single line containing the names of all data columns, separated by arbitrarily many spaces. Column names are not allowed to contain spaces. The file can contain unnecessary columns, but for this plugin it needs to at least provide the columns named &lsquo;depth&rsquo; and &lsquo;vs\_to\_density&rsquo;. Note that the data lines in the file need to be sorted in order of increasing depth from 0 to the maximal depth in the model domain. Points in the model that are outside of the provided depth range will be assigned the maximum or minimum depth values, respectively. Points do not need to be equidistant, but the computation of properties is optimized in speed if they are.
If the plugin is used in 2d it will use an equatorial slice of the seismic tomography model.

&lsquo;SAVANI perturbation&rsquo;: An initial temperature field in which the temperature is perturbed following the SAVANI shear wave velocity model by Auer and others, which can be downloaded here \url{http://n.ethz.ch/~auerl/savani.tar.bz2}. Information on the vs model can be found in Auer, L., Boschi, L., Becker, T.W., Nissen-Meyer, T. \& Giardini, D., 2014. Savani: A variable resolution whole-mantle model of anisotropic shear velocity variations based on multiple data sets. Journal of Geophysical Research: Solid Earth 119.4 (2014): 3006-3034. The scaling between the shear wave perturbation and the density perturbation can be constant and set by the user with the &rsquo;Vs to density scaling&rsquo; parameter or depth-dependent and read in from a file. To convert density the user can specify the &rsquo;Thermal expansion coefficient in initial temperature scaling&rsquo; parameter. The scaling is as follows: $\delta \ln \rho (r,\theta,\phi) = \xi \cdot \delta \ln v_s(r,\theta, \phi)$ and $\delta T(r,\theta,\phi) = - \frac{1}{\alpha} \delta \ln \rho(r,\theta,\phi)$. $\xi$ is the &lsquo;vs to density scaling&rsquo; parameter and $\alpha$ is the &rsquo;Thermal expansion coefficient in initial temperature scaling&rsquo; parameter. The temperature perturbation is added to an otherwise constant temperature (incompressible model) or adiabatic reference profile (compressible model).If a depth is specified in &rsquo;Remove temperature heterogeneity down to specified depth&rsquo;, there is no temperature perturbation prescribed down to that depth.
Note the required file format if the vs to density scaling is read in from a file: The first lines may contain any number of comments if they begin with &rsquo;#&rsquo;, but one of these lines needs to contain the number of points in the reference state as for example &rsquo;# POINTS: 3&rsquo;. Following the comment lines there has to be a single line containing the names of all data columns, separated by arbitrarily many spaces. Column names are not allowed to contain spaces. The file can contain unnecessary columns, but for this plugin it needs to at least provide the columns named &lsquo;depth&rsquo; and &lsquo;vs\_to\_density&rsquo;. Note that the data lines in the file need to be sorted in order of increasing depth from 0 to the maximal depth in the model domain. Points in the model that are outside of the provided depth range will be assigned the maximum or minimum depth values, respectively. Points do not need to be equidistant, but the computation of properties is optimized in speed if they are.

&lsquo;adiabatic&rsquo;: Temperature is prescribed as an adiabatic profile with upper and lower thermal boundary layers, whose ages are given as input parameters. Note that this plugin uses the &rsquo;Adiabatic conditions model&rsquo; to compute the adiabat. Thus, the results depend on variables defined outside of this specific subsection; e.g. the globally defined &rsquo;Adiabatic surface temperature&rsquo;, and the variables defined in the &rsquo;Material model&rsquo; section including densities, heat capacities and thermal expansivities.

&lsquo;adiabatic boundary&rsquo;: An initial temperature condition that allows for discretizing the model domain into two layers separated by a user-defined isothermal boundary. The user includes an input ascii data file that is formatted as 3 columns of &lsquo;longitude(radians)&rsquo;, &lsquo;colatitude(radians)&rsquo;, and &lsquo;isotherm depth(meters)&rsquo;, where &lsquo;isotherm depth&rsquo; represents the depth of an initial temperature of 1673.15 K (by default). The first lines in the data file may contain any number of comments if they begin with &lsquo;#&rsquo;, but one of these lines needs to contain the number of grid points in each dimension as for example &lsquo;# POINTS: 69 121&rsquo;. Note that the coordinates need to be sorted in a specific order: the &lsquo;longitude&rsquo; coordinate needs to ascend first, followed by the &lsquo;colatitude&rsquo; coordinate in order to assign the correct data (isotherm depth) to the prescribed coordinates. The temperature is defined from the surface (273.15 K) to the isotherm depth (1673.15 K) as a linear gradient. Below the isotherm depth the temperature increases approximately adiabatically (0.0005 K per meter). This plugin should work for all geometry models, but is currently only tested for spherical models.

&lsquo;ascii data&rsquo;: Implementation of a model in which the initial temperature is derived from files containing data in ascii format. Note the required format of the input data: The first lines may contain any number of comments if they begin with &lsquo;#&rsquo;, but one of these lines needs to contain the number of grid points in each dimension as for example &lsquo;# POINTS: 3 3&rsquo;. The order of the data columns has to be &lsquo;x&rsquo;, &lsquo;y&rsquo;, &lsquo;Temperature [K]&rsquo; in a 2d model and  &lsquo;x&rsquo;, &lsquo;y&rsquo;, &lsquo;z&rsquo;, &lsquo;Temperature [K]&rsquo; in a 3d model, which means that there has to be a single column containing the temperature. Note that the data in the input files need to be sorted in a specific order: the first coordinate needs to ascend first, followed by the second and the third at last in order to assign the correct data to the prescribed coordinates. If you use a spherical model, then the assumed grid changes. &lsquo;x&rsquo; will be replaced by the radial distance of the point to the bottom of the model, &lsquo;y&rsquo; by the azimuth angle and &lsquo;z&rsquo; by the polar angle measured positive from the north pole. The grid will be assumed to be a latitude-longitude grid. Note that the order of spherical coordinates is &lsquo;r&rsquo;, &lsquo;phi&rsquo;, &lsquo;theta&rsquo; and not &lsquo;r&rsquo;, &lsquo;theta&rsquo;, &lsquo;phi&rsquo;, since this allows for dimension independent expressions.

&lsquo;ascii data layered&rsquo;: Implementation of a model in which the initial temperature is derived from files containing data in ascii format. Each file defines a surface on which temperature is defined. Between the surfaces, the temperatures can be chosen to be constant (with a value defined by the nearest shallower surface), or linearly interpolated between surfaces. Note the required format of the input ascii data file: The first lines may contain any number of comments if they begin with &lsquo;#&rsquo;, but one of these lines needs to contain the number of grid points in each dimension as for example &lsquo;# POINTS: 3 3&rsquo;. The order of the data columns has to be &lsquo;x&rsquo;, &lsquo;y&rsquo;, &lsquo;Temperature [K]&rsquo; in a 2d model and &lsquo;x&rsquo;, &lsquo;y&rsquo;, &lsquo;z&rsquo;, &lsquo;Temperature [K]&rsquo; in a 3d model; i.e. the last two columns always contain the position of the isotherm along the vertical direction, and the temperature at that point. The first column needs to ascend first, followed by the second in order to assign the correct data to the prescribed coordinates. If you use a spherical model, then the assumed grid changes. &lsquo;x&rsquo; will be replaced by the azimuth angle and &lsquo;y&rsquo; (if 3d) by the polar angle measured positive from the north pole. The last column will be the distance of the point from the origin (i.e. radial position). The grid in this case will be a latitude-longitude grid. Note that the order of spherical coordinates in 3d is &lsquo;phi&rsquo;, &lsquo;theta&rsquo;, &lsquo;r&rsquo;, &lsquo;T&rsquo;and not &lsquo;theta&rsquo;, &lsquo;phi&rsquo;, &lsquo;r&rsquo;, &lsquo;T&rsquo; as this is more consistent with other ASPECT plugins. Outside of the region defined by the grid, the plugin will use the value at the edge of the region.

&lsquo;ascii profile&rsquo;: Implementation of a model in which the initial temperature is read from a file that provides these values as a function of depth. Note the required format of the input data: The first lines may contain any number of comments if they begin with &lsquo;#&rsquo;, but one of these lines needs to contain the number of points in the temperature profile, for example &lsquo;# POINTS: 10&rsquo;. Following the comment lines, there has to be a single line containing the names of all data columns, separated by arbitrarily many spaces. Column names are not allowed to contain spaces. The file can contain unnecessary columns, but for this plugin it needs to at least provide columns named &lsquo;depth&rsquo; and&lsquo;temperature&rsquo;.Note that the data lines in the file need to be sorted in order of increasing depth from 0 to the maximal depth in the model domain. Points in the model that are outside of the provided depth range will be assigned the maximum or minimum depth values, respectively. Points do not need to be equidistant, but the computation of properties is optimized in speed if they are.

&lsquo;continental geotherm&rsquo;: This is a temperature initial condition that computes a continental geotherm based on the solution of the steady-state conductive equation $k\frac{d^2 T}{dy^2}+\rho H = 0$ as described in e.g. Turcotte and Schubert, Ch. 4.6, or Chapman (1986). As boundary conditions, we take the surface temperature and the temperature of the Lithosphere-Asthenosphere Boundary (LAB).
The geotherm is computed for a homogeneous lithosphere composed of an upper crust, lower crust and mantle layer. The crustal layers are assumed to have a constant radioactive heating, and all layers are assumed to have a constant thermal conductivity. Layer thicknesses, surface temperature and LAB temperature should be specified by the user. For consistency, the density, heat production and thermal conductivity of each layer are read from the visco plastic material model and the compositional heating model.
For any depths below the depth of the LAB, a unrealistically high temperature is returned, such that this plugin can be combined with another temperature plugin through the &rsquo;minimum&rsquo; operator.
Note that the current implementation only works for a 3-layer lithosphere, even though in principle the heat conduction equation can be solved for any number of layers. The naming of the compositional fields that represent the layers is also very specific, namely &lsquo;upper\_crust&rsquo;, &lsquo;lower\_crust&rsquo;, and &lsquo;lithospheric\_mantle&rsquo;.
Make sure the top and bottom temperatures of the lithosphere agree with temperatures set in for example the temperature boundary conditions.

&lsquo;function&rsquo;: Specify the initial temperature in terms of an explicit formula. The format of these functions follows the syntax understood by the muparser library, see {ref}`sec:run-aspect:parameters-overview:muparser-format`.

&lsquo;harmonic perturbation&rsquo;: An initial temperature field in which the temperature is perturbed following a harmonic function (spherical harmonic or sine depending on geometry and dimension) in lateral and radial direction from an otherwise constant temperature (incompressible model) or adiabatic reference profile (compressible model).

&lsquo;inclusion shape perturbation&rsquo;: An initial temperature field in which there is an inclusion in a constant-temperature box field. The size, shape, gradient, position, and temperature of the inclusion are defined by parameters.

&lsquo;lithosphere mask&rsquo;: Implementation of a model in which the initial temperature is set to a specified lithosphere temperature above the lithosphere-asthenosphere boundary (specified by an ascii file or maximum lithosphere depth value). Below this the initial temperature is set as NaN.  Note the required format of the input data file: The first lines may contain any number of comments if they begin with &rsquo;#&rsquo;, but one of these lines needs to contain the number of grid points in each dimension as for example &rsquo;# POINTS: 3 3&rsquo;. For a spherical model, the order of the data columns has to be &rsquo;phi&rsquo;, &rsquo;theta&rsquo;, &rsquo;depth (m)&rsquo;, where phi is the azimuth angle and theta is the polar angle measured positive from the north pole. This plug-in can be combined with another using the &rsquo;replace if valid&rsquo; operator.

&lsquo;mandelbox&rsquo;: Fractal-shaped temperature field.

&lsquo;patch on S40RTS&rsquo;: Implementation of a model in which the initial temperature is derived from a file containing shear wave velocity perturbations in ascii format (e.g. a high resolution upper mantle tomography) combined with S40RTS. Note the required format of the input ascii input data: The first lines may contain any number of comments if they begin with &rsquo;#&rsquo;, but one of these lines needs to contain the number of grid points in each dimension as for example &rsquo;# POINTS: 3 3 3&rsquo;. The order of the data columns has to be  &lsquo;x&rsquo;, &lsquo;y&rsquo;, &lsquo;z&rsquo;, &rsquo;Vs Perturbation&rsquo; in a 3d model, which means that there has to be a single column containing the temperature. Note that the data in the input files need to be sorted in a specific order: the first coordinate needs to ascend first, followed by the second and the third at last in order to assign the correct data to the prescribed coordinates. In the spherical model data will be handled as Cartesian, however, &lsquo;x&rsquo; will be replaced by the radial distance of the point to the bottom of the model, &lsquo;y&rsquo; by the azimuth angle and &lsquo;z&rsquo; by the polar angle measured positive from the north pole. The grid will be assumed to be a latitude-longitude grid. Note that the order of spherical coordinates is &lsquo;r&rsquo;, &lsquo;phi&rsquo;, &lsquo;theta&rsquo; and not &lsquo;r&rsquo;, &lsquo;theta&rsquo;, &lsquo;phi&rsquo;, since this allows for dimension independent expressions. See S40RTS documentation for details on input parameters in the S40RTS perturbation subsection. The boundary between the two tomography models is smoothed using a depth weighted combination of Vs values within the region of smoothing.

&lsquo;perturbed box&rsquo;: An initial temperature field in which the temperature is perturbed slightly from an otherwise constant value equal to one. The perturbation is chosen in such a way that the initial temperature is constant to one along the entire boundary.

&lsquo;polar box&rsquo;: An initial temperature field in which the temperature is perturbed slightly from an otherwise constant value equal to one. The perturbation is such that there are two poles on opposing corners of the box.

&lsquo;prescribed temperature&rsquo;: This model fixes the initial temperature to the prescribed temperature outputs computed by the material model. This only works if the material model implements prescribed temperature outputs.

&lsquo;random Gaussian perturbation&rsquo;: An initial temperature field in which the temperature is perturbed from a temperature of zero following a given number of Gaussian perturbations placed randomly throughout the model domain. The number, width, and maximum magnitude of the perturbations can be chosen as model parameters. This plugin is meant to be used in combination with another initial temperature model that determines the background temperature (such as the &rsquo;function&rsquo; or the &rsquo;adiabatic&rsquo; plugin) using the &rsquo;add&rsquo; operator to combine them.

&lsquo;spherical gaussian perturbation&rsquo;: An initial temperature field in which the temperature is perturbed by a single Gaussian added to an otherwise spherically symmetric state. Additional parameters are read from the parameter file in subsection &rsquo;Spherical gaussian perturbation&rsquo;.

&lsquo;spherical hexagonal perturbation&rsquo;: An initial temperature field in which the temperature is perturbed following an $N$-fold pattern in a specified direction from an otherwise spherically symmetric state. The class&rsquo;s name comes from previous versions when the only option was $N=6$.

&lsquo;world builder&rsquo;: Specify the initial temperature through the World Builder. More information on the World Builder can be found at \url{https://geodynamicworldbuilder.github.io}. Make sure to specify the location of the World Builder file in the parameter &rsquo;World builder file&rsquo;.

**Warning**: This parameter provides an old and deprecated way of specifying initial temperature models and shouldn&rsquo;t be used. Please use &rsquo;List of model names&rsquo; instead.

(parameters:Initial_20temperature_20model/Adiabatic)=
## **Subsection:** Initial temperature model / Adiabatic
(parameters:Initial_20temperature_20model/Adiabatic/Age_20bottom_20boundary_20layer)=
### __Parameter name:__ Age bottom boundary layer
**Default value:** 0.

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The age of the lower thermal boundary layer, used for the calculation of the half-space cooling model temperature. Units: years if the &rsquo;Use years in output instead of seconds&rsquo; parameter is set; seconds otherwise.

(parameters:Initial_20temperature_20model/Adiabatic/Age_20top_20boundary_20layer)=
### __Parameter name:__ Age top boundary layer
**Default value:** 0.

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The age of the upper thermal boundary layer, used for the calculation of the half-space cooling model temperature. Units: years if the &rsquo;Use years in output instead of seconds&rsquo; parameter is set; seconds otherwise.

(parameters:Initial_20temperature_20model/Adiabatic/Amplitude)=
### __Parameter name:__ Amplitude
**Default value:** 0.

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The amplitude (in K) of the initial spherical temperature perturbation at the bottom of the model domain. This perturbation will be added to the adiabatic temperature profile, but not to the bottom thermal boundary layer. Instead, the maximum of the perturbation and the bottom boundary layer temperature will be used.

(parameters:Initial_20temperature_20model/Adiabatic/Cooling_20model)=
### __Parameter name:__ Cooling model
**Default value:** half-space cooling

**Pattern:** [Selection half-space cooling|plate cooling ]

**Documentation:** Whether to use the half space cooling model or the plate cooling model

(parameters:Initial_20temperature_20model/Adiabatic/Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/initial-temperature/adiabatic/

**Pattern:** [DirectoryName]

**Documentation:** The name of a directory that contains the model data. This path may either be absolute (if starting with a &lsquo;/&rsquo;) or relative to the current directory. The path may also include the special text &lsquo;$ASPECT_SOURCE_DIR&rsquo; which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the &lsquo;data/&rsquo; subdirectory of ASPECT.

(parameters:Initial_20temperature_20model/Adiabatic/Data_20file_20name)=
### __Parameter name:__ Data file name
**Default value:** adiabatic.txt

**Pattern:** [Anything]

**Documentation:** The file name of the model data.

(parameters:Initial_20temperature_20model/Adiabatic/Lithosphere_20thickness)=
### __Parameter name:__ Lithosphere thickness
**Default value:** 125e3

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Thickness of the lithosphere for plate cooling model. \si{\m}

(parameters:Initial_20temperature_20model/Adiabatic/Position)=
### __Parameter name:__ Position
**Default value:** center

**Pattern:** [Selection center ]

**Documentation:** Where the initial temperature perturbation should be placed. If &lsquo;center&rsquo; is given, then the perturbation will be centered along a &lsquo;midpoint&rsquo; of some sort of the bottom boundary. For example, in the case of a box geometry, this is the center of the bottom face; in the case of a spherical shell geometry, it is along the inner surface halfway between the bounding radial lines.

(parameters:Initial_20temperature_20model/Adiabatic/Radius)=
### __Parameter name:__ Radius
**Default value:** 0.

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The Radius (in m) of the initial spherical temperature perturbation at the bottom of the model domain.

(parameters:Initial_20temperature_20model/Adiabatic/Scale_20factor)=
### __Parameter name:__ Scale factor
**Default value:** 1.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Scalar factor, which is applied to the model data. You might want to use this to scale the input to a reference model. Another way to use this factor is to convert units of the input files. For instance, if you provide velocities in cm/yr set this factor to 0.01.

(parameters:Initial_20temperature_20model/Adiabatic/Subadiabaticity)=
### __Parameter name:__ Subadiabaticity
**Default value:** 0.

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** If this value is larger than 0, the initial temperature profile will not be adiabatic, but subadiabatic. This value gives the maximal deviation from adiabaticity. Set to 0 for an adiabatic temperature profile. Units: \si{\kelvin}.

The function object in the Function subsection represents the compositional fields that will be used as a reference profile for calculating the thermal diffusivity. This function is one-dimensional and depends only on depth. The format of this functions follows the syntax understood by the muparser library, see {ref}`sec:run-aspect:parameters-overview:muparser-format`.

(parameters:Initial_20temperature_20model/Adiabatic/Top_20boundary_20layer_20age_20model)=
### __Parameter name:__ Top boundary layer age model
**Default value:** constant

**Pattern:** [Selection constant|function|ascii data ]

**Documentation:** How to define the age of the top thermal boundary layer. Options are: &rsquo;constant&rsquo; for a constant age specified by the parameter &rsquo;Age top boundary layer&rsquo;; &rsquo;function&rsquo; for an analytical function describing the age as specified in the subsection &rsquo;Age function&rsquo;; and &rsquo;ascii data&rsquo; to use an &rsquo;ascii data&rsquo; file specified by the parameter &rsquo;Data file name&rsquo;.

(parameters:Initial_20temperature_20model/Adiabatic/Age_20function)=
## **Subsection:** Initial temperature model / Adiabatic / Age function
(parameters:Initial_20temperature_20model/Adiabatic/Age_20function/Coordinate_20system)=
### __Parameter name:__ Coordinate system
**Default value:** cartesian

**Pattern:** [Selection cartesian|spherical ]

**Documentation:** A selection that determines the assumed coordinate system for the function variables. Allowed values are &lsquo;cartesian&rsquo;, &lsquo;spherical&rsquo;, and &lsquo;depth&rsquo;. &lsquo;spherical&rsquo; coordinates are interpreted as r,phi or r,phi,theta in 2d/3d respectively with theta being the polar angle. &lsquo;depth&rsquo; will create a function, in which only the first parameter is non-zero, which is interpreted to be the depth of the point.

(parameters:Initial_20temperature_20model/Adiabatic/Age_20function/Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**

**Pattern:** [Anything]

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form &lsquo;var1=value1, var2=value2, ...&rsquo;.

A typical example would be to set this runtime parameter to &lsquo;pi=3.1415926536&rsquo; and then use &lsquo;pi&rsquo; in the expression of the actual formula. (That said, for convenience this class actually defines both &lsquo;pi&rsquo; and &lsquo;Pi&rsquo; by default, but you get the idea.)

(parameters:Initial_20temperature_20model/Adiabatic/Age_20function/Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0

**Pattern:** [Anything]

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as &lsquo;sin&rsquo; or &lsquo;cos&rsquo;. In addition, it may contain expressions like &lsquo;if(x>0, 1, -1)&rsquo; where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon.

(parameters:Initial_20temperature_20model/Adiabatic/Age_20function/Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,y,t

**Pattern:** [Anything]

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are &lsquo;x&rsquo; (in 1d), &lsquo;x,y&rsquo; (in 2d) or &lsquo;x,y,z&rsquo; (in 3d) for spatial coordinates and &lsquo;t&rsquo; for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to &lsquo;r,phi,theta,t&rsquo; and then use these variable names in your function expression.

(parameters:Initial_20temperature_20model/Adiabatic/Function)=
## **Subsection:** Initial temperature model / Adiabatic / Function
(parameters:Initial_20temperature_20model/Adiabatic/Function/Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**

**Pattern:** [Anything]

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form &lsquo;var1=value1, var2=value2, ...&rsquo;.

A typical example would be to set this runtime parameter to &lsquo;pi=3.1415926536&rsquo; and then use &lsquo;pi&rsquo; in the expression of the actual formula. (That said, for convenience this class actually defines both &lsquo;pi&rsquo; and &lsquo;Pi&rsquo; by default, but you get the idea.)

(parameters:Initial_20temperature_20model/Adiabatic/Function/Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0

**Pattern:** [Anything]

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as &lsquo;sin&rsquo; or &lsquo;cos&rsquo;. In addition, it may contain expressions like &lsquo;if(x>0, 1, -1)&rsquo; where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon.

(parameters:Initial_20temperature_20model/Adiabatic/Function/Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,t

**Pattern:** [Anything]

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are &lsquo;x&rsquo; (in 1d), &lsquo;x,y&rsquo; (in 2d) or &lsquo;x,y,z&rsquo; (in 3d) for spatial coordinates and &lsquo;t&rsquo; for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to &lsquo;r,phi,theta,t&rsquo; and then use these variable names in your function expression.

(parameters:Initial_20temperature_20model/Adiabatic_20boundary)=
## **Subsection:** Initial temperature model / Adiabatic boundary
(parameters:Initial_20temperature_20model/Adiabatic_20boundary/Adiabatic_20temperature_20gradient)=
### __Parameter name:__ Adiabatic temperature gradient
**Default value:** 0.0005

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The value of the adiabatic temperature gradient. Units: \si{\kelvin\per\meter}.

(parameters:Initial_20temperature_20model/Adiabatic_20boundary/Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/initial-temperature/adiabatic-boundary/

**Pattern:** [DirectoryName]

**Documentation:** The name of a directory that contains the model data. This path may either be absolute (if starting with a &lsquo;/&rsquo;) or relative to the current directory. The path may also include the special text &lsquo;$ASPECT_SOURCE_DIR&rsquo; which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the &lsquo;data/&rsquo; subdirectory of ASPECT.

(parameters:Initial_20temperature_20model/Adiabatic_20boundary/Data_20file_20name)=
### __Parameter name:__ Data file name
**Default value:** adiabatic_boundary.txt

**Pattern:** [Anything]

**Documentation:** The file name of the model data.

(parameters:Initial_20temperature_20model/Adiabatic_20boundary/Isotherm_20temperature)=
### __Parameter name:__ Isotherm temperature
**Default value:** 1673.15

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The value of the isothermal boundary temperature. Units: \si{\kelvin}.

(parameters:Initial_20temperature_20model/Adiabatic_20boundary/Scale_20factor)=
### __Parameter name:__ Scale factor
**Default value:** 1.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Scalar factor, which is applied to the model data. You might want to use this to scale the input to a reference model. Another way to use this factor is to convert units of the input files. For instance, if you provide velocities in cm/yr set this factor to 0.01.

(parameters:Initial_20temperature_20model/Adiabatic_20boundary/Surface_20temperature)=
### __Parameter name:__ Surface temperature
**Default value:** 273.15

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The value of the surface temperature. Units: \si{\kelvin}.

(parameters:Initial_20temperature_20model/Ascii_20data_20model)=
## **Subsection:** Initial temperature model / Ascii data model
(parameters:Initial_20temperature_20model/Ascii_20data_20model/Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/initial-temperature/ascii-data/test/

**Pattern:** [DirectoryName]

**Documentation:** The name of a directory that contains the model data. This path may either be absolute (if starting with a &lsquo;/&rsquo;) or relative to the current directory. The path may also include the special text &lsquo;$ASPECT_SOURCE_DIR&rsquo; which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the &lsquo;data/&rsquo; subdirectory of ASPECT.

(parameters:Initial_20temperature_20model/Ascii_20data_20model/Data_20file_20name)=
### __Parameter name:__ Data file name
**Default value:** initial_isotherm_500K_box_3d.txt

**Pattern:** [Anything]

**Documentation:** The file name of the model data.

(parameters:Initial_20temperature_20model/Ascii_20data_20model/Data_20file_20names)=
### __Parameter name:__ Data file names
**Default value:** initial_isotherm_500K_box_3d.txt

**Pattern:** [List of <[Anything]> of length 0...4294967295 (inclusive)]

**Documentation:** The file names of the model data (comma separated).

(parameters:Initial_20temperature_20model/Ascii_20data_20model/First_20point_20on_20slice)=
### __Parameter name:__ First point on slice
**Default value:** 0.0,1.0,0.0

**Pattern:** [Anything]

**Documentation:** Point that determines the plane in which the 2d slice lies in. This variable is only used if &rsquo;Slice dataset in 2d plane&rsquo; is true. The slice will go through this point, the point defined by the parameter &rsquo;Second point on slice&rsquo;, and the center of the model domain. After the rotation, this first point will lie along the (0,1,0) axis of the coordinate system. The coordinates of the point have to be given in Cartesian coordinates.

(parameters:Initial_20temperature_20model/Ascii_20data_20model/Interpolation_20scheme)=
### __Parameter name:__ Interpolation scheme
**Default value:** linear

**Pattern:** [Selection piecewise constant|linear ]

**Documentation:** Method to interpolate between layer boundaries. Select from piecewise constant or linear. Piecewise constant takes the value from the nearest layer boundary above the data point. The linear option interpolates linearly between layer boundaries. Above and below the domain given by the layer boundaries, the values aregiven by the top and bottom layer boundary.

(parameters:Initial_20temperature_20model/Ascii_20data_20model/Scale_20factor)=
### __Parameter name:__ Scale factor
**Default value:** 1.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Scalar factor, which is applied to the model data. You might want to use this to scale the input to a reference model. Another way to use this factor is to convert units of the input files. For instance, if you provide velocities in cm/yr set this factor to 0.01.

(parameters:Initial_20temperature_20model/Ascii_20data_20model/Second_20point_20on_20slice)=
### __Parameter name:__ Second point on slice
**Default value:** 1.0,0.0,0.0

**Pattern:** [Anything]

**Documentation:** Second point that determines the plane in which the 2d slice lies in. This variable is only used if &rsquo;Slice dataset in 2d plane&rsquo; is true. The slice will go through this point, the point defined by the parameter &rsquo;First point on slice&rsquo;, and the center of the model domain. The coordinates of the point have to be given in Cartesian coordinates.

(parameters:Initial_20temperature_20model/Ascii_20data_20model/Slice_20dataset_20in_202D_20plane)=
### __Parameter name:__ Slice dataset in 2D plane
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Whether to use a 2d data slice of a 3d data file or the entire data file. Slicing a 3d dataset is only supported for 2d models.

(parameters:Initial_20temperature_20model/Ascii_20profile)=
## **Subsection:** Initial temperature model / Ascii profile
(parameters:Initial_20temperature_20model/Ascii_20profile/Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/initial-temperature/ascii-profile/tests/

**Pattern:** [DirectoryName]

**Documentation:** The name of a directory that contains the model data. This path may either be absolute (if starting with a &lsquo;/&rsquo;) or relative to the current directory. The path may also include the special text &lsquo;$ASPECT_SOURCE_DIR&rsquo; which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the &lsquo;data/&rsquo; subdirectory of ASPECT.

(parameters:Initial_20temperature_20model/Ascii_20profile/Data_20file_20name)=
### __Parameter name:__ Data file name
**Default value:** simple_test.txt

**Pattern:** [Anything]

**Documentation:** The file name of the model data.

(parameters:Initial_20temperature_20model/Ascii_20profile/Scale_20factor)=
### __Parameter name:__ Scale factor
**Default value:** 1.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Scalar factor, which is applied to the model data. You might want to use this to scale the input to a reference model. Another way to use this factor is to convert units of the input files. For instance, if you provide velocities in cm/yr set this factor to 0.01.

(parameters:Initial_20temperature_20model/Continental_20geotherm)=
## **Subsection:** Initial temperature model / Continental geotherm
(parameters:Initial_20temperature_20model/Continental_20geotherm/Layer_20thicknesses)=
### __Parameter name:__ Layer thicknesses
**Default value:** 30000.

**Pattern:** [List of <[Double 0...MAX_DOUBLE (inclusive)]> of length 0...4294967295 (inclusive)]

**Documentation:** List of the 3 thicknesses of the lithospheric layers &rsquo;upper\_crust&rsquo;, &rsquo;lower\_crust&rsquo; and &rsquo;mantle\_lithosphere&rsquo;. If only one thickness is given, then the same thickness is used for all layers. Units: \si{meter}.

(parameters:Initial_20temperature_20model/Continental_20geotherm/Lithosphere_2dAsthenosphere_20boundary_20isotherm)=
### __Parameter name:__ Lithosphere-Asthenosphere boundary isotherm
**Default value:** 1673.15

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The value of the isotherm that is assumed at the Lithosphere-Asthenosphere boundary. Units: \si{\kelvin}.

(parameters:Initial_20temperature_20model/Continental_20geotherm/Surface_20temperature)=
### __Parameter name:__ Surface temperature
**Default value:** 273.15

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The value of the surface temperature. Units: \si{\kelvin}.

(parameters:Initial_20temperature_20model/Function)=
## **Subsection:** Initial temperature model / Function
(parameters:Initial_20temperature_20model/Function/Coordinate_20system)=
### __Parameter name:__ Coordinate system
**Default value:** cartesian

**Pattern:** [Selection cartesian|spherical|depth ]

**Documentation:** A selection that determines the assumed coordinate system for the function variables. Allowed values are &lsquo;cartesian&rsquo;, &lsquo;spherical&rsquo;, and &lsquo;depth&rsquo;. &lsquo;spherical&rsquo; coordinates are interpreted as r,phi or r,phi,theta in 2d/3d respectively with theta being the polar angle. &lsquo;depth&rsquo; will create a function, in which only the first parameter is non-zero, which is interpreted to be the depth of the point.

(parameters:Initial_20temperature_20model/Function/Function_20constants)=
### __Parameter name:__ Function constants
**Default value:**

**Pattern:** [Anything]

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form &lsquo;var1=value1, var2=value2, ...&rsquo;.

A typical example would be to set this runtime parameter to &lsquo;pi=3.1415926536&rsquo; and then use &lsquo;pi&rsquo; in the expression of the actual formula. (That said, for convenience this class actually defines both &lsquo;pi&rsquo; and &lsquo;Pi&rsquo; by default, but you get the idea.)

(parameters:Initial_20temperature_20model/Function/Function_20expression)=
### __Parameter name:__ Function expression
**Default value:** 0

**Pattern:** [Anything]

**Documentation:** The formula that denotes the function you want to evaluate for particular values of the independent variables. This expression may contain any of the usual operations such as addition or multiplication, as well as all of the common functions such as &lsquo;sin&rsquo; or &lsquo;cos&rsquo;. In addition, it may contain expressions like &lsquo;if(x>0, 1, -1)&rsquo; where the expression evaluates to the second argument if the first argument is true, and to the third argument otherwise. For a full overview of possible expressions accepted see the documentation of the muparser library at http://muparser.beltoforion.de/.

If the function you are describing represents a vector-valued function with multiple components, then separate the expressions for individual components by a semicolon.

(parameters:Initial_20temperature_20model/Function/Variable_20names)=
### __Parameter name:__ Variable names
**Default value:** x,y,t

**Pattern:** [Anything]

**Documentation:** The names of the variables as they will be used in the function, separated by commas. By default, the names of variables at which the function will be evaluated are &lsquo;x&rsquo; (in 1d), &lsquo;x,y&rsquo; (in 2d) or &lsquo;x,y,z&rsquo; (in 3d) for spatial coordinates and &lsquo;t&rsquo; for time. You can then use these variable names in your function expression and they will be replaced by the values of these variables at which the function is currently evaluated. However, you can also choose a different set of names for the independent variables at which to evaluate your function expression. For example, if you work in spherical coordinates, you may wish to set this input parameter to &lsquo;r,phi,theta,t&rsquo; and then use these variable names in your function expression.

(parameters:Initial_20temperature_20model/Harmonic_20perturbation)=
## **Subsection:** Initial temperature model / Harmonic perturbation
(parameters:Initial_20temperature_20model/Harmonic_20perturbation/Lateral_20wave_20number_20one)=
### __Parameter name:__ Lateral wave number one
**Default value:** 3

**Pattern:** [Integer range -2147483648...2147483647 (inclusive)]

**Documentation:** Doubled first lateral wave number of the harmonic perturbation. Equals the spherical harmonic degree in 3d spherical shells. In all other cases one equals half of a sine period over the model domain. This allows for single up-/downswings. Negative numbers reverse the sign of the perturbation but are not allowed for the spherical harmonic case.

(parameters:Initial_20temperature_20model/Harmonic_20perturbation/Lateral_20wave_20number_20two)=
### __Parameter name:__ Lateral wave number two
**Default value:** 2

**Pattern:** [Integer range -2147483648...2147483647 (inclusive)]

**Documentation:** Doubled second lateral wave number of the harmonic perturbation. Equals the spherical harmonic order in 3d spherical shells. In all other cases one equals half of a sine period over the model domain. This allows for single up-/downswings. Negative numbers reverse the sign of the perturbation.

(parameters:Initial_20temperature_20model/Harmonic_20perturbation/Magnitude)=
### __Parameter name:__ Magnitude
**Default value:** 1.0

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The magnitude of the Harmonic perturbation.

(parameters:Initial_20temperature_20model/Harmonic_20perturbation/Reference_20temperature)=
### __Parameter name:__ Reference temperature
**Default value:** 1600.0

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The reference temperature that is perturbed by the harmonic function. Only used in incompressible models.

(parameters:Initial_20temperature_20model/Harmonic_20perturbation/Vertical_20wave_20number)=
### __Parameter name:__ Vertical wave number
**Default value:** 1

**Pattern:** [Integer range -2147483648...2147483647 (inclusive)]

**Documentation:** Doubled radial wave number of the harmonic perturbation.  One equals half of a sine period over the model domain.  This allows for single up-/downswings. Negative numbers  reverse the sign of the perturbation.

(parameters:Initial_20temperature_20model/Inclusion_20shape_20perturbation)=
## **Subsection:** Initial temperature model / Inclusion shape perturbation
(parameters:Initial_20temperature_20model/Inclusion_20shape_20perturbation/Ambient_20temperature)=
### __Parameter name:__ Ambient temperature
**Default value:** 1.0

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** The background temperature for the temperature field.

(parameters:Initial_20temperature_20model/Inclusion_20shape_20perturbation/Center_20X)=
### __Parameter name:__ Center X
**Default value:** 0.5

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** The X coordinate for the center of the shape.

(parameters:Initial_20temperature_20model/Inclusion_20shape_20perturbation/Center_20Y)=
### __Parameter name:__ Center Y
**Default value:** 0.5

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** The Y coordinate for the center of the shape.

(parameters:Initial_20temperature_20model/Inclusion_20shape_20perturbation/Center_20Z)=
### __Parameter name:__ Center Z
**Default value:** 0.5

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** The Z coordinate for the center of the shape. This is only necessary for three-dimensional fields.

(parameters:Initial_20temperature_20model/Inclusion_20shape_20perturbation/Inclusion_20gradient)=
### __Parameter name:__ Inclusion gradient
**Default value:** constant

**Pattern:** [Selection gaussian|linear|constant ]

**Documentation:** The gradient of the inclusion to be generated.

(parameters:Initial_20temperature_20model/Inclusion_20shape_20perturbation/Inclusion_20shape)=
### __Parameter name:__ Inclusion shape
**Default value:** circle

**Pattern:** [Selection square|circle ]

**Documentation:** The shape of the inclusion to be generated.

(parameters:Initial_20temperature_20model/Inclusion_20shape_20perturbation/Inclusion_20temperature)=
### __Parameter name:__ Inclusion temperature
**Default value:** 0.0

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** The temperature of the inclusion shape. This is only the true temperature in the case of the constant gradient. In all other cases, it gives one endpoint of the temperature gradient for the shape.

(parameters:Initial_20temperature_20model/Inclusion_20shape_20perturbation/Shape_20radius)=
### __Parameter name:__ Shape radius
**Default value:** 1.0

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The radius of the inclusion to be generated. For shapes with no radius (e.g. square), this will be the width, and for shapes with no width, this gives a general guideline for the size of the shape.

(parameters:Initial_20temperature_20model/Lithosphere_20Mask)=
## **Subsection:** Initial temperature model / Lithosphere Mask
(parameters:Initial_20temperature_20model/Lithosphere_20Mask/Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/initial-temperature/lithosphere-mask/

**Pattern:** [DirectoryName]

**Documentation:** The path to the LAB depth data file

(parameters:Initial_20temperature_20model/Lithosphere_20Mask/Depth_20specification_20method)=
### __Parameter name:__ Depth specification method
**Default value:** Value

**Pattern:** [Selection File|Value ]

**Documentation:** Method that is used to specify the depth of the lithosphere-asthenosphere boundary.

(parameters:Initial_20temperature_20model/Lithosphere_20Mask/LAB_20depth_20filename)=
### __Parameter name:__ LAB depth filename
**Default value:** LAB_CAM2016.txt

**Pattern:** [FileName (Type: input)]

**Documentation:** File from which the lithosphere-asthenosphere boundary depth data is read.

(parameters:Initial_20temperature_20model/Lithosphere_20Mask/Lithosphere_20temperature)=
### __Parameter name:__ Lithosphere temperature
**Default value:** 1600.

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The initial temperature within lithosphere, applied abovethe maximum lithosphere depth.

(parameters:Initial_20temperature_20model/Lithosphere_20Mask/Maximum_20lithosphere_20depth)=
### __Parameter name:__ Maximum lithosphere depth
**Default value:** 200000.0

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Units: \si{\meter}.The maximum depth of the lithosphere. The model will be NaNs below this depth.

(parameters:Initial_20temperature_20model/Patch_20on_20S40RTS)=
## **Subsection:** Initial temperature model / Patch on S40RTS
(parameters:Initial_20temperature_20model/Patch_20on_20S40RTS/Maximum_20grid_20depth)=
### __Parameter name:__ Maximum grid depth
**Default value:** 700000.0

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The maximum depth of the Vs ascii grid. The model will read in  Vs from S40RTS below this depth.

(parameters:Initial_20temperature_20model/Patch_20on_20S40RTS/Remove_20temperature_20heterogeneity_20down_20to_20specified_20depth)=
### __Parameter name:__ Remove temperature heterogeneity down to specified depth
**Default value:** -1.7976931348623157e+308

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** This will set the heterogeneity prescribed by the Vs ascii grid and S40RTS to zero down to the specified depth (in meters). Note that your resolution has to be adequate to capture this cutoff. For example if you specify a depth of 660 km, but your closest spherical depth layers are only at 500 km and 750 km (due to a coarse resolution) it will only zero out heterogeneities down to 500 km. Similar caution has to be taken when using adaptive meshing.

(parameters:Initial_20temperature_20model/Patch_20on_20S40RTS/Smoothing_20length_20scale)=
### __Parameter name:__ Smoothing length scale
**Default value:** 200000.0

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The depth range (above maximum grid depth) over which to smooth. The boundary is smoothed using a depth weighted combination of Vs values from the ascii grid and S40RTS at each point in the region of smoothing.

(parameters:Initial_20temperature_20model/Patch_20on_20S40RTS/Ascii_20data_20model)=
## **Subsection:** Initial temperature model / Patch on S40RTS / Ascii data model
(parameters:Initial_20temperature_20model/Patch_20on_20S40RTS/Ascii_20data_20model/Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/initial-temperature/patch-on-S40RTS/test/

**Pattern:** [DirectoryName]

**Documentation:** The name of a directory that contains the model data. This path may either be absolute (if starting with a &lsquo;/&rsquo;) or relative to the current directory. The path may also include the special text &lsquo;$ASPECT_SOURCE_DIR&rsquo; which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the &lsquo;data/&rsquo; subdirectory of ASPECT.

(parameters:Initial_20temperature_20model/Patch_20on_20S40RTS/Ascii_20data_20model/Data_20file_20name)=
### __Parameter name:__ Data file name
**Default value:** upper_shell_3d.txt

**Pattern:** [Anything]

**Documentation:** The file name of the model data.

(parameters:Initial_20temperature_20model/Patch_20on_20S40RTS/Ascii_20data_20model/Scale_20factor)=
### __Parameter name:__ Scale factor
**Default value:** 1.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Scalar factor, which is applied to the model data. You might want to use this to scale the input to a reference model. Another way to use this factor is to convert units of the input files. For instance, if you provide velocities in cm/yr set this factor to 0.01.

(parameters:Initial_20temperature_20model/Random_20Gaussian_20perturbation)=
## **Subsection:** Initial temperature model / Random Gaussian perturbation
(parameters:Initial_20temperature_20model/Random_20Gaussian_20perturbation/Maximum_20magnitude)=
### __Parameter name:__ Maximum magnitude
**Default value:** 25.0

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The maximum magnitude of the Gaussian perturbation. For each perturbation, a random magnitude between plus and minus the maximum magnitude will be chosen. Units: \si{\kelvin}.

(parameters:Initial_20temperature_20model/Random_20Gaussian_20perturbation/Number_20of_20perturbations)=
### __Parameter name:__ Number of perturbations
**Default value:** 100

**Pattern:** [Integer range -2147483648...2147483647 (inclusive)]

**Documentation:** Total number of perturbations to be introduced into the model. Perturbations will be placed at random locations within the model domain.

(parameters:Initial_20temperature_20model/Random_20Gaussian_20perturbation/Width)=
### __Parameter name:__ Width
**Default value:** 1000.0

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The Gaussian RMS width of the perturbations. Units: \si{\meter}.

(parameters:Initial_20temperature_20model/S40RTS_20perturbation)=
## **Subsection:** Initial temperature model / S40RTS perturbation
(parameters:Initial_20temperature_20model/S40RTS_20perturbation/Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/initial-temperature/S40RTS/

**Pattern:** [DirectoryName]

**Documentation:** The path to the model data.

(parameters:Initial_20temperature_20model/S40RTS_20perturbation/Initial_20condition_20file_20name)=
### __Parameter name:__ Initial condition file name
**Default value:** S40RTS.sph

**Pattern:** [Anything]

**Documentation:** The file name of the spherical harmonics coefficients from Ritsema et al.

(parameters:Initial_20temperature_20model/S40RTS_20perturbation/Maximum_20degree)=
### __Parameter name:__ Maximum degree
**Default value:** 20

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** The maximum degree the users specify when reading the data file of spherical harmonic coefficients, which must be smaller than the maximum degree the data file stored. This parameter will be used only if &rsquo;Specify a lower maximum degree&rsquo; is set to true.

(parameters:Initial_20temperature_20model/S40RTS_20perturbation/Reference_20temperature)=
### __Parameter name:__ Reference temperature
**Default value:** 1600.0

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The reference temperature that is perturbed by the spherical harmonic functions. Only used in incompressible models.

(parameters:Initial_20temperature_20model/S40RTS_20perturbation/Remove_20degree_200_20from_20perturbation)=
### __Parameter name:__ Remove degree 0 from perturbation
**Default value:** true

**Pattern:** [Bool]

**Documentation:** Option to remove the degree zero component from the perturbation, which will ensure that the laterally averaged temperature for a fixed depth is equal to the background temperature.

(parameters:Initial_20temperature_20model/S40RTS_20perturbation/Remove_20temperature_20heterogeneity_20down_20to_20specified_20depth)=
### __Parameter name:__ Remove temperature heterogeneity down to specified depth
**Default value:** -1.7976931348623157e+308

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** This will set the heterogeneity prescribed by S20RTS or S40RTS to zero down to the specified depth (in meters). Note that your resolution has to be adequate to capture this cutoff. For example if you specify a depth of 660 km, but your closest spherical depth layers are only at 500 km and 750 km (due to a coarse resolution) it will only zero out heterogeneities down to 500 km. Similar caution has to be taken when using adaptive meshing.

(parameters:Initial_20temperature_20model/S40RTS_20perturbation/Specify_20a_20lower_20maximum_20degree)=
### __Parameter name:__ Specify a lower maximum degree
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Option to use a lower maximum degree when reading the data file of spherical harmonic coefficients. This is probably used for the faster tests or when the users only want to see the spherical harmonic pattern up to a certain degree.

(parameters:Initial_20temperature_20model/S40RTS_20perturbation/Spline_20knots_20depth_20file_20name)=
### __Parameter name:__ Spline knots depth file name
**Default value:** Spline_knots.txt

**Pattern:** [Anything]

**Documentation:** The file name of the spline knot locations from Ritsema et al.

(parameters:Initial_20temperature_20model/S40RTS_20perturbation/Thermal_20expansion_20coefficient_20in_20initial_20temperature_20scaling)=
### __Parameter name:__ Thermal expansion coefficient in initial temperature scaling
**Default value:** 2e-5

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The value of the thermal expansion coefficient $\beta$. Units: \si{\per\kelvin}.

(parameters:Initial_20temperature_20model/S40RTS_20perturbation/Use_20thermal_20expansion_20coefficient_20from_20material_20model)=
### __Parameter name:__ Use thermal expansion coefficient from material model
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Option to take the thermal expansion coefficient from the material model instead of from what is specified in this section.

(parameters:Initial_20temperature_20model/S40RTS_20perturbation/Vs_20to_20density_20scaling)=
### __Parameter name:__ Vs to density scaling
**Default value:** 0.25

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** This parameter specifies how the perturbation in shear wave velocity as prescribed by S20RTS or S40RTS is scaled into a density perturbation. See the general description of this model for more detailed information.

(parameters:Initial_20temperature_20model/S40RTS_20perturbation/Vs_20to_20density_20scaling_20method)=
### __Parameter name:__ Vs to density scaling method
**Default value:** constant

**Pattern:** [Selection file|constant ]

**Documentation:** Method that is used to specify how the vs-to-density scaling varies with depth.

(parameters:Initial_20temperature_20model/S40RTS_20perturbation/Ascii_20data_20vs_20to_20density_20model)=
## **Subsection:** Initial temperature model / S40RTS perturbation / Ascii data vs to density model
(parameters:Initial_20temperature_20model/S40RTS_20perturbation/Ascii_20data_20vs_20to_20density_20model/Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/initial-temperature/S40RTS/

**Pattern:** [DirectoryName]

**Documentation:** The name of a directory that contains the model data. This path may either be absolute (if starting with a &lsquo;/&rsquo;) or relative to the current directory. The path may also include the special text &lsquo;$ASPECT_SOURCE_DIR&rsquo; which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the &lsquo;data/&rsquo; subdirectory of ASPECT.

(parameters:Initial_20temperature_20model/S40RTS_20perturbation/Ascii_20data_20vs_20to_20density_20model/Data_20file_20name)=
### __Parameter name:__ Data file name
**Default value:** vs_to_density_Steinberger.txt

**Pattern:** [Anything]

**Documentation:** The file name of the model data.

(parameters:Initial_20temperature_20model/S40RTS_20perturbation/Ascii_20data_20vs_20to_20density_20model/Scale_20factor)=
### __Parameter name:__ Scale factor
**Default value:** 1.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Scalar factor, which is applied to the model data. You might want to use this to scale the input to a reference model. Another way to use this factor is to convert units of the input files. For instance, if you provide velocities in cm/yr set this factor to 0.01.

(parameters:Initial_20temperature_20model/SAVANI_20perturbation)=
## **Subsection:** Initial temperature model / SAVANI perturbation
(parameters:Initial_20temperature_20model/SAVANI_20perturbation/Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/initial-temperature/SAVANI/

**Pattern:** [DirectoryName]

**Documentation:** The path to the model data.

(parameters:Initial_20temperature_20model/SAVANI_20perturbation/Initial_20condition_20file_20name)=
### __Parameter name:__ Initial condition file name
**Default value:** savani.dlnvs.60.m.ab

**Pattern:** [Anything]

**Documentation:** The file name of the spherical harmonics coefficients from Auer et al.

(parameters:Initial_20temperature_20model/SAVANI_20perturbation/Maximum_20degree)=
### __Parameter name:__ Maximum degree
**Default value:** 20

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** The maximum degree the users specify when reading the data file of spherical harmonic coefficients, which must be smaller than the maximum degree the data file stored. This parameter will be used only if &rsquo;Specify a lower maximum degree&rsquo; is set to true.

(parameters:Initial_20temperature_20model/SAVANI_20perturbation/Reference_20temperature)=
### __Parameter name:__ Reference temperature
**Default value:** 1600.0

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The reference temperature that is perturbed by the spherical harmonic functions. Only used in incompressible models.

(parameters:Initial_20temperature_20model/SAVANI_20perturbation/Remove_20degree_200_20from_20perturbation)=
### __Parameter name:__ Remove degree 0 from perturbation
**Default value:** true

**Pattern:** [Bool]

**Documentation:** Option to remove the degree zero component from the perturbation, which will ensure that the laterally averaged temperature for a fixed depth is equal to the background temperature.

(parameters:Initial_20temperature_20model/SAVANI_20perturbation/Remove_20temperature_20heterogeneity_20down_20to_20specified_20depth)=
### __Parameter name:__ Remove temperature heterogeneity down to specified depth
**Default value:** -1.7976931348623157e+308

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** This will set the heterogeneity prescribed by SAVANI to zero down to the specified depth (in meters). Note that your resolution has to be adequate to capture this cutoff. For example if you specify a depth of 660 km, but your closest spherical depth layers are only at 500 km and 750 km (due to a coarse resolution) it will only zero out heterogeneities down to 500 km. Similar caution has to be taken when using adaptive meshing.

(parameters:Initial_20temperature_20model/SAVANI_20perturbation/Specify_20a_20lower_20maximum_20degree)=
### __Parameter name:__ Specify a lower maximum degree
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Option to use a lower maximum degree when reading the data file of spherical harmonic coefficients. This is probably used for the faster tests or when the users only want to see the spherical harmonic pattern up to a certain degree.

(parameters:Initial_20temperature_20model/SAVANI_20perturbation/Spline_20knots_20depth_20file_20name)=
### __Parameter name:__ Spline knots depth file name
**Default value:** Spline_knots.txt

**Pattern:** [Anything]

**Documentation:** The file name of the spline knots taken from the 28 spherical layers of SAVANI tomography model.

(parameters:Initial_20temperature_20model/SAVANI_20perturbation/Thermal_20expansion_20coefficient_20in_20initial_20temperature_20scaling)=
### __Parameter name:__ Thermal expansion coefficient in initial temperature scaling
**Default value:** 2e-5

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The value of the thermal expansion coefficient $\beta$. Units: \si{\per\kelvin}.

(parameters:Initial_20temperature_20model/SAVANI_20perturbation/Use_20thermal_20expansion_20coefficient_20from_20material_20model)=
### __Parameter name:__ Use thermal expansion coefficient from material model
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Option to take the thermal expansion coefficient from the material model instead of from what is specified in this section.

(parameters:Initial_20temperature_20model/SAVANI_20perturbation/Vs_20to_20density_20scaling)=
### __Parameter name:__ Vs to density scaling
**Default value:** 0.25

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** This parameter specifies how the perturbation in shear wave velocity as prescribed by SAVANI is scaled into a density perturbation. See the general description of this model for more detailed information.

(parameters:Initial_20temperature_20model/SAVANI_20perturbation/Vs_20to_20density_20scaling_20method)=
### __Parameter name:__ Vs to density scaling method
**Default value:** constant

**Pattern:** [Selection file|constant ]

**Documentation:** Method that is used to specify how the vs-to-density scaling varies with depth.

(parameters:Initial_20temperature_20model/SAVANI_20perturbation/Ascii_20data_20vs_20to_20density_20model)=
## **Subsection:** Initial temperature model / SAVANI perturbation / Ascii data vs to density model
(parameters:Initial_20temperature_20model/SAVANI_20perturbation/Ascii_20data_20vs_20to_20density_20model/Data_20directory)=
### __Parameter name:__ Data directory
**Default value:** $ASPECT_SOURCE_DIR/data/initial-temperature/S40RTS/

**Pattern:** [DirectoryName]

**Documentation:** The name of a directory that contains the model data. This path may either be absolute (if starting with a &lsquo;/&rsquo;) or relative to the current directory. The path may also include the special text &lsquo;$ASPECT_SOURCE_DIR&rsquo; which will be interpreted as the path in which the ASPECT source files were located when ASPECT was compiled. This interpretation allows, for example, to reference files located in the &lsquo;data/&rsquo; subdirectory of ASPECT.

(parameters:Initial_20temperature_20model/SAVANI_20perturbation/Ascii_20data_20vs_20to_20density_20model/Data_20file_20name)=
### __Parameter name:__ Data file name
**Default value:** vs_to_density_Steinberger.txt

**Pattern:** [Anything]

**Documentation:** The file name of the model data.

(parameters:Initial_20temperature_20model/SAVANI_20perturbation/Ascii_20data_20vs_20to_20density_20model/Scale_20factor)=
### __Parameter name:__ Scale factor
**Default value:** 1.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Scalar factor, which is applied to the model data. You might want to use this to scale the input to a reference model. Another way to use this factor is to convert units of the input files. For instance, if you provide velocities in cm/yr set this factor to 0.01.

(parameters:Initial_20temperature_20model/Spherical_20gaussian_20perturbation)=
## **Subsection:** Initial temperature model / Spherical gaussian perturbation
(parameters:Initial_20temperature_20model/Spherical_20gaussian_20perturbation/Amplitude)=
### __Parameter name:__ Amplitude
**Default value:** 0.01

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The amplitude of the perturbation.

(parameters:Initial_20temperature_20model/Spherical_20gaussian_20perturbation/Angle)=
### __Parameter name:__ Angle
**Default value:** 0.

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The angle where the center of the perturbation is placed.

(parameters:Initial_20temperature_20model/Spherical_20gaussian_20perturbation/Filename_20for_20initial_20geotherm_20table)=
### __Parameter name:__ Filename for initial geotherm table
**Default value:** initial-geotherm-table

**Pattern:** [FileName (Type: input)]

**Documentation:** The file from which the initial geotherm table is to be read. The format of the file is defined by what is read in source/initial\_temperature/spherical\_shell.cc.

(parameters:Initial_20temperature_20model/Spherical_20gaussian_20perturbation/Non_2ddimensional_20depth)=
### __Parameter name:__ Non-dimensional depth
**Default value:** 0.7

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The non-dimensional radial distance where the center of the perturbation is placed.

(parameters:Initial_20temperature_20model/Spherical_20gaussian_20perturbation/Sigma)=
### __Parameter name:__ Sigma
**Default value:** 0.2

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The standard deviation of the Gaussian perturbation.

(parameters:Initial_20temperature_20model/Spherical_20gaussian_20perturbation/Sign)=
### __Parameter name:__ Sign
**Default value:** 1.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** The sign of the perturbation.

(parameters:Initial_20temperature_20model/Spherical_20hexagonal_20perturbation)=
## **Subsection:** Initial temperature model / Spherical hexagonal perturbation
(parameters:Initial_20temperature_20model/Spherical_20hexagonal_20perturbation/Angular_20mode)=
### __Parameter name:__ Angular mode
**Default value:** 6

**Pattern:** [Integer range -2147483648...2147483647 (inclusive)]

**Documentation:** The number of convection cells with which to perturb the system.

(parameters:Initial_20temperature_20model/Spherical_20hexagonal_20perturbation/Rotation_20offset)=
### __Parameter name:__ Rotation offset
**Default value:** -45.

**Pattern:** [Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]

**Documentation:** Amount of clockwise rotation in degrees to apply to the perturbations. Default is set to -45 in order to provide backwards compatibility.
