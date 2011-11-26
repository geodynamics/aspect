/**
 * @defgroup Simulator The main class implementing the mantle convection problem
 *
 * This group contains the primary class implementing the mantle convection
 * problem as well as other things that may be directly related to it.
 */

/**
 * @defgroup Postprocessing Postprocessing solutions
 *
 * This group contains all classes, namespaces and functions that have to do
 * with postprocessing the solution of the simulation at the end of each time
 * step. This includes generating graphical output, computing integrals, or
 * other statistics.
 */

/**
 * @defgroup MaterialModels Describing the properties of the fluid
 *
 * This group contains all classes, namespaces and functions that have to do
 * with modeling the fluid under consideration. This includes handling all input
 * parameters related to the material, as well as describing how density, viscosity,
 * and many other parameters depend on pressure, temperature, and/or position within
 * the fluid.
 */

/**
 * @defgroup GeometryModels Describing the properties of the domain
 *
 * A module for the definition of properties of the geometry. This 
 * primarily includes the definition of the shape of the domain
 * (e.g. whether it is a full spherical shell, a quadrant/octant,
 * a description of the geoid, etc. The classes and functions of this
 * module also describe which kinds of boundary conditions hold on
 * the different parts of the boundary of the geometry.
 */

/**
 * @defgroup GravityModels Describing the properties of the gravity
 *
 * A module for the definition of gravity. Gravity is described as
 * a class that provides a function that returns the gravity vector
 * direction and magnitude as a function of location. The module also
 * contains the various functions necessary to add different gravity
 * models that can then be selected in the input file.
 */

/**
 * @defgroup InitialConditionsModels Describing initial conditions
 *
 * A module for the definition of functions and classes that have to
 * do with initial conditions for the temperature and, if necessary,
 * other quantities.
 */

/**
 * @defgroup BoundaryTemperatures Describing temperature boundary conditions
 *
 * A module for the definition of functions and classes that have to
 * do with describing boundary values for the temperature.
 */
