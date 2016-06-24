/**
 * @defgroup Simulator The main class implementing the mantle convection
 * problem
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
 * @defgroup Visualization Visualizing solutions
 *
 * This group contains all classes, namespaces and functions that have to do
 * with visualizing the solution of the simulation at the end of each time
 * step. Visualization is a part of postprocessing solutions.
 *
 * @ingroup Postprocessing
 */

/**
 * @defgroup MaterialModels Describing the material model
 *
 * This group contains all classes, namespaces and functions that have to do
 * with modeling the material properties of the fluid under consideration, as
 * well as all of the other variables in use such as compositional
 * fields. This includes handling all input parameters related to the
 * material, as well as describing how density, viscosity, and many other
 * parameters depend on pressure, temperature, and/or position within the
 * fluid.
 */

/**
 * @defgroup GeometryModels Describing the properties of the domain
 *
 * A module for the definition of properties of the geometry. This primarily
 * includes the definition of the shape of the domain (e.g. whether it is a
 * full spherical shell, a quadrant/octant, a description of the geoid, etc.
 * The classes and functions of this module also describe which kinds of
 * boundary conditions hold on the different parts of the boundary of the
 * geometry.
 */

/**
 * @defgroup GravityModels Describing the properties of the gravity
 *
 * A module for the definition of gravity. Gravity is described as a class
 * that provides a function that returns the gravity vector direction and
 * magnitude as a function of location. The module also contains the various
 * functions necessary to add different gravity models that can then be
 * selected in the input file.
 */

/**
 * @defgroup HeatingModels Describing terms on the right hand side of the temperature equation
 *
 * A module for the definition of terms that appear as heating or cooling terms
 * on the right hand side of the temperature equation. This includes, in
 * particular, terms such as shear and adiabatic heating, phase changes,
 * and radiogenic sources.
 */

/**
 * @defgroup MeshRefinement Describing how to refine the mesh
 *
 * A module for the definition of algorithms and methods to refine the
 * mesh adaptively.
 */

/**
 * @defgroup TerminationCriteria Describing when to terminate a simulation
 *
 * A module for the definition of algorithms and methods that determine
 * when to stop a simulation.
 */



/**
 * @defgroup InitialConditions Describing initial conditions
 *
 * A group for everything that has to do with the description of the
 * initial state of a model.
 */

/**
 * @defgroup InitialConditionsModels Describing initial conditions for the temperature
 *
 * A module for the definition of functions and classes that have to do with
 * initial conditions for the temperature.
 *
 * @ingroup InitialConditions
 */

/**
 * @defgroup CompositionalInitialConditionsModels Describing initial conditions for compositional fields
 *
 * A module for the definition of functions and classes that have to do with
 * initial conditions for the temperature.
 *
 * @ingroup InitialConditions
 */


/**
 * @defgroup AdiabaticConditions Adiabatic conditions
 *
 * A group for all things that have to do with adiabatic conditions.
 *
 * @ingroup InitialConditions
 */

/**
 * @defgroup PrescribedStokesSolution Describing prescribed velocity fields
 *
 * A group for all things that have to do with prescribing a velocity
 * field for simulations that only compute advected temperature and
 * compositional fields, but do not compute a Stokes solution.
 *
 * @ingroup InitialConditions
 */



/**
 * @defgroup BoundaryConditions Describing boundary conditions
 *
 * A group for everything that has to do with the description of boundary
 * values for the problem.
 */


/**
 * @defgroup BoundaryTemperatures Describing temperature boundary conditions
 *
 * A module for the definition of functions and classes that have to do with
 * describing boundary values for the temperature.
 *
 * @ingroup BoundaryConditions
 */

/**
 * @defgroup BoundaryCompositions Describing boundary conditions for compositional fields
 *
 * A module for the definition of functions and classes that have to do with
 * describing boundary values for compositional fields.
 *
 * @ingroup BoundaryConditions
 */

/**
 * @defgroup VelocityBoundaryConditions Describing boundary conditions for the velocity field
 *
 * A module for the definition of functions and classes that have to do with
 * describing boundary values for the velocity field.
 *
 * @ingroup BoundaryConditions
 */

/**
 * @defgroup TractionBoundaryConditionsModels Describing traction boundary conditions for the velocity field
 *
 * A module for the definition of functions and classes that have to do with
 * describing traction boundary values for the velocity field.
 *
 * @ingroup BoundaryConditions
 */

/**
 * @defgroup Particle Describing advected particles
 *
 * A module that contains everything related to particles.
 */

/**
 * @defgroup ParticleGenerators Describing the method to generate particles
 *
 * A module for the definition of functions and classes that have to do with
 * describing a method to generate particles in the model domain.
 *
 * @ingroup Particle
 */

/**
 * @defgroup ParticleIntegrators Describing the method to integrate particles
 *
 * A module for the definition of functions and classes that have to do with
 * integrating the particle positions through time.
 *
 * @ingroup Particle
 */

/**
 * @defgroup ParticleInterpolators Describing the method to interpolate particle properties to arbitrary positions
 *
 * A module for the definition of functions and classes that have to do with
 * interpolating particle properties to arbitrary positions.
 *
 * @ingroup Particle
 */

/**
 * @defgroup ParticleOutput Describing the output method of particles
 *
 * A module for the definition of functions and classes that have to do with
 * describing a method to output the particles.
 *
 * @ingroup Particle
 */

/**
 * @defgroup ParticleProperties Describing the properties of particles
 *
 * A module for the definition of functions and classes that have to do with
 * describing the properties of the advected particles.
 *
 * @ingroup Particle
 */
