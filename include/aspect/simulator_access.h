/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


#ifndef _aspect_simulator_access_h
#define _aspect_simulator_access_h

#include <aspect/global.h>
#include <aspect/parameters.h>
#include <aspect/introspection.h>

#include <deal.II/base/table_handler.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/lac/affine_constraints.h>

namespace WorldBuilder
{
  class World;
}

namespace aspect
{
  // forward declarations:
  template <int dim> class Simulator;
  template <int dim> struct SimulatorSignals;
  template <int dim> class LateralAveraging;
  template <int dim> struct RotationProperties;

  namespace GravityModel
  {
    template <int dim> class Interface;
  }

  namespace HeatingModel
  {
    template <int dim> class Manager;
  }

  namespace MaterialModel
  {
    template <int dim> class Interface;
  }

  namespace InitialTemperature
  {
    template <int dim> class Manager;
    template <int dim> class Interface;
  }

  namespace BoundaryTemperature
  {
    template <int dim> class Manager;
    template <int dim> class Interface;
  }

  namespace BoundaryHeatFlux
  {
    template <int dim> class Interface;
  }

  namespace BoundaryConvectiveHeating
  {
    template <int dim> class Manager;
    template <int dim> class Interface;
  }

  namespace BoundaryComposition
  {
    template <int dim> class Manager;
    template <int dim> class Interface;
  }

  namespace BoundaryTraction
  {
    template <int dim> class Manager;
    template <int dim> class Interface;
  }

  namespace BoundaryVelocity
  {
    template <int dim> class Manager;
    template <int dim> class Interface;
  }

  namespace InitialComposition
  {
    template <int dim> class Manager;
    template <int dim> class Interface;
  }

  namespace InitialTopographyModel
  {
    template <int dim> class Interface;
  }

  namespace MeshRefinement
  {
    template <int dim> class Manager;
  }

  namespace AdiabaticConditions
  {
    template <int dim> class Interface;
  }

  namespace Postprocess
  {
    template <int dim> class Manager;
  }

  template <int dim> class MeltHandler;
  template <int dim> class VolumeOfFluidHandler;

  namespace MeshDeformation
  {
    template <int dim> class MeshDeformationHandler;
  }

  template <int dim> class NewtonHandler;

  template <int dim> class StokesMatrixFreeHandler;

  namespace Particle
  {
    template <int dim> class Manager;
  }

  namespace TimeStepping
  {
    template <int dim> class Manager;
  }

  /**
   * SimulatorAccess is a base class for different plugins like postprocessors.
   * It provides access to the various variables of the main class that
   * plugins may want to use in their evaluations, such as solution vectors,
   * the current time, time step sizes, material models, or the triangulations
   * and DoFHandlers that correspond to solutions.
   *
   * This class is the interface between plugins and the main simulator class.
   * Using this insulation layer, the plugins need not know anything about the
   * internal details of the simulation class.
   *
   * Every Postprocessor is required to derive from SimulatorAccess. It is
   * optional for other plugins like MaterialModel, GravityModel, etc..
   *
   * Since the functions providing access to details of the simulator class
   * are meant to be used only by derived classes of this class (rather than
   * becoming part of the public interface of these classes), the functions of
   * this class are made @p protected.
   *
   * @ingroup Simulator
   */
  template <int dim>
  class SimulatorAccess
  {
    public:
      /**
       * Default constructor. Initialize the SimulatorAccess object without
       * a reference to a particular Simulator object. You will later have
       * to call initialize() to provide this reference to the Simulator
       * object.
       */
      SimulatorAccess ();

      /**
       * Create a SimulatorAccess object that is already initialized for
       * a particular Simulator.
       */
      SimulatorAccess (const Simulator<dim> &simulator_object);

      /**
       * Destructor. Does nothing but is virtual so that derived classes
       * destructors are also virtual.
       */
      virtual ~SimulatorAccess () = default;

      /**
       * Initialize this class for a given simulator. This function is marked
       * as virtual so that derived classes can do something upon
       * initialization as well, for example look up and cache data; derived
       * classes should call this function from the base class as well,
       * however.
       *
       * @param simulator_object A reference to the main simulator object.
       */
      virtual void initialize_simulator (const Simulator<dim> &simulator_object);

      /** @name Accessing variables that identify overall properties of the simulator */
      /** @{ */

      /**
       * Return a reference to an introspection object that describes overall
       * properties of the simulator. In particular, it provides symbolic
       * names for extractors and component masks for each variable, etc, and
       * thereby reduces the need for implicit knowledge throughout the code
       * base.
       */
      const Introspection<dim> &
      introspection () const;

      /**
       * Return a reference to the Simulator itself. Note that you can not
       * access any members or functions of the Simulator. This function
       * exists so that any class with SimulatorAccess can create other
       * objects with SimulatorAccess (because initializing them requires a
       * reference to the Simulator).
       */
      const Simulator<dim> &
      get_simulator () const;

      /**
       * Return a reference to the parameters object that describes all run-time
       * parameters used in the current simulation.
       */
      const Parameters<dim> &
      get_parameters () const;

      /**
       * Get Access to the structure containing the signals of the simulator.
       */
      SimulatorSignals<dim> &
      get_signals() const;

      /**
       * Return the MPI communicator for this simulation.
       */
      MPI_Comm
      get_mpi_communicator () const;

      /**
       * Return the timer object for this simulation. Since the timer is
       * mutable in the Simulator class, this allows plugins to define their
       * own sections in the timer to measure the time spent in sections of
       * their code.
       */
      TimerOutput &
      get_computing_timer () const;

      /**
       * Return a reference to the stream object that only outputs something
       * on one processor in a parallel program and simply ignores output put
       * into it on all other processors.
       */
      const ConditionalOStream &
      get_pcout () const;

      /**
       * Return the current simulation time in seconds.
       */
      double get_time () const;

      /**
       * Return the size of the current time step.
       */
      double
      get_timestep () const;

      /**
       * Return the size of the last time step.
       */
      double
      get_old_timestep () const;

      /**
       * Return the current number of a time step.
       */
      unsigned int
      get_timestep_number () const;

      /**
       * Return a reference to the manager of the time stepping strategies.
       * This can then be used, for example, to check whether a checkpoint needs to
       * be made upon termination.
       */
      const TimeStepping::Manager<dim> &
      get_timestepping_manager() const;

      /**
       * Return the current nonlinear iteration number of a time step.
       */
      unsigned int
      get_nonlinear_iteration () const;

      /**
       * Return a reference to the triangulation in use by the simulator
       * object.
       */
      const parallel::distributed::Triangulation<dim> &
      get_triangulation () const;

      /**
       * Return the global volume of the computational domain.
       */
      double
      get_volume () const;

      /**
       * Return a reference to the mapping used to describe the boundary of
       * the domain.
       */
      const Mapping<dim> &
      get_mapping () const;

      /**
       * Return the directory specified in the input parameter file to be the
       * place where output files are to be placed. The string is terminated
       * by a directory separator (i.e., '/').
       */
      std::string
      get_output_directory () const;

      /**
       * Return whether we use the adiabatic heating term.
       */
      bool
      include_adiabatic_heating () const;

      /**
       * Return whether we use the latent heat term.
       */
      bool
      include_latent_heat () const;

      /**
       * Return whether we solve the equations for melt transport.
       */
      bool
      include_melt_transport () const;

      /**
       * Return the stokes velocity degree.
       */
      int
      get_stokes_velocity_degree () const;

      /**
       * Return the adiabatic surface temperature.
       */
      double
      get_adiabatic_surface_temperature () const;

      /**
       * Return the adiabatic surface pressure.
       */
      double
      get_surface_pressure () const;

      /**
       * Return whether things like velocities should be converted from the
       * seconds in the MKS system to years. The value of this flag is set by
       * the corresponding entry in the input parameter file.
       */
      bool
      convert_output_to_years () const;

      /**
       * Return the number of the current pre refinement step.
       * This can be useful for plugins that want to function differently in
       * the initial adaptive refinements and later on.
       * This will be not initialized before Simulator<dim>::run() is called.
       * It iterates upward from 0 to parameters.initial_adaptive_refinement
       * during the initial adaptive refinement steps, and equals
       * std::numeric_limits<unsigned int>::max() afterwards.
       */
      unsigned int
      get_pre_refinement_step () const;

      /**
       * Return the number of compositional fields specified in the input
       * parameter file that will be advected along with the flow field.
       */
      unsigned int
      n_compositional_fields () const;

      /**
       * Return the simulation start time in seconds. For most simulations,
       * this will be zero -- but it may be set to something different in
       * input files.
       */
      double
      get_start_time () const;

      /**
       * Return the simulation end time in seconds.
       */
      double
      get_end_time () const;

      /**
       * Compute the error indicators in the same way they are normally used
       * for mesh refinement. The mesh is not refined when doing so, but the
       * indicators can be used when generating graphical output to check why
       * mesh refinement is proceeding as it is.
       */
      void
      get_refinement_criteria(Vector<float> &estimated_error_per_cell) const;

      /**
       * Returns the entropy viscosity on each locally owned cell as it is
       * used to stabilize the temperature equation.
       *
       * @param viscosity_per_cell Output vector with as many entries as
       * active cells. Each entry corresponding to a locally owned active
       * cell index will contain the artificial viscosity for this cell.
       * @param skip_interior_cells A boolean flag. If set to true the function
       * will only compute the artificial viscosity in cells at boundaries.
       */
      void
      get_artificial_viscosity(Vector<float> &viscosity_per_cell,
                               const bool skip_interior_cells = false) const;

      /**
       * Returns the entropy viscosity on each locally owned cell as it is
       * used to stabilize the composition equation.
       */
      void
      get_artificial_viscosity_composition(Vector<float> &viscosity_per_cell,
                                           const unsigned int compositional_variable) const;
      /** @} */


      /** @name Accessing variables that identify the solution of the problem */
      /** @{ */

      /**
       * Return a reference to the vector that has the current linearization
       * point of the entire system, i.e. the velocity and pressure variables
       * as well as the temperature and compositional fields. This vector is
       * associated with the DoFHandler object returned by get_dof_handler().
       * This vector is only different from the one returned by get_solution()
       * during the solver phase.
       *
       * @note In general the vector is a distributed vector; however, it
       * contains ghost elements for all locally relevant degrees of freedom.
       */
      const LinearAlgebra::BlockVector &
      get_current_linearization_point () const;

      /**
       * Return a reference to the vector that has the current solution of the
       * entire system, i.e. the velocity and pressure variables as well as
       * the temperature and compositional fields.
       * This vector is associated with the DoFHandler object returned by
       * get_dof_handler().
       *
       * @note In general the vector is a distributed vector; however, it
       * contains ghost elements for all locally relevant degrees of freedom.
       */
      const LinearAlgebra::BlockVector &
      get_solution () const;

      /**
       * Return a reference to the vector that has the solution of the entire
       * system at the previous time step. This vector is associated with the
       * DoFHandler object returned by get_stokes_dof_handler().
       *
       * @note In general the vector is a distributed vector; however, it
       * contains ghost elements for all locally relevant degrees of freedom.
       */
      const LinearAlgebra::BlockVector &
      get_old_solution () const;

      /**
       * Return a reference to the vector that has the solution of the entire
       * system at the second-to-last time step. This vector is associated with the
       * DoFHandler object returned by get_stokes_dof_handler().
       *
       * @note In general the vector is a distributed vector; however, it
       * contains ghost elements for all locally relevant degrees of freedom.
       */
      const LinearAlgebra::BlockVector &
      get_old_old_solution () const;

      /**
       * Return a reference to the vector that has the reactions computed by the
       * operator splitting scheme in the current time step.
       *
       * @note In general the vector is a distributed vector; however, it
       * contains ghost elements for all locally relevant degrees of freedom.
       */
      const LinearAlgebra::BlockVector &
      get_reaction_vector () const;

      /**
       * Return a reference to the vector that has the mesh velocity for
       * simulations with mesh deformation.
       *
       * @note In general the vector is a distributed vector; however, it
       * contains ghost elements for all locally relevant degrees of freedom.
       */
      const LinearAlgebra::BlockVector &
      get_mesh_velocity () const;

      /**
       * Return a reference to the DoFHandler that is used to discretize the
       * variables at the current time step.
       */
      const DoFHandler<dim> &
      get_dof_handler () const;

      /**
       * Return a reference to the finite element that is
       * used to discretize the variables at the current time step.
       * This is the finite element for the entire, coupled problem, i.e.,
       * it contains sub-elements for velocity, pressure, temperature and all
       * other variables in this problem (e.g., compositional variables, if
       * used in this simulation).
       */
      const FiniteElement<dim> &
      get_fe () const;

      /**
       * Return a reference to the system matrix at the current time step.
       */
      const LinearAlgebra::BlockSparseMatrix &
      get_system_matrix () const;

      /**
       * Return a reference to the system preconditioner matrix at the current time step.
       */
      const LinearAlgebra::BlockSparseMatrix &
      get_system_preconditioner_matrix () const;

      /** @} */


      /** @name Accessing variables that identify aspects of the simulation */
      /** @{ */

      /**
       * Return a pointer to the material model to access functions like
       * density().
       */
      const MaterialModel::Interface<dim> &
      get_material_model () const;

      /**
       * Return a pointer to the gravity model description.
       */
      const GravityModel::Interface<dim> &
      get_gravity_model () const;

      /**
       * Return a pointer to the initial topography model.
       */
      const InitialTopographyModel::Interface<dim> &
      get_initial_topography_model () const;

      /**
       * Return a pointer to the geometry model.
       */
      const GeometryModel::Interface<dim> &
      get_geometry_model () const;


      /**
       * Return a pointer to the object that describes the adiabatic
       * conditions.
       */
      const AdiabaticConditions::Interface<dim> &
      get_adiabatic_conditions () const;

      /**
       * Return whether the current model has a boundary temperature object
       * set. This is useful because a simulation does not actually have to
       * declare any boundary temperature model, for example if all
       * boundaries are insulating. In such cases, there is no
       * boundary temperature model that can provide, for example,
       * a minimal and maximal temperature on the boundary.
       */
      bool has_boundary_temperature () const;

      /**
       * Return a reference to the manager of the boundary temperature models.
       * This can then, for example, be used to get the names of the boundary temperature
       * models used in a computation, or to compute the boundary temperature
       * for a given position.
       */
      const BoundaryTemperature::Manager<dim> &
      get_boundary_temperature_manager () const;

      /**
       * Return a reference to the manager of the boundary convective heating models.
       * This can then, for example, be used to get the names of the boundary convective
       * heating models used in a computation, or to compute the boundary conditions
       * for a given position.
       */
      const BoundaryConvectiveHeating::Manager<dim> &
      get_boundary_convective_heating_manager () const;

      /**
       * Return a reference to the object that describes heat flux
       * boundary conditions.
       */
      const BoundaryHeatFlux::Interface<dim> &
      get_boundary_heat_flux () const;

      /**
       * Return whether the current model has a boundary composition object
       * set. This is useful because a simulation does not actually have to
       * declare any boundary composition model, for example if all
       * boundaries are reflecting. In such cases, there is no
       * boundary composition model.
       */
      bool has_boundary_composition () const;

      /**
       * Return an reference to the manager of the boundary composition models.
       * This can then, for example, be used to get the names of the boundary composition
       * models used in a computation, or to compute the boundary composition
       * for a given position.
       */
      const BoundaryComposition::Manager<dim> &
      get_boundary_composition_manager () const;

      /**
       * Return an reference to the manager of the boundary traction models.
       * This can then, for example, be used to get the names of the boundary traction
       * models used in a computation, or to compute the boundary traction
       * for a given position.
       */
      const BoundaryTraction::Manager<dim> &
      get_boundary_traction_manager () const;

      /**
       * Return a reference to the manager of the initial temperature models.
       * This can then, for example, be used to get the names of the initial temperature
       * models used in a computation, or to compute the initial temperature
       * for a given position.
       *
       * While the Simulator class creates a shared pointer to an initial
       * temperature manager before the first time step, it releases
       * the pointer once it no longer needs access to the initial
       * compositions. As a consequence, you can only call this function
       * during the first time step.
       *
       * If the Simulator's shared pointer were the only
       * one that points to the initial temperature manager object, that
       * would also destroy the object pointed to. However, plugin classes
       * can have member variables that are *also* shared pointers to
       * these manager objects, and if you initialize such a shared
       * pointer from the result of this function -- typically in the
       * `initialize()` function of a plugin class -- then the Simulator
       * giving up its shared pointer does not actually destroy the
       * manager object but extends its lifetime until the last plugin
       * that has a pointer to it is destroyed itself. As a consequence,
       * if you need access to the initial temperature in a plugin, you
       * will need to keep a shared pointer to it around for as long
       * as you need it.
       */
      std::shared_ptr<const InitialTemperature::Manager<dim>>
      get_initial_temperature_manager_pointer () const;

      /**
       * Return a reference to the manager of the initial temperature model.
       * This can then, for example, be used to get the names of the initial temperature
       * models used in a computation.
       *
       * While the Simulator class creates a shared pointer to an initial
       * temperature manager before the first time step, it releases
       * the pointer once it no longer needs access to the initial
       * temperature. As a consequence, you can only call this function
       * during the first time step. If a plugin needs access to the initial
       * temperature at a later time, it has to store its own shared
       * pointer to that object, and that is what can be achieved using
       * the get_initial_temperature_manager_pointer() function above.
       */
      const InitialTemperature::Manager<dim> &
      get_initial_temperature_manager () const;

      /**
       * Return a pointer to the manager of the initial composition model.
       * This can then, for example, be used to get the names of the initial composition
       * models used in a computation.
       *
       * While the Simulator class creates a shared pointer to an initial
       * composition manager before the first time step, it releases
       * the pointer once it no longer needs access to the initial
       * compositions. As a consequence, you can only call this function
       * during the first time step.
       *
       * If the Simulator's shared pointer were the only
       * one that points to the initial composition manager object, that
       * would also destroy the object pointed to. However, plugin classes
       * can have member variables that are *also* shared pointers to
       * these manager objects, and if you initialize such a shared
       * pointer from the result of this function -- typically in the
       * `initialize()` function of a plugin class -- then the Simulator
       * giving up its shared pointer does not actually destroy the
       * manager object but extends its lifetime until the last plugin
       * that has a pointer to it is destroyed itself. As a consequence,
       * if you need access to the initial compositions in a plugin, you
       * will need to keep a shared pointer to it around for as long
       * as you need it.
       */
      std::shared_ptr<const InitialComposition::Manager<dim>>
      get_initial_composition_manager_pointer () const;

      /**
       * Return a reference to the manager of the initial composition model.
       * This can then, for example, be used to get the names of the initial composition
       * models used in a computation.
       *
       * While the Simulator class creates a shared pointer to an initial
       * composition manager before the first time step, it releases
       * the pointer once it no longer needs access to the initial
       * compositions. As a consequence, you can only call this function
       * during the first time step. If a plugin needs access to the initial
       * composition at a later time, it has to store its own shared
       * pointer to that object, and that is what can be achieved using
       * the get_initial_composition_manager_pointer() function above.
       */
      const InitialComposition::Manager<dim> &
      get_initial_composition_manager () const;

      /**
       * Return a set of boundary indicators that describes which of the
       * boundaries have a fixed temperature.
       */
      const std::set<types::boundary_id> &
      get_fixed_temperature_boundary_indicators () const;

      /**
       * Return a set of boundary indicators that describes which of the
       * boundaries have a fixed heat flux.
       */
      const std::set<types::boundary_id> &
      get_fixed_heat_flux_boundary_indicators () const;

      /**
       * Return a set of boundary indicators that describes which of the
       * boundaries have a convective heating, i.e. Robin boundary condition.
       */
      const std::set<types::boundary_id> &
      get_fixed_convective_heating_boundary_indicators () const;

      /**
       * Return a set of boundary indicators that describes which of the
       * boundaries have a fixed composition.
       */
      const std::set<types::boundary_id> &
      get_fixed_composition_boundary_indicators () const;

      /**
       * Return a set of boundary indicators that describes which of the
       * boundaries have a mesh deformation boundary condition. Note that
       * it does not specify which boundaries have which mesh deformation
       * condition, only which boundaries have a mesh deformation condition.
       */
      const std::set<types::boundary_id> &
      get_mesh_deformation_boundary_indicators () const;

      /**
       * Return an reference to the manager of the boundary velocity models.
       * This can then, for example, be used to get the names of the boundary velocity
       * models used in a computation, or to compute the boundary velocity
       * for a given position.
       */
      const BoundaryVelocity::Manager<dim> &
      get_boundary_velocity_manager () const;

      /**
       * Return a pointer to the manager of the heating model.
       * This can then, for example, be used to get the names of the heating models
       * used in a computation.
       */
      const HeatingModel::Manager<dim> &
      get_heating_model_manager () const;

      /**
       * Return a reference to the manager of the mesh refinement strategies.
       * this can then, for example, be used to get the names of the active refinement
       * strategies for such purposes as confirming that a particular one has
       * been included.
       */
      const MeshRefinement::Manager<dim> &
      get_mesh_refinement_manager () const;

      /**
       * Return a reference to the melt handler.
       */
      const MeltHandler<dim> &
      get_melt_handler () const;

      /**
       * Return a reference to the VolumeOfFluid handler.
       */
      const VolumeOfFluidHandler<dim> &
      get_volume_of_fluid_handler () const;

      /**
       * Return a reference to the Newton handler that controls the Newton
       * iteration to resolve nonlinearities.
       */
      const NewtonHandler<dim> &
      get_newton_handler () const;

#ifdef ASPECT_WITH_WORLD_BUILDER
      /**
       * Return a reference to the world builder that controls the setup of
       * initial conditions.
       *
       * This call will only succeed if ASPECT was configured to use
       * the WorldBuilder.
       *
       * While the Simulator class creates a shared pointer to a
       * WorldBuilder object before the first time step, it releases
       * the pointer once it no longer needs access to the initial
       * conditions. As a consequence, you can only call this function
       * during the first time step. If a plugin needs access to the object
       * so returned at a later time, it has to store its own shared
       * pointer to that object, and that is what can be achieved using
       * the get_world_builder_pointer() function below.
       */
      const WorldBuilder::World &
      get_world_builder () const;

      /**
       * This function is to get_world_builder() what
       * get_initial_temperature_manager_pointer() is to
       * the get_initial_temperature_manager() function: It returns a
       * shared pointer so that objects that still need access to the
       * WorldBuilder object after the Simulator class has released
       * it, can extend the lifetime of the object pointed to by
       * keeping a shared pointer to it.
       */
      std::shared_ptr<const WorldBuilder::World>
      get_world_builder_pointer () const;
#endif
      /**
       * Return a reference to the mesh deformation handler. This function will
       * throw an exception if mesh deformation is not activated.
       */
      const MeshDeformation::MeshDeformationHandler<dim> &
      get_mesh_deformation_handler () const;

      /**
       * Return a reference to the lateral averaging object owned
       * by the simulator, which can be used to query lateral averages
       * of various quantities at depth slices.
       */
      const LateralAveraging<dim> &
      get_lateral_averaging () const;

      /**
       * Return a pointer to the object that describes the DoF
       * constraints for the time step we are currently solving.
       */
      const AffineConstraints<double> &
      get_current_constraints () const;

      /**
       * Return whether the Simulator object has been completely initialized
       * and has started to run its time stepping loop.
       *
       * This function is useful to determine in a plugin whether some
       * of the information one can query about the Simulator can be trusted
       * because it has already been set up completely. For example,
       * while the Simulator is being
       * set up, plugins may already have access to it via the current
       * SimulatorAccess object, but data such as the current time, the
       * time step number, etc, may all still be in a state that is not
       * reliable since it may not have been initialized at that time. (As
       * an example, at the very beginning of the Simulator object's existence,
       * the time step number is set to numbers::invalid_unsigned_int, and
       * only when the time step loop is started is it set to a valid
       * value). Similar examples are that at some point the Simulator
       * sets the solution vector to the correct size, but only at a later
       * time (though before the time stepping starts), the *contents* of
       * the solution vector are set based on the initial conditions
       * specified in the input file.
       *
       * Only when this function returns @p true is all of the information
       * returned by the SimulatorAccess object reliable and correct.
       *
       * @note This function returns @p true starting with the moment where the
       *   Simulator starts the time stepping loop. However, it may
       *   temporarily revert to returning @p false if, for example,
       *   the Simulator does the initial mesh refinement steps where
       *   it starts the time loop, but then goes back to
       *   initialization steps (mesh refinement, interpolation of initial
       *   conditions, etc.) before re-starting the time loop.
       */
      bool simulator_is_past_initialization () const;

      /**
       * Return the value used for rescaling the pressure in the linear
       * solver.
       */
      double
      get_pressure_scaling () const;

      /**
       * Return whether we need to apply a compatibility modification
       * to the pressure right hand side. See documentation of
       * Simulator<dim>::do_pressure_rhs_compatibility_modification for more
       * information.
       */
      bool
      pressure_rhs_needs_compatibility_modification() const;

      /**
       * Return whether the model uses a prescribed Stokes solution.
       */
      bool
      model_has_prescribed_stokes_solution () const;

      /**
       * A convenience function that copies the values of the compositional
       * fields at the quadrature point q given as input parameter to the
       * output vector composition_values_at_q_point.
       */
      static
      void
      get_composition_values_at_q_point (const std::vector<std::vector<double>> &composition_values,
                                         const unsigned int                      q,
                                         std::vector<double>                    &composition_values_at_q_point);

      /**
       * Return a writable reference to the statistics object into which
       * you can store additional data that then shows up in the
       * <code>output_dir/statistics</code> file.
       *
       * Postprocessor objects get a reference to this object automatically
       * when called, but other plugins may not. They do not usually
       * produce output anyway, but through this function they can still
       * record information as necessary.
       * @return
       */
      TableHandler &get_statistics_object() const;

      /**
       * Return a reference to the melt handler.
       */
      const Postprocess::Manager<dim> &
      get_postprocess_manager () const;

      /**
       * Returns the number of active particle managers.
       */
      unsigned int
      n_particle_managers() const;

      /**
       * Returns a const reference to a single particle manager given the
       * index.
       */
      const Particle::Manager<dim> &
      get_particle_manager(const unsigned int particle_manager_index) const;

      /**
       * Returns a reference to a single particle manager, in case anyone wants to
       * change something within the particle manager. Use with care, usually
       * you want to only let the functions within the particle subsystem
       * change member variables of the particle manager.
       */
      Particle::Manager<dim> &
      get_particle_manager(const unsigned int particle_manager_index);

      /**
       * Return true if using the block GMG Stokes solver.
       */
      bool is_stokes_matrix_free();

      /**
       * Return a reference to the StokesMatrixFreeHandler that controls the
       * matrix-free Stokes solver.
       */
      const StokesMatrixFreeHandler<dim> &
      get_stokes_matrix_free () const;

      /**
       * Compute the angular momentum and other rotation properties
       * of the velocities in the given solution vector.
       *
       * @param use_constant_density determines whether to use a constant
       * density (which corresponds to computing a net rotation instead of net
       * angular momentum).
       * @param solution Solution vector to compute the properties for.
       * @param limit_to_top_faces allows to only compute the net angular momentum
       * (or net rotation) of the top surface.
       */
      RotationProperties<dim>
      compute_net_angular_momentum(const bool use_constant_density,
                                   const LinearAlgebra::BlockVector &solution,
                                   const bool limit_to_top_faces = false) const;

      /**
      * Eliminate the nullspace of the velocity in the given vector. Both
      * vectors are expected to contain the current solution.
      *
      * @param solution The locally relevant vector for the whole
      * finite element, this vector will be filled at the end.
      * @param distributed_stokes_solution only contains velocity and pressure and
      * only locally owned elements.
      */
      void remove_nullspace(LinearAlgebra::BlockVector &solution,
                            LinearAlgebra::BlockVector &distributed_stokes_solution) const;

      /**
       * Adjust the pressure variable (which is only determined up to
       * a constant by the equations) by adding a constant to it in
       * such a way that the pressure on the surface or within the
       * entire volume has a known average value. See the documentation
       * of the normalize_pressure() function in the Simulator class
       * for more information.
       *
       * @return This function returns the pressure adjustment by value.
       * This is so that its negative can later be used again in
       * denormalize_pressure().
       */
      double normalize_pressure(LinearAlgebra::BlockVector &vector) const;

      /**
       * Invert the action of the normalize_pressure() function above. This
       * means that we move from a pressure that satisfies the pressure
       * normalization (e.g., has a zero average pressure, or a zero average
       * surface pressure) to one that does not actually satisfy this
       * normalization. See the function denormalize_pressure() in the
       * Simulator class for more information.
       *
       * This function modifies @p vector in-place.
       */
      void denormalize_pressure(const double                      pressure_adjustment,
                                LinearAlgebra::BlockVector       &vector) const;

      /** @} */

    private:
      /**
       * A pointer to the simulator object to which we want to get access.
       */
      const Simulator<dim> *simulator;
  };
}


#endif
