/*
  Copyright (C) 2011, 2012 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/
/*  $Id$  */


#ifndef __aspect__simulator_h
#define __aspect__simulator_h

#include <deal.II/base/timer.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/symmetric_tensor.h>

#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_precondition.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>

#include <aspect/global.h>
#include <aspect/material_model/interface.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/boundary_temperature/interface.h>
#include <aspect/initial_conditions/interface.h>
#include <aspect/compositional_initial_conditions/interface.h>
#include <aspect/velocity_boundary_conditions/interface.h>
#include <aspect/postprocess/interface.h>
#include <aspect/adiabatic_conditions.h>



namespace aspect
{
  using namespace dealii;

  namespace internal
  {
    namespace Assembly
    {
      namespace Scratch
      {
        template <int dim>      struct StokesPreconditioner;
        template <int dim>      struct StokesSystem;
        template <int dim>      struct TemperatureSystem;
        template <int dim>      struct CompositionSystem;
      }

      namespace CopyData
      {
        template <int dim>      struct StokesPreconditioner;
        template <int dim>      struct StokesSystem;
        template <int dim>      struct TemperatureSystem;
        template <int dim>      struct CompositionSystem;
      }
    }
  }

  /**
   * SimulatorAccess is base class for different plugins like postprocessors. It provides access to
   * the various variables of the main class that plugins may want to use
   * in their evaluations, such as solution vectors, the current time, time step
   * sizes, material models, or the triangulations and DoFHandlers that correspond to solutions.
   *
   * This class is the interface between plugins and the main simulator
   * class. Using this insulation layer, the plugins need not know anything
   * about the internal details of the simulation class.
   *
   * Every Postprocessor is required to derive from SimulatorAccess. It is optional
   * for other plugins like MaterialModel, GravityModel, etc..
   *
   * Since the functions providing access to details of the simulator class
   * are meant to be used only by derived classes of this class (rather
   * than becoming part of the public interface of these classes), the
   * functions of this class are made @p protected.
   *
   * @ingroup ?
   */
  template <int dim>
  class SimulatorAccess
  {
    public:
      /**
       * Destructor. Does nothing but is virtual so that derived classes
       * destructors are also virtual.
       **/
      virtual
      ~SimulatorAccess ();

      /**
       * Initialize this class for a given simulator. This function
       * is marked as virtual so that derived classes can do something
       * upon initialization as well, for example look up and cache
       * data; derived classes should call this function from the base
       * class as well, however.
       *
       * @param simulator A reference to the main simulator object.
       **/
      virtual void initialize (const Simulator<dim> &simulator);

    protected:
      /** @name Accessing variables that identify overall properties of the simulator */
      /** @{ */

      /**
       * Return the MPI communicator for this simulation.
       */
      MPI_Comm
      get_mpi_communicator () const;

      /**
       * Return a reference to the stream object that only outputs something on one
      * processor in a parallel program and simply ignores output put into it on
      * all other processors.
       */
      const ConditionalOStream &
      get_pcout () const;

      /**
       * Return the current simulation time in seconds.
       */
      double get_time () const;

      /**
       * Return the size of the last time step.
       */
      double
      get_timestep () const;

      /**
       * Return the current number of a time step.
       */
      unsigned int
      get_timestep_number () const;

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
       * Return a reference to the mapping used to describe the boundary
      * of the domain.
       */
      const Mapping<dim> &
      get_mapping () const;

      /**
       * Return the directory specified in the input parameter file to be
       * the place where output files are to be placed. The string is
       * terminated by a directory separator (i.e., '/').
       */
      std::string
      get_output_directory () const;

      /**
       * Return whether we use the adiabatic heating term.
       */
      bool
      include_adiabatic_heating () const;

      /**
       * Return the adiabatic surface temperature.
       */
      double
      get_adiabatic_surface_temperature () const;

      /**
      * Return whether things like velocities should be converted from
      * the seconds in the MKS system to years. The value of this flag
      * is set by the corresponding entry in the input parameter file.
      */
      bool
      convert_output_to_years () const;

      /**
      * Return the number of compositional fields specified in the input
      * parameter file that will be advected along with the flow field.
      */
      unsigned int
      n_compositional_fields () const;

      /**
      * Compute the error indicators in the same way they are normally used
      * for mesh refinement. The mesh is not refined when doing so, but the
      * indicators can be used when generating graphical output to check
      * why mesh refinement is proceeding as it is.
      */
      void
      get_refinement_criteria(Vector<float> &estimated_error_per_cell) const;
      /** @} */


      /** @name Accessing variables that identify the solution of the problem */
      /** @{ */


      /**
       * Return a reference to the vector that has the current
       * solution of the entire system, i.e. the velocity and
       * pressure variables as well as the temperature.  This vector
       * is associated with the DoFHandler object returned by
       * get_dof_handler().
       *
       * @note In general the vector is a distributed vector; however, it
       * contains ghost elements for all locally relevant degrees of freedom.
       */
      const LinearAlgebra::BlockVector &
      get_solution () const;

      /**
       * Return a reference to the vector that has the solution
       * of the entire system at the previous time step.
       * This vector is associated with the DoFHandler object returned by
       * get_stokes_dof_handler().
       *
       * @note In general the vector is a distributed vector; however, it
       * contains ghost elements for all locally relevant degrees of freedom.
       */
      const LinearAlgebra::BlockVector &
      get_old_solution () const;

      /**
       * Return a reference to the DoFHandler that is used to discretize
       * the variables at the current time step.
       */
      const DoFHandler<dim> &
      get_dof_handler () const;

      /**
       * Fill the argument with a set of depth averages of the current
       * temperature field. The function fills a
       * vector that contains average field values over slices of the
       * domain of same depth. The function resizes the output vector
       * to match the number of depth slices.
       *
       * @param values The output vector of depth averaged values.
       */
      void
      get_depth_average_temperature(std::vector<double> &values) const;

      /**
       * Fill the argument with a set of depth averages of the current
       * compositional fields. See get_depth_average_temperature.
       *
       * @param composition_index The index of the compositional field whose
       * matrix we want to assemble (0 <= composition_index < number of
       * compositional fields in this problem).
       *
       * @param values The output vector of depth averaged values.
       */
      void
      get_depth_average_composition(const unsigned int composition_index, std::vector<double> &values) const;

      /**
       * Compute a depth average of the current viscosity
       */
      void
      get_depth_average_viscosity(std::vector<double> &values) const;

      /**
       * Compute a depth average of the current velocity magnitude
       */
      void
      get_depth_average_velocity_magnitude(std::vector<double> &values) const;

      /**
       * Compute a depth average of the current sinking velocity
       */
      void
      get_depth_average_sinking_velocity(std::vector<double> &values) const;

      /**
       * Compute the seismic shear wave speed, Vs anomaly per element
       */
      void
      get_Vs_anomaly(Vector<float> &values) const;

      /**
       * Compute the seismic pressure wave speed, Vp anomaly per element
       */
      void
      get_Vp_anomaly(Vector<float> &values) const;

      /**
       * Compute a depth average of the seismic shear wave speed: Vs
       */
      void
      get_depth_average_Vs(std::vector<double> &values) const;
      /** @} */

      /**
       * Compute a depth average of the seismic pressure wave speed: Vp
       */
      void
      get_depth_average_Vp(std::vector<double> &values) const;
      /** @} */


      /** @name Accessing variables that identify aspects of the simulation */
      /** @{ */

      /**
       * Return a pointer to the material model to access functions like density().
       */
      const MaterialModel::Interface<dim> &
      get_material_model () const;

      /**
       * Return a pointer to the gravity model description.
       */
      const GravityModel::Interface<dim> &
      get_gravity_model () const;

      /**
       * Return a pointer to the geometry model.
       */
      const GeometryModel::Interface<dim> &
      get_geometry_model () const;


      /**
       * Return a pointer to the object that describes the adiabatic conditions.
       */
      const AdiabaticConditions<dim> &
      get_adiabatic_conditions () const;

      /**
       * Return a pointer to the object that describes the temperature boundary
       * values.
       */
      const BoundaryTemperature::Interface<dim> &
      get_boundary_temperature () const;

      /**
       * Copy the values of the compositional fields at the quadrature point
       * q given as input parameter to the output vector
       * composition_values_at_q_point.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      void
      get_composition_values_at_q_point (const std::vector<std::vector<double>> &composition_values,
                                         const unsigned int                      q,
                                         std::vector<double>                    &composition_values_at_q_point) const;

      /** @} */

    private:
      /**
       * A pointer to the simulator object to which we want to
       * get access.
       */
      SmartPointer<const Simulator<dim> > simulator;
  };

  /**
   * This is the main class of ASPECT. It implements the overall simulation
   * algorithm using the numerical methods discussed in the papers and manuals
   * that accompany ASPECT.
   *
   * @ingroup Simulator
   */
  template <int dim>
  class Simulator : public Subscriptor
  {
    private:
      /**
       * A structure that contains enum values that identify the nonlinear solver in use.
       */
      struct NonlinearSolver
      {
        enum Kind
        {
          IMPES,
          iterated_IMPES,
          iterated_Stokes
        };
      };

    public:
      /**
       * A structure that holds run-time parameters. These parameters are all
       * declared for the ParameterHandler class in the declare_parameters()
       * member function, and read in the parse_parameters() function.
       *
       * Each of the member variables of this class corresponds to a parameter
       * declared for the ParameterHandler class. Rather than duplicating the
       * documentation of each of these parameters for the member variables
       * here, please refer to the documentation of run-time parameters in
       * the ASPECT manual for more information.
       *
       * @ingroup Simulator
       */
      struct Parameters
      {
        /**
         * Constructor. Fills the values of member functions from the
         * given parameter object.
         *
         * @param prm The parameter object that has previously been filled
         * with content by reading an input file.
         **/
        Parameters (ParameterHandler &prm);

        /**
         * Declare the run-time parameters this class takes,
         * and call the respective <code>declare_parameters</code>
         * functions of the namespaces that describe geometries,
         * material models, etc.
         *
         * @param prm The object in which the run-time parameters
         * are to be declared.
         **/
        static
        void declare_parameters (ParameterHandler &prm);

        /**
         * Read run-time parameters from an object that has previously
         * parsed an input file.
         *
         * @param prm The object from which to obtain the run-time parameters.
         **/
        void parse_parameters (ParameterHandler &prm);

        /**
         * @name Global parameters
         * @{
         */
        typedef typename NonlinearSolver::Kind NonlinearSolverKind;

        NonlinearSolverKind            nonlinear_solver;

        bool                           resume_computation;
        double                         start_time;
        double                         end_time;
        double                         CFL_number;
        bool                           convert_to_years;
        std::string                    output_directory;
        double                         surface_pressure;
        double                         adiabatic_surface_temperature;
        unsigned int                   timing_output_frequency;
        double                         linear_solver_tolerance;
        double                         temperature_solver_tolerance;
        double                         composition_solver_tolerance;
        /**
         * @}
         */

        /**
         * @name Parameters that have to do with terms in the model
         * @{
         */
        bool                           include_shear_heating;
        bool                           include_adiabatic_heating;
        double                         radiogenic_heating_rate;
        std::set<types::boundary_id_t> fixed_temperature_boundary_indicators;
        std::set<types::boundary_id_t> zero_velocity_boundary_indicators;
        std::set<types::boundary_id_t> tangential_velocity_boundary_indicators;
        std::map<types::boundary_id_t, std::string> prescribed_velocity_boundary_indicators;
        /**
         * @}
         */

        /**
         * @name Parameters that have to do with mesh refinement
         * @{
         */
        unsigned int                   initial_global_refinement;
        unsigned int                   initial_adaptive_refinement;
        double                         refinement_fraction;
        double                         coarsening_fraction;
        std::string                    refinement_strategy;
        std::vector<double>            additional_refinement_times;
        unsigned int                   adaptive_refinement_interval;
        bool                           run_postprocessors_on_initial_refinement;
        /**
         * @}
         */

        /**
         * @name Parameters that have to do with the stabilization of transport equations
         * @{
         */
        double                         stabilization_alpha;
        double                         stabilization_c_R;
        double                         stabilization_beta;
        /**
         * @}
         */

        /**
         * @name Parameters that have to do with checkpointing
         * @{
         */
        int                            checkpoint_time_secs;
        int                            checkpoint_steps;
        /**
         * @}
         */

        /**
         * @name Parameters that have to do with spatial discretizations
         * @{
         */
        unsigned int                   stokes_velocity_degree;
        bool                           use_locally_conservative_discretization;
        unsigned int                   temperature_degree;
        unsigned int                   composition_degree;
        std::string                    pressure_normalization;
        /**
         * @}
         */

        /**
         * @name Parameters that have to do with compositional field
         * @{
         */
        unsigned int                   n_compositional_fields;
        std::vector<unsigned int>      normalized_fields;
        /**
         * @}
         */
      };

      /**
       * Constructor.
       *
       * @param mpi_communicator The MPI communicator on which this
       * class is to work. The class creates a clone of the actual
       * communicator to make its communications private from the
       * rest of the world.
       *
       * @param prm The run-time parameter object from which this class
       * obtains its settings.
       *
       * This function is implemented in
       * <code>source/simulator/core.cc</code>.
       */
      Simulator (const MPI_Comm mpi_communicator,
                 ParameterHandler &prm);

      /**
       * Destructor. Destroy what needs to be destroyed after waiting
       * for all threads that may still be doing something in the
       * background.
       */
      ~Simulator ();

      /**
       * Declare the run-time parameters this class takes,
       * and call the respective <code>declare_parameters</code>
       * functions of the namespaces that describe geometries,
       * material models, etc.
       *
       * @param prm The object in which the run-time parameters
       * are to be declared.
       *
       * This function is implemented in
       * <code>source/simulator/parameters.cc</code>.
       */
      static
      void declare_parameters (ParameterHandler &prm);

      /**
       * The function that runs the overall algorithm. It
       * contains the loop over all time steps as well as the logic
       * of what to do when before the loop starts and within the time
       * loop.
       *
       * This function is implemented in
       * <code>source/simulator/core.cc</code>.
       */
      void run ();

    private:
      /**
       * @name Top-level functions in the overall flow of the numerical algorithm
       * @{
       */

      /**
       * The function that sets up the DoFHandler objects, It also sets up the
       * various partitioners and computes those constraints on the Stokes
       * variable and temperature that are the same between all time steps.
       *
       * This function is implemented in
       * <code>source/simulator/core.cc</code>.
       */
      void setup_dofs ();

      /**
       * A function that is responsible for initializing the temperature/compositional
       * field before the first time step. The temperature field then serves as the
       * temperature from which the velocity is computed during the first time
       * step, and is subsequently overwritten by the temperature field one gets
       * by advancing by one time step.
       *
       * This function is implemented in
       * <code>source/simulator/initial_conditions.cc</code>.
       */
      void set_initial_temperature_and_compositional_fields ();

      /**
       * A function that initializes the pressure variable before the first
       * time step. It does so by either interpolating (for continuous pressure
       * finite elements) or projecting (for discontinuous elements) the adiabatic
       * pressure computed from the material model.
       *
       * Note that the pressure so set is overwritten by the pressure in fact
       * computed during the first time step. We need this function, however, so
       * that the evaluation of pressure-dependent coefficients (e.g. pressure
       * dependent densities or thermal coefficients) during the first
       * time step has some useful pressure to start with.
       *
       * This function is implemented in
       * <code>source/simulator/initial_conditions.cc</code>.
       */
      void compute_initial_pressure_field ();

      /**
       * Do some housekeeping at the beginning of each time step. This includes
       * generating some screen output, adding some information to the statistics
       * file, and interpolating time-dependent boundary conditions specific to
       * this particular time step (the time independent boundary conditions, for
       * example for hanging nodes or for tangential flow, are computed only
       * once per mesh in setup_dofs()).
       *
       * This function is implemented in
       * <code>source/simulator/core.cc</code>.
       */
      void start_timestep ();

      /**
       * Do the various steps necessary to assemble and solve the things
       * necessary in each time step.
       *
       * This function is implemented in
       * <code>source/simulator/core.cc</code>.
       */
      void solve_timestep ();

      /**
       * Initiate the assembly of the Stokes preconditioner matrix via
       * assemble_stokes_preconditoner(), then set up the data structures
       * to actually build a preconditioner from this matrix.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void build_stokes_preconditioner ();

      /**
       * Initialize the preconditioner for the temperature equation.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void build_temperature_preconditioner ();

      /**
       * Initialize preconditioner for the composition equation.
       *
       * @param composition_index The index of the compositional field whose
       * preconditioner we want to build (0 <= composition_index < number of
       * compositional fields in this problem).
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void build_composition_preconditioner (unsigned int composition_index);

      /**
       * Initiate the assembly of the Stokes matrix and right hand side.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void assemble_stokes_system ();

      /**
       * Initiate the assembly of the temperature matrix and right hand side.
       * This function does not build a preconditioner for the matrix as one
       * may want to re-use a preconditioner initialized using a previously
       * computed matrix.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void assemble_temperature_system ();

      /**
       * Initiate the assembly of the composition matrix and right hand side
       * and build a preconditioner for the matrix.
       *
       * @param composition_index The index of the compositional field whose
       * matrix we want to assemble (0 <= composition_index < number of
       * compositional fields in this problem).
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void assemble_composition_system (unsigned int composition_index);

      /**
       * Solve one block of the the temperature/composition linear system. Return the initial
       * nonlinear residual, i.e., if the linear system to be solved is $Ax=b$, then
       * return $\|Ax_0-b\|$ where $x_0$ is the initial guess for the solution variable
       * and is taken from the current_linearization_point member variable.
       *
       * @param index The index of the block to be solved:
       * 0                              temperature
       * 1...n_compositional_fields     compositional field
       *
       * This function is implemented in
       * <code>source/simulator/solver.cc</code>.
       */
      double solve_temperature_or_composition (unsigned int index);

      /**
       * Solve the Stokes linear system. Return the initial nonlinear residual,
       * i.e., if the linear system to be solved is $Ax=b$, then return $\|Ax_0-b\|$
       * where $x_0$ is the initial guess for the solution variable and is taken from
       * the current_linearization_point member variable.
       *
       * This function is implemented in
       * <code>source/simulator/solver.cc</code>.
       */
      double solve_stokes ();

      /**
       * Handles assembly and solving of the temperature
       * and Stokes system, handling the nonlinear system
       * in different ways.
       */
      void solve_system ();

      /**
       * This function is called at the end of every time step. It
       * runs all the postprocessors that have been listed in the input
       * parameter file (see the manual) in turn. In particular, this
       * usually includes generating graphical output every few time steps.
       *
       * The function also updates the statistics output file at the end of
       * each time step.
       *
       * This function is implemented in
       * <code>source/simulator/core.cc</code>.
       */
      void postprocess ();

      /**
       * Compute error indicators based on a variety of criteria to determine which
       * cells to refine and coarsen, which is done in refine_mesh().
       *
       * This function is implemented in
       * <code>source/simulator/core.cc</code>.
       */
      void compute_refinement_criterion (Vector<float> &estimated_error_per_cell) const;

      /**
       * Refine the mesh according to error indicators calculated by
       * compute_refinement_criterion(), set up all necessary data structures
       * on this new mesh, and interpolate the old solutions onto the new
       * mesh.
       *
       * This function is implemented in
       * <code>source/simulator/core.cc</code>.
       */
      void refine_mesh (const unsigned int max_grid_level);

      /**
       * @}
       */

      /**
       * @name Functions used in saving the state of the program and restarting from a saved state
       * @{
       */
      /**
       * Save the state of this program to a set of files in the output
       * directory. In reality, however, only some variables are stored (in
       * particular the mesh, the solution vectors, etc) whereas others can
       * either be re-generated (matrices, DoFHandler objects, etc) or are
       * read from the input parameter file. See the manual for more
       * information.
       *
       * This function is implemented in
       * <code>source/simulator/checkpoint_restart.cc</code>.
       */
      void create_snapshot();

      /**
       * Restore the state of this program from a set of files in the output
       * directory. In reality, however, only some variables are stored (in
       * particular the mesh, the solution vectors, etc) whereas others can
       * either be re-generated (matrices, DoFHandler objects, etc) or are
       * read from the input parameter file. See the manual for more
       * information. This function only restores those variables that can
       * neither be re-generated from other information nor are read from
       * the input parameter file.
       *
       * This function is implemented in
       * <code>source/simulator/checkpoint_restart.cc</code>.
       */
      void resume_from_snapshot();

      /**
       * Save a number of variables using BOOST serialization mechanism.
       *
       * This function is implemented in
       * <code>source/simulator/checkpoint_restart.cc</code>.
       */
      template <class Archive>
      void serialize (Archive &ar, const unsigned int version);
      /**
       * @}
       */

      /**
       * @name Functions used in setting up linear systems
       * @{
       */
      /**
       * Set up the size and structure of the matrix used to store the
       * elements of the linear system.
       *
       * This function is implemented in
       * <code>source/simulator/core.cc</code>.
       */
      void setup_system_matrix (const std::vector<IndexSet> &system_partitioning);

      /**
       * Set up the size and structure of the matrix used to store the
       * elements of the matrix that is used to build the preconditioner for
       * the system.
       *
       * This function is implemented in
       * <code>source/simulator/core.cc</code>.
       */
      void setup_system_preconditioner (const std::vector<IndexSet> &system_partitioning);

      /**
       * @}
       */

      /**
       * @name Functions used in the assembly of linear systems
       * @{
       */
      /**
       * Initiate the assembly of the preconditioner for the Stokes system.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void assemble_stokes_preconditioner ();

      /**
       * Compute the integrals for the preconditioner for the Stokes system
       * on a single cell.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void
      local_assemble_stokes_preconditioner (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                            internal::Assembly::Scratch::StokesPreconditioner<dim> &scratch,
                                            internal::Assembly::CopyData::StokesPreconditioner<dim> &data);

      /**
       * Copy the contribution to the preconditioner for the Stokes system
       * from a single cell into the global matrix that stores these elements.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void
      copy_local_to_global_stokes_preconditioner (const internal::Assembly::CopyData::StokesPreconditioner<dim> &data);

      /**
       * Compute the integrals for the Stokes matrix and right hand side
       * on a single cell.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void
      local_assemble_stokes_system (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                    internal::Assembly::Scratch::StokesSystem<dim>  &scratch,
                                    internal::Assembly::CopyData::StokesSystem<dim> &data);

      /**
       * Copy the contribution to the Stokes system
       * from a single cell into the global matrix that stores these elements.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void
      copy_local_to_global_stokes_system (const internal::Assembly::CopyData::StokesSystem<dim> &data);

      /**
       * Compute the integrals for the temperature matrix and right hand side
       * on a single cell.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void
      local_assemble_temperature_system (const std::pair<double,double> global_T_range,
                                         const double                   global_max_velocity,
                                         const double                   global_entropy_variation,
                                         const typename DoFHandler<dim>::active_cell_iterator &cell,
                                         internal::Assembly::Scratch::TemperatureSystem<dim>  &scratch,
                                         internal::Assembly::CopyData::TemperatureSystem<dim> &data);

      /**
       * Copy the contribution to the temperature system
       * from a single cell into the global matrix that stores these elements.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void
      copy_local_to_global_temperature_system (const internal::Assembly::CopyData::TemperatureSystem<dim> &data);

      /**
       * Compute the integrals for the composition matrix and right hand side
       * on a single cell.
       *
       * @param composition_index The index of the compositional field whose
       * local matrix we want to assemble (0 <= composition_index < number of
       * compositional fields in this problem).
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void
      local_assemble_composition_system (const unsigned int             composition_index,
                                         const std::pair<double,double> global_C_range,
                                         const double                   global_max_velocity,
                                         const double                   global_entropy_variation,
                                         const typename DoFHandler<dim>::active_cell_iterator &cell,
                                         internal::Assembly::Scratch::CompositionSystem<dim>  &scratch,
                                         internal::Assembly::CopyData::CompositionSystem<dim> &data);

      /**
       * Copy the contribution to the composition system
       * from a single cell into the global matrix that stores these elements.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void
      copy_local_to_global_composition_system (const internal::Assembly::CopyData::CompositionSystem<dim> &data);

      /**
       * @}
       */

      /**
       * @name Helper functions
       * @{
       */
      /**
       * This routine adjusts the second block of the right hand side of a
       * Stokes system  (containing the term that comes from compressibility,
       * so that the system becomes
       * compatible: $0=\int div u = \int g$. The vector to adjust is given as the
       * argument of this function. This function makes use of the
       * helper vector pressure_shape_function_integrals that contains
       * $h_i=(q_i,1)$ with the pressure functions $q_i$
       * and we adjust the right hand side $g$ by $h_i \int g / |\Omega|$.
       *
       * The purpose of this function is described in the second paper on
       * the numerical methods in Aspect.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      void make_pressure_rhs_compatible(LinearAlgebra::BlockVector &vector);

      /**
       * Computes a running average of a vector.
       * http://en.wikipedia.org/wiki/Moving_average
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      void compute_running_average(std::vector<double> &values, const int npoints) const;

      /**
       * Compute a depth average of the current temperature/composition. The function
       * fills a vector that contains average temperatures/compositions over slices of the
       * domain of same depth. The function resizes the output vector to match
       * the number of depth slices.
       *
       * @param index The index of the block to be solved:
       * 0                              temperature
       * 1...n_compositional_fields     compositional field
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      void compute_depth_average_field(const unsigned int index,
				       std::vector<double> &values) const;

      /**
       * Compute a depth average of the current temperature. The function
       * fills a vector that contains average temperatures over slices of the
       * domain of same depth. The function resizes the output vector to match
       * the number of depth slices.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      void compute_depth_average_viscosity(std::vector<double> &values) const;

      /**
       * Compute a depth average of the current velocity magnitude.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      void compute_depth_average_velocity_magnitude(std::vector<double> &values) const;

      /**
       * Compute a depth average of the current sinking velocity.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      void compute_depth_average_sinking_velocity(std::vector<double> &values) const;

      /**
       * Compute a depth average of the seismic shear wave speed, Vs
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      void compute_depth_average_Vs(std::vector<double> &values) const;

      /**
       * Compute a depth average of the seismic pressure wave speed, Vp
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      void compute_depth_average_Vp(std::vector<double> &values) const;

      /**
       * Compute the seismic shear wave speed, Vs anomaly per element.
       * we compute the anomaly by computing a smoothed (over 200 km or so) laterally averaged
       * temperature profile and associated seismic velocity that is then subtracted from the
       * seismic velocity at the current pressure temperature conditions
      */
      void compute_Vs_anomaly(Vector<float> &values) const;

      /**
       * Compute the seismic pressure wave speed, Vp anomaly per element.
       * we compute the anomaly by computing a smoothed (over 200 km or so) laterally averaged
       * temperature profile and associated seismic velocity that is then subtracted from the
       * seismic velocity at the current pressure temperature conditions
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      void compute_Vp_anomaly(Vector<float> &values) const;

      /**
       * Adjust the pressure variable
       * (which is only determined up to a constant) by adding a constant to it
       * in such a way that the pressure on the surface has a known average value.
       * Whether a face is part of the surface is determined by asking whether its
       * depth of its midpoint (as determined by the geometry model) is less than
       * 1/3*1/sqrt(dim-1)*diameter of the face. For reasonably curved boundaries,
       * this rules out side faces that are perpendicular ot the surface boundary
       * but includes those faces that are along the boundary even if the real
       * boundary is curved.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      void normalize_pressure(LinearAlgebra::BlockVector &vector);

      /**
       * Compute the maximal velocity throughout the domain. This is needed
       * to compute the size of the time step.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      double get_maximal_velocity (const LinearAlgebra::BlockVector &solution) const;

      /**
       * Compute the variation (i.e., the difference between maximal and
       * minimal value) of the entropy $(T-\bar T)^2$ where $\bar T$ is the
       * average temperature throughout the domain given as argument to
       * this function.
       *
       * This function is used in computing the artificial diffusion
       * stabilization term.
       *
       * @param index The index of the field we want to calculate the entropy
       * variation of:
       * 0                              temperature
       * 1...n_compositional_fields     compositional field
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      double get_entropy_variation (const double average_temperature, const unsigned int index) const;

      /**
       * Compute the minimal and maximal temperature througout the domain from a
       * solution vector extrapolated from the previous time steps. This is needed
       * to compute the artificial diffusion stabilization terms.
       *
       * @param index The index of the field we want to calculate the entropy
       * variation of:
       * 0                              temperature
       * 1...n_compositional_fields     compositional field
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      std::pair<double,double> get_extrapolated_temperature_or_composition_range (const unsigned int index) const;

      /**
       * Compute the size of the next time step from the mesh size and
       * the velocity on each cell. The computed time step has to satisfy
       * the CFL number chosen in the input parameter file on each cell
       * of the mesh.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      double compute_time_step () const;

      /**
       * Compute the artificial diffusion coefficient value on a cell
       * given the values and gradients of the solution passed as
       * arguments.
       *
       * @param index The index of the field we need the artificial
       * diffusion coefficient for:
       * 0                              temperature
       * 1...n_compositional_fields     compositional field
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      double
      compute_viscosity(const std::vector<double>          &old_temperature,
                        const std::vector<double>          &old_old_temperature,
                        const std::vector<Tensor<1,dim> >  &old_temperature_grads,
                        const std::vector<Tensor<1,dim> >  &old_old_temperature_grads,
                        const std::vector<double>          &old_temperature_laplacians,
                        const std::vector<double>          &old_old_temperature_laplacians,
                        const std::vector<Tensor<1,dim> >  &old_velocity_values,
                        const std::vector<Tensor<1,dim> >  &old_old_velocity_values,
                        const std::vector<SymmetricTensor<2,dim> >  &old_strain_rates,
                        const std::vector<SymmetricTensor<2,dim> >  &old_old_strain_rates,
                        const std::vector<double>          &old_pressure,
                        const std::vector<double>          &old_old_pressure,
                        const std::vector<std::vector<double>> &old_composition,
                        const std::vector<std::vector<double>> &old_old_composition,
                        const double                        global_u_infty,
                        const double                        global_T_variation,
                        const double                        average_temperature,
                        const double                        global_entropy_variation,
                        const std::vector<Point<dim> >     &evaluation_points,
                        const double                        cell_diameter,
                        const unsigned int                  index) const;

      /**
       * Compute the residual of the temperature equation to be used
       * for the artificial diffusion coefficient value on a cell
       * given the values and gradients of the solution
       * passed as arguments.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void
      compute_temperature_system_residual(const std::vector<double>          &old_temperature,
                                          const std::vector<double>          &old_old_temperature,
                                          const std::vector<Tensor<1,dim> >  &old_temperature_grads,
                                          const std::vector<Tensor<1,dim> >  &old_old_temperature_grads,
                                          const std::vector<double>          &old_temperature_laplacians,
                                          const std::vector<double>          &old_old_temperature_laplacians,
                                          const std::vector<Tensor<1,dim> >  &old_velocity_values,
                                          const std::vector<Tensor<1,dim> >  &old_old_velocity_values,
                                          const std::vector<SymmetricTensor<2,dim> >  &old_strain_rates,
                                          const std::vector<SymmetricTensor<2,dim> >  &old_old_strain_rates,
                                          const std::vector<double>          &old_pressure,
                                          const std::vector<double>          &old_old_pressure,
                                          const std::vector<std::vector<double>> &old_composition,
                                          const std::vector<std::vector<double>> &old_old_composition,
                                          const double                        average_temperature,
                                          const std::vector<Point<dim> >     &evaluation_points,
                                          double                             &max_residual,
                                          double                             &max_velocity) const;

      /**
       * Compute the residual of the composition equation to be used
       * for the artificial diffusion coefficient value on a cell
       * given the values and gradients of the solution
       * passed as arguments.
       *
       * @param composition_index The index of the compositional field whose
       * residual we want to get (0 <= composition_index < number of
       * compositional fields in this problem).
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void
      compute_composition_system_residual(const std::vector<double>          &old_composition,
                                          const std::vector<double>          &old_old_composition,
                                          const std::vector<Tensor<1,dim> >  &old_composition_grads,
                                          const std::vector<Tensor<1,dim> >  &old_old_composition_grads,
                                          const std::vector<double>          &old_composition_laplacians,
                                          const std::vector<double>          &old_old_composition_laplacians,
                                          const std::vector<Tensor<1,dim> >  &old_velocity_values,
                                          const std::vector<Tensor<1,dim> >  &old_old_velocity_values,
                                          const double                        average_composition,
                                          const unsigned int                  composition_index,
                                          double                             &max_residual,
                                          double                             &max_velocity) const;
      /**
       * Extract the values of temperature, pressure, composition and
       * optional strain rate for the current linearization point.
       * These values are stored as input arguments for the material
       * model. The compositional fields are extracted with the
       * individual compositional fields as outer vectors and the values
       * at each quadrature point as inner vectors, but the material
       * model needs it the other way round. Hence, this vector of vectors
       * is transposed.
       *
       * @param compute_strainrate If the strain rate should be computed
       * or not
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void
      compute_material_model_input_values (const TrilinosWrappers::MPI::BlockVector                    &current_linearization_point,
                                           const FEValues<dim,dim>                                     &input_finite_element_values,
                                           const bool                                                   compute_strainrate,
                                           typename MaterialModel::Interface<dim>::MaterialModelInputs &material_model_inputs) const;

      /**
       * Copy the values of the compositional fields at the quadrature point
       * q given as input parameter to the output vector
       * composition_values_at_q_point.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      void
      extract_composition_values_at_q_point (const std::vector<std::vector<double>> &composition_values,
                                             const unsigned int                      q,
                                             std::vector<double>                    &composition_values_at_q_point) const;

      /**
       * Generate and output some statistics like timing information and memory consumption.
       * Whether this function does anything or not is controlled through the
       * variable aspect::output_parallel_statistics.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      void output_program_stats();

      /**
       * This function is called at the end of each time step and writes the
       * statistics object that contains data like the current time, the number
       * of linear solver iterations, and whatever the postprocessors have
       * generated, to disk.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      void output_statistics();
      /**
       * @}
       */

      /**
       * @name Variables that have to do with input, output and parallel communication
       * @{
       */
      Parameters                          parameters;

      MPI_Comm                            mpi_communicator;

      ConditionalOStream                  pcout;

      /**
       * An object that stores a bunch of statistics such as the number of
       * linear solver iterations, the time corresponding to each time
       * step, etc, as well as whatever the various postprocessors want
       * to put into it.
       *
       * This variable is written to disk after every time step, by the
       * Simulator::output_statistics() function.
       */
      TableHandler                        statistics;

      Postprocess::Manager<dim>           postprocess_manager;
      TimerOutput                         computing_timer;

      /**
       * In output_statistics(), where we output the statistics object above,
       * we do the actual writing on a separate thread. This variable is
       * the handle we get for this thread so that we can wait for it
       * to finish, either if we want to write the statistics object for
       * the next thread, or if we want to terminate altogether.
       */
      Threads::Thread<>                   output_statistics_thread;
      /**
       * @}
       */

      /**
       * @name Variables that describe the physical setup of the problem
       * @{
       */
      const std::auto_ptr<GeometryModel::Interface<dim> >       geometry_model;
      const std::auto_ptr<MaterialModel::Interface<dim> >       material_model;
      const std::auto_ptr<GravityModel::Interface<dim> >        gravity_model;
      const std::auto_ptr<BoundaryTemperature::Interface<dim> > boundary_temperature;
      std::auto_ptr<    InitialConditions::Interface<dim> >   initial_conditions;
      std::auto_ptr<const CompositionalInitialConditions::Interface<dim> >   compositional_initial_conditions;
      std::auto_ptr<const AdiabaticConditions<dim> >            adiabatic_conditions;
      std::map<types::boundary_id_t,std_cxx1x::shared_ptr<VelocityBoundaryConditions::Interface<dim> > > velocity_boundary_conditions;
      /**
       * @}
       */

      /**
       * @name Variables that describe the time discretization
       * @{
       */
      double                                                    time;
      double                                                    time_step;
      double                                                    old_time_step;
      unsigned int                                              timestep_number;
      /**
       * @}
       */

      /**
       * @name Variables related to checkpointing
       * @{
       */
      time_t                                                    last_checkpoint_time;
      /**
       * @}
       */

      /**
       * @name Variables that describe the spatial discretization
       * @{
       */
      parallel::distributed::Triangulation<dim>                 triangulation;
      double                                                    global_Omega_diameter;
      double                                                    global_volume;

      const MappingQ<dim>                                       mapping;

      const FESystem<dim>                                       finite_element;

      DoFHandler<dim>                                           dof_handler;

      ConstraintMatrix                                          constraints;
      ConstraintMatrix                                          current_constraints;

      double                                                    pressure_scaling;

      /**
       * @}
       */


      /**
       * @name Variables that describe the linear systems and solution vectors
       * @{
       */
      LinearAlgebra::BlockSparseMatrix                          system_matrix;
      LinearAlgebra::BlockSparseMatrix                          system_preconditioner_matrix;

      LinearAlgebra::BlockVector                                solution;
      LinearAlgebra::BlockVector                                old_solution;
      LinearAlgebra::BlockVector                                old_old_solution;
      LinearAlgebra::BlockVector                                system_rhs;

      TrilinosWrappers::MPI::BlockVector                        current_linearization_point;

      // only used if is_compressible()
      LinearAlgebra::BlockVector                                pressure_shape_function_integrals;



      std_cxx1x::shared_ptr<LinearAlgebra::PreconditionAMG>     Amg_preconditioner;
      std_cxx1x::shared_ptr<LinearAlgebra::PreconditionILU>     Mp_preconditioner;
      std_cxx1x::shared_ptr<LinearAlgebra::PreconditionILU>     T_preconditioner;
//TODO: use n_compositional_field separate preconditioners
      std_cxx1x::shared_ptr<LinearAlgebra::PreconditionILU>     C_preconditioner;

      bool                                                      rebuild_stokes_matrix;
      bool                                                      rebuild_stokes_preconditioner;
      /**
       * @}
       */

      friend class boost::serialization::access;
      friend class SimulatorAccess<dim>;
  };
}


#endif
