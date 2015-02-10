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
#include <deal.II/base/tensor_function.h>

#include <aspect/global.h>
#include <aspect/simulator_access.h>
#include <aspect/material_model/interface.h>
#include <aspect/heating_model/interface.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/boundary_temperature/interface.h>
#include <aspect/boundary_composition/interface.h>
#include <aspect/initial_conditions/interface.h>
#include <aspect/compositional_initial_conditions/interface.h>
#include <aspect/velocity_boundary_conditions/interface.h>
#include <aspect/mesh_refinement/interface.h>
#include <aspect/termination_criteria/interface.h>
#include <aspect/postprocess/interface.h>
#include <aspect/adiabatic_conditions/interface.h>

#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>


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
        template <int dim>      struct AdvectionSystem;
      }

      namespace CopyData
      {
        template <int dim>      struct StokesPreconditioner;
        template <int dim>      struct StokesSystem;
        template <int dim>      struct AdvectionSystem;

      }
    }
  }

  /**
   * This is the main class of ASPECT. It implements the overall simulation
   * algorithm using the numerical methods discussed in the papers and manuals
   * that accompany ASPECT.
   *
   * @ingroup Simulator
   */
  template <int dim>
  class Simulator
  {
    private:
      /**
       * A structure that contains enum values that identify the nonlinear
       * solver in use.
       */
      struct NonlinearSolver
      {
        enum Kind
        {
          IMPES,
          iterated_IMPES,
          iterated_Stokes,
          Stokes_only
        };
      };

      struct NullspaceRemoval
      {
        enum Kind
        {
          none = 0,
          net_translation_x = 0x1,
          net_translation_y = 0x2,
          net_translation_z = 0x4,
          net_translation   = 0x1+0x2+0x4,
          linear_momentum_x = 0x8,
          linear_momentum_y = 0x10,
          linear_momentum_z = 0x20,
          linear_momentum   = 0x8+0x10+0x20,
          net_rotation      = 0x40,
          angular_momentum  = 0x80
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
       * here, please refer to the documentation of run-time parameters in the
       * ASPECT manual for more information.
       *
       * @ingroup Simulator
       */
      struct Parameters
      {
        /**
         * Constructor. Fills the values of member functions from the given
         * parameter object.
         *
         * @param prm The parameter object that has previously been filled
         * with content by reading an input file.
         *
         * @param mpi_communicator The MPI communicator we will use for this
         * simulation. We need this when calling parse_parameters() so that we
         * can verify some of the input arguments.
         */
        Parameters (ParameterHandler &prm,
                    MPI_Comm mpi_communicator);

        /**
         * Declare the run-time parameters this class takes, and call the
         * respective <code>declare_parameters</code> functions of the
         * namespaces that describe geometries, material models, etc.
         *
         * @param prm The object in which the run-time parameters are to be
         * declared.
         */
        static
        void declare_parameters (ParameterHandler &prm);

        /**
         * Read run-time parameters from an object that has previously parsed
         * an input file. This reads all parameters that do not require
         * knowledge of the geometry model we use. There is a separate
         * function parse_geometry_dependent_parameters() that is called as
         * soon as the geometry object has been created and that can translate
         * between the symbolic names for boundary components that the
         * geometry model publishes and the boundary indicators used
         * internally.
         *
         * @param prm The object from which to obtain the run-time parameters.
         *
         * @param mpi_communicator The MPI communicator we will use for this
         * simulation. We need this when calling parse_parameters() so that we
         * can verify some of the input arguments.
         */
        void parse_parameters (ParameterHandler &prm,
                               const MPI_Comm mpi_communicator);

        /**
         * Read those run-time parameters from a ParameterHandler object that
         * depend on knowing which geometry object we use. This function
         * complements parse_parameters() but is only called once the geometry
         * object has been created. This function is separate because we allow
         * the use of symbolic names in defining which boundary components
         * have which boundary conditions, and the names one can specify there
         * are not available until after the geometry object has been created.
         *
         * This function is called from the GeometryModel::create_geometry()
         * function.
         *
         * @param prm The object from which to obtain the run-time parameters.
         * @param geometry_model The geometry model that provides boundary
         * names etc.
         */
        void parse_geometry_dependent_parameters (ParameterHandler &prm,
                                                  const GeometryModel::Interface<dim> &geometry_model);

        /**
         * @name Global parameters
         * @{
         */
        typedef typename NonlinearSolver::Kind NonlinearSolverKind;

        NonlinearSolverKind            nonlinear_solver;

        double                         nonlinear_tolerance;
        bool                           resume_computation;
        double                         start_time;
        double                         CFL_number;
        double                         maximum_time_step;
        bool                           use_conduction_timestep;
        bool                           convert_to_years;
        std::string                    output_directory;
        double                         surface_pressure;
        double                         adiabatic_surface_temperature;
        unsigned int                   timing_output_frequency;
        bool                           use_direct_stokes_solver;
        double                         linear_stokes_solver_tolerance;
        unsigned int                   max_nonlinear_iterations;
        unsigned int                   n_cheap_stokes_solver_steps;
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
        bool                           include_latent_heat;
        double                         radiogenic_heating_rate;
        std::set<types::boundary_id> fixed_temperature_boundary_indicators;
        std::set<types::boundary_id> fixed_composition_boundary_indicators;
        std::set<types::boundary_id> zero_velocity_boundary_indicators;
        std::set<types::boundary_id> tangential_velocity_boundary_indicators;
        /**
         * map from boundary id to a pair "components", "velocity boundary
         * type", where components is of the format "[x][y][z]" and the
         * velocity type is mapped to one of the plugins of velocity boundary
         * conditions (e.g. "function")
         */
        std::map<types::boundary_id, std::pair<std::string,std::string> > prescribed_velocity_boundary_indicators;
        /**
         * Selection of operations to perform to remove nullspace from
         * velocity field.
         */
        typename NullspaceRemoval::Kind nullspace_removal;
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
        unsigned int                   min_grid_level;
        std::vector<double>            additional_refinement_times;
        unsigned int                   adaptive_refinement_interval;
        bool                           run_postprocessors_on_initial_refinement;
        /**
         * @}
         */

        /**
         * @name Parameters that have to do with the stabilization of
         * transport equations
         * @{
         */
        unsigned int                   stabilization_alpha;
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
         * @name Parameters that have to do with compositional fields
         * @{
         */
        unsigned int                   n_compositional_fields;
        std::vector<std::string>       names_of_compositional_fields;
        std::vector<unsigned int>      normalized_fields;
        /**
         * @}
         */
        /**
         * @name Parameters that have to do with free surface
         * @{
         */
        bool                           free_surface_enabled;
        std::set<types::boundary_id> free_surface_boundary_indicators;
        /**
         * @}
         */

      };

      /**
       * Constructor.
       *
       * @param mpi_communicator The MPI communicator on which this class is
       * to work. The class creates a clone of the actual communicator to make
       * its communications private from the rest of the world.
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
       * Destructor. Destroy what needs to be destroyed after waiting for all
       * threads that may still be doing something in the background.
       */
      ~Simulator ();

      /**
       * Declare the run-time parameters this class takes, and call the
       * respective <code>declare_parameters</code> functions of the
       * namespaces that describe geometries, material models, etc.
       *
       * @param prm The object in which the run-time parameters are to be
       * declared.
       *
       * This function is implemented in
       * <code>source/simulator/parameters.cc</code>.
       */
      static
      void declare_parameters (ParameterHandler &prm);

      /**
       * The function that runs the overall algorithm. It contains the loop
       * over all time steps as well as the logic of what to do when before
       * the loop starts and within the time loop.
       *
       * This function is implemented in
       * <code>source/simulator/core.cc</code>.
       */
      void run ();

    private:

      /**
       * A structure that is used as an argument to functions that can work on
       * both the temperature and the compositional variables and that need to
       * be told which one of the two, as well as on which of the
       * compositional variables.
       */
      struct AdvectionField
      {
        /**
         * An enum indicating whether the identified variable is the
         * temperature or one of the compositional fields.
         */
        enum FieldType { temperature_field, compositional_field };

        /**
         * A variable indicating whether the identified variable is the
         * temperature or one of the compositional fields.
         */
        const FieldType    field_type;

        /**
         * A variable identifying which of the compositional fields is
         * selected. This variable is meaningless if the temperature is
         * selected.
         */
        const unsigned int compositional_variable;

        /**
         * Constructor.
         * @param field_type Determines whether this variable should select
         * the temperature field or a compositional field.
         * @param compositional_variable The number of the compositional field
         * if the first argument in fact chooses a compositional variable.
         * Meaningless if the first argument equals temperature.
         *
         * This function is implemented in
         * <code>source/simulator/helper_functions.cc</code>.
         */
        AdvectionField (const FieldType field_type,
                        const unsigned int compositional_variable = numbers::invalid_unsigned_int);

        /**
         * A static function that creates an object identifying the
         * temperature.
         *
         * This function is implemented in
         * <code>source/simulator/helper_functions.cc</code>.
         */
        static
        AdvectionField temperature ();

        /**
         * A static function that creates an object identifying given
         * compositional field.
         *
         * This function is implemented in
         * <code>source/simulator/helper_functions.cc</code>.
         */
        static
        AdvectionField composition (const unsigned int compositional_variable);

        /**
         * Return whether this object refers to the temperature field.
         */
        bool
        is_temperature () const;

        /**
         * Look up the component index for this temperature or compositional
         * field. See Introspection::component_indices for more information.
         */
        unsigned int component_index(const Introspection<dim> &introspection) const;

        /**
         * Look up the block index for this temperature or compositional
         * field. See Introspection::block_indices for more information.
         */
        unsigned int block_index(const Introspection<dim> &introspection) const;

        /**
         * Look up the base element within the larger composite finite element
         * we used for everything, for this temperature or compositional field
         * See Introspection::base_elements for more information.
         */
        unsigned int base_element(const Introspection<dim> &introspection) const;
      };


      /**
       * A class that is empty but that can be used as a member variable and
       * whose constructor will be run in the order in which the member
       * variables are initialized. Because this class has a constructor that
       * takes a function object that it will execute whenever the member
       * variable is initialized, this allows running arbitrary actions in
       * between member variable initializers, for example if some member
       * variable is partially initialized at point A within the member
       * variable initializer list, its initialization can only be finalized
       * after point B (because it depends on what another member variable
       * decides to do), but needs to be finished by point C within the member
       * initialization. In such a case, one may have a member variable of the
       * current time placed in the list of member variables such that it is
       * initialized at point B, and then initialize it using a function
       * object that performs the finalization of initialization.
       */
      struct IntermediaryConstructorAction
      {
        IntermediaryConstructorAction (std_cxx1x::function<void ()> action);
      };

      /**
       * @name Top-level functions in the overall flow of the numerical
       * algorithm
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
       * This function initializes the variables of the introspection object.
       * It is called by setup_dofs() right after distributing degrees of
       * freedom since this is when all of the information is available.
       *
       * This function is implemented in
       * <code>source/simulator/core.cc</code>.
       */
      void setup_introspection ();

      /**
       * A function that is responsible for initializing the
       * temperature/compositional field before the first time step. The
       * temperature field then serves as the temperature from which the
       * velocity is computed during the first time step, and is subsequently
       * overwritten by the temperature field one gets by advancing by one
       * time step.
       *
       * This function is implemented in
       * <code>source/simulator/initial_conditions.cc</code>.
       */
      void set_initial_temperature_and_compositional_fields ();

      /**
       * A function that initializes the pressure variable before the first
       * time step. It does so by either interpolating (for continuous
       * pressure finite elements) or projecting (for discontinuous elements)
       * the adiabatic pressure computed from the material model.
       *
       * Note that the pressure so set is overwritten by the pressure in fact
       * computed during the first time step. We need this function, however,
       * so that the evaluation of pressure-dependent coefficients (e.g.
       * pressure dependent densities or thermal coefficients) during the
       * first time step has some useful pressure to start with.
       *
       * This function is implemented in
       * <code>source/simulator/initial_conditions.cc</code>.
       */
      void compute_initial_pressure_field ();

      /**
       * Given the 'constraints' member that contains all constraints that are
       * independent of the time (e.g., hanging node constraints, tangential
       * flow constraints, etc), copy it over to 'current_constraints' and add
       * to the latter all constraints that do depend on time such as
       * temperature or velocity Dirichlet boundary conditions. This function
       * is therefore called at the beginning of every time step in
       * start_timestep(), but also when setting up the initial values.
       *
       * This function is implemented in
       * <code>source/simulator/core.cc</code>.
       */
      void compute_current_constraints ();

      /**
       * Do some housekeeping at the beginning of each time step. This
       * includes generating some screen output, adding some information to
       * the statistics file, and interpolating time-dependent boundary
       * conditions specific to this particular time step (the time
       * independent boundary conditions, for example for hanging nodes or for
       * tangential flow, are computed only once per mesh in setup_dofs()).
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
       * assemble_stokes_preconditoner(), then set up the data structures to
       * actually build a preconditioner from this matrix.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void build_stokes_preconditioner ();

      /**
       * Initialize the preconditioner for the advection equation of field
       * index.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void build_advection_preconditioner (const AdvectionField &advection_field,
                                           std_cxx1x::shared_ptr<aspect::LinearAlgebra::PreconditionILU> &preconditioner);

      /**
       * Initiate the assembly of the Stokes matrix and right hand side.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void assemble_stokes_system ();

      /**
       * Initiate the assembly of one advection matrix and right hand side and
       * build a preconditioner for the matrix.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void assemble_advection_system (const AdvectionField &advection_field);

      /**
       * Solve one block of the the temperature/composition linear system.
       * Return the initial nonlinear residual, i.e., if the linear system to
       * be solved is $Ax=b$, then return $\|Ax_0-b\|$ where $x_0$ is the
       * initial guess for the solution variable and is taken from the
       * current_linearization_point member variable.
       *
       * This function is implemented in
       * <code>source/simulator/solver.cc</code>.
       */
      double solve_advection (const AdvectionField &advection_field);

      /**
       * Solve the Stokes linear system. Return the initial nonlinear
       * residual, i.e., if the linear system to be solved is $Ax=b$, then
       * return $\|Ax_0-b\|$ where $x_0$ is the initial guess for the solution
       * variable and is taken from the current_linearization_point member
       * variable. For the purpose of this function, this residual is computed
       * only the velocity and pressure equations (i.e., for the 2x2 block
       * system involving the velocity and pressure variables).
       *
       * This function is implemented in
       * <code>source/simulator/solver.cc</code>.
       */
      double solve_stokes ();

      /**
       * This function is called at the end of every time step. It runs all
       * the postprocessors that have been listed in the input parameter file
       * (see the manual) in turn. In particular, this usually includes
       * generating graphical output every few time steps.
       *
       * The function also updates the statistics output file at the end of
       * each time step.
       *
       * This function is implemented in
       * <code>source/simulator/core.cc</code>.
       */
      void postprocess ();

      /**
       * Refine the mesh according to error indicators calculated by
       * compute_refinement_criterion(), set up all necessary data structures
       * on this new mesh, and interpolate the old solutions onto the new
       * mesh.
       *
       * @param[in] max_grid_level The maximum refinement level of the mesh.
       * This is the sum of the initial global refinement and the initial
       * adaptive refinement (as provided by the user in the input file) and
       * in addition it gets increased by one at each additional refinement
       * time.
       *
       * This function is implemented in
       * <code>source/simulator/core.cc</code>.
       */
      void refine_mesh (const unsigned int max_grid_level);

      /**
       * @}
       */

      /**
       * @name Functions used in saving the state of the program and
       * restarting from a saved state
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
       * neither be re-generated from other information nor are read from the
       * input parameter file.
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
       * Compute the integrals for the preconditioner for the Stokes system on
       * a single cell.
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
       * Compute the integrals for the Stokes matrix and right hand side on a
       * single cell.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void
      local_assemble_stokes_system (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                    internal::Assembly::Scratch::StokesSystem<dim>  &scratch,
                                    internal::Assembly::CopyData::StokesSystem<dim> &data);

      /**
       * Copy the contribution to the Stokes system from a single cell into
       * the global matrix that stores these elements.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void
      copy_local_to_global_stokes_system (const internal::Assembly::CopyData::StokesSystem<dim> &data);

      /**
       * Compute the integrals for one advection matrix and right hand side on
       * a single cell.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void
      local_assemble_advection_system (const AdvectionField &advection_field,
                                       const std::pair<double,double> global_field_range,
                                       const double                   global_max_velocity,
                                       const double                   global_entropy_variation,
                                       const typename DoFHandler<dim>::active_cell_iterator &cell,
                                       internal::Assembly::Scratch::AdvectionSystem<dim>  &scratch,
                                       internal::Assembly::CopyData::AdvectionSystem<dim> &data);

      /**
       * Compute the heating term for the advection system index. Currently
       * the heating term is 0 for compositional fields, but this can be
       * changed in the future to allow for interactions between compositional
       * fields.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      double compute_heating_term(const internal::Assembly::Scratch::AdvectionSystem<dim>  &scratch,
                                  typename MaterialModel::Interface<dim>::MaterialModelInputs &material_model_inputs,
                                  typename MaterialModel::Interface<dim>::MaterialModelOutputs &material_model_outputs,
                                  const AdvectionField &advection_field,
                                  const unsigned int q) const;


      /**
       * Copy the contribution to the advection system from a single cell into
       * the global matrix that stores these elements.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void
      copy_local_to_global_advection_system (const internal::Assembly::CopyData::AdvectionSystem<dim> &data);

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
       * so that the system becomes compatible: $0=\int div u = \int g$. The
       * vector to adjust is given as the argument of this function. This
       * function makes use of the helper vector
       * pressure_shape_function_integrals that contains $h_i=(q_i,1)$ with
       * the pressure functions $q_i$ and we adjust the right hand side $g$ by
       * $h_i \int g / |\Omega|$.
       *
       * The purpose of this function is described in the second paper on the
       * numerical methods in Aspect.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      void make_pressure_rhs_compatible(LinearAlgebra::BlockVector &vector);

      /**
       * Fills a vector with the artificial viscosity for the temperature on
       * each local cell.
       */
      void get_artificial_viscosity (Vector<float> &viscosity_per_cell) const;

      /**
       * Internal routine to compute the depth average of a certain quantitiy.
       *
       * The functor @p fctr must be an object of a user defined type that can
       * be arbitrary but has to satisfy certain requirements. In essence,
       * this class type needs to implement the following interface of member
       * functions:
       * @code
       * template <int dim>
       * class Functor
       * {
       *   public:
       *     // operator() will have @p in and @p out filled out if @p true
       *     bool need_material_properties() const;
       *
       *     // called once at the beginning with the number of quadrature points
       *     void setup(const unsigned int q_points);
       *
       *     // fill @p output for each quadrature point
       *     void operator()(const typename MaterialModel::Interface<dim>::MaterialModelInputs &in,
       *        const typename MaterialModel::Interface<dim>::MaterialModelOutputs &out,
       *        FEValues<dim> &fe_values,
       *        const LinearAlgebra::BlockVector &solution,
       *        std::vector<double> &output);
       * };
       * @endcode
       *
       * @param values The output vector of depth averaged values. The
       * function takes the pre-existing size of this vector as the number of
       * depth slices.
       * @param fctr Instance of a class satisfying the signature above.
       */
      template<class FUNCTOR>
      void compute_depth_average(std::vector<double> &values,
                                 FUNCTOR &fctr) const;

      /**
       * Compute a depth average of the current temperature/composition. The
       * function fills a vector that contains average
       * temperatures/compositions over slices of the domain of same depth.
       * The function resizes the output vector to match the number of depth
       * slices.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       *
       * @param advection_field Temperature or compositional field to average.
       * @param values The output vector of depth averaged values. The
       * function takes the pre-existing size of this vector as the number of
       * depth slices.
       */
      void compute_depth_average_field(const AdvectionField &advection_field,
                                       std::vector<double> &values) const;

      /**
       * Compute a depth average of the current viscosity. The function fills
       * a vector that contains average viscosities over slices of the domain
       * of same depth. The function resizes the output vector to match the
       * number of depth slices.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       *
       * @param values The output vector of depth averaged values. The
       * function takes the pre-existing size of this vector as the number of
       * depth slices.
       */
      void compute_depth_average_viscosity(std::vector<double> &values) const;

      /**
       * Compute a depth average of the current velocity magnitude.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       *
       * @param values The output vector of depth averaged values. The
       * function takes the pre-existing size of this vector as the number of
       * depth slices.
       */
      void compute_depth_average_velocity_magnitude(std::vector<double> &values) const;

      /**
       * Compute a depth average of the current sinking velocity.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       *
       * @param values The output vector of depth averaged values. The
       * function takes the pre-existing size of this vector as the number of
       * depth slices.
       */
      void compute_depth_average_sinking_velocity(std::vector<double> &values) const;

      /**
       * Compute a depth average of the seismic shear wave speed, Vs
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       *
       * @param values The output vector of depth averaged values. The
       * function takes the pre-existing size of this vector as the number of
       * depth slices.
       */
      void compute_depth_average_Vs(std::vector<double> &values) const;

      /**
       * Compute a depth average of the seismic pressure wave speed, Vp
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       *
       * @param values The output vector of depth averaged values. The
       * function takes the pre-existing size of this vector as the number of
       * depth slices.
       */
      void compute_depth_average_Vp(std::vector<double> &values) const;

      /**
       * Compute the seismic shear wave speed, Vs anomaly per element. we
       * compute the anomaly by computing a smoothed (over 200 km or so)
       * laterally averaged temperature profile and associated seismic
       * velocity that is then subtracted from the seismic velocity at the
       * current pressure temperature conditions
       *
       * @param values The output vector of depth averaged values. The
       * function takes the pre-existing size of this vector as the number of
       * depth slices.
       */
      void compute_Vs_anomaly(Vector<float> &values) const;

      /**
       * Compute the seismic pressure wave speed, Vp anomaly per element. we
       * compute the anomaly by computing a smoothed (over 200 km or so)
       * laterally averaged temperature profile and associated seismic
       * velocity that is then subtracted from the seismic velocity at the
       * current pressure temperature conditions
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       *
       * @param values The output vector of depth averaged values. The
       * function takes the pre-existing size of this vector as the number of
       * depth slices.
       */
      void compute_Vp_anomaly(Vector<float> &values) const;

      /**
       * Adjust the pressure variable (which is only determined up to a
       * constant) by adding a constant to it in such a way that the pressure
       * on the surface has a known average value. Whether a face is part of
       * the surface is determined by asking whether its depth of its midpoint
       * (as determined by the geometry model) is less than
       * 1/3*1/sqrt(dim-1)*diameter of the face. For reasonably curved
       * boundaries, this rules out side faces that are perpendicular ot the
       * surface boundary but includes those faces that are along the boundary
       * even if the real boundary is curved.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      void normalize_pressure(LinearAlgebra::BlockVector &vector);

      /**
       * Invert the action of the function above.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      void denormalize_pressure(LinearAlgebra::BlockVector &vector);


      /**
       * Interpolate the given function onto the velocity FE space and write
       * it into the given vector.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      void interpolate_onto_velocity_system(const TensorFunction<1,dim> &func,
                                            LinearAlgebra::Vector &vec);

      /**
       * Set up data structures for null space removal. Called after every
       * mesh refinement.
       *
       * This function is implemented in
       * <code>source/simulator/nullspace.cc</code>.
       */
      void setup_nullspace_removal();

      /**
       * Eliminate the nullspace of the velocity in the given vector. Both
       * vectors are expected to contain the up to date data.
       *
       * @param relevant_dst locally relevant vector for the whole FE, will be
       * filled at the end.
       * @param tmp_distributed_stokes only contains velocity and pressure.
       *
       * This function is implemented in
       * <code>source/simulator/nullspace.cc</code>.
       */
      void remove_nullspace(LinearAlgebra::BlockVector &relevant_dst,
                            LinearAlgebra::BlockVector &tmp_distributed_stokes);

      /**
       * Remove the angular momentum of the given vector
       */
      void remove_net_angular_momentum( LinearAlgebra::BlockVector &relevant_dst, LinearAlgebra::BlockVector &tmp_distributed_stokes);

      /**
       * Remove the linear momentum of the given vector
       */
      void remove_net_linear_momentum( LinearAlgebra::BlockVector &relevant_dst, LinearAlgebra::BlockVector &tmp_distributed_stokes);

      /**
       * Compute the maximal velocity throughout the domain. This is needed to
       * compute the size of the time step.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      double get_maximal_velocity (const LinearAlgebra::BlockVector &solution) const;

      /**
       * Compute the variation (i.e., the difference between maximal and
       * minimal value) of the entropy $(T-\bar T)^2$ where $\bar T$ is the
       * average temperature throughout the domain given as argument to this
       * function.
       *
       * This function is used in computing the artificial diffusion
       * stabilization term.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      double get_entropy_variation (const double average_value,
                                    const AdvectionField &advection_field) const;

      /**
       * Compute the minimal and maximal temperature througout the domain from
       * a solution vector extrapolated from the previous time steps. This is
       * needed to compute the artificial diffusion stabilization terms.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      std::pair<double,double>
      get_extrapolated_advection_field_range (const AdvectionField &advection_field) const;

      /**
       * Compute the size of the next time step from the mesh size and the
       * velocity on each cell. The computed time step has to satisfy the CFL
       * number chosen in the input parameter file on each cell of the mesh.
       * If specified in the parameter file, the time step will be the minimum
       * of the convection *and* conduction time steps. Also returns whether
       * the timestep is dominated by convection (true) or conduction (false).
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      std::pair<double,bool> compute_time_step () const;

      /**
       * Compute the artificial diffusion coefficient value on a cell given
       * the values and gradients of the solution passed as arguments.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      double
      compute_viscosity(internal::Assembly::Scratch::AdvectionSystem<dim> &scratch,
                        const double                        global_u_infty,
                        const double                        global_T_variation,
                        const double                        average_temperature,
                        const double                        global_entropy_variation,
                        const double                        cell_diameter,
                        const AdvectionField     &advection_field) const;

      /**
       * Compute the residual of one advection equation to be used for the
       * artificial diffusion coefficient value on a cell given the values and
       * gradients of the solution passed as arguments.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void
      compute_advection_system_residual(internal::Assembly::Scratch::AdvectionSystem<dim> &scratch,
                                        const double                        average_field,
                                        const AdvectionField     &advection_field,
                                        double                             &max_residual,
                                        double                             &max_velocity,
                                        double                             &max_density,
                                        double                             &max_specific_heat) const;

      /**
       * Extract the values of temperature, pressure, composition and optional
       * strain rate for the current linearization point. These values are
       * stored as input arguments for the material model. The compositional
       * fields are extracted with the individual compositional fields as
       * outer vectors and the values at each quadrature point as inner
       * vectors, but the material model needs it the other way round. Hence,
       * this vector of vectors is transposed.
       *
       * @param[in] input_solution A solution vector (or linear combination of
       * such vectors) with as many entries as there are degrees of freedom in
       * the mesh. It will be evaluated on the cell with which the FEValues
       * object was last re-initialized.
       * @param[in] input_finite_element_values The FEValues object that
       * describes the finite element space in use and that is used to
       * evaluate the solution values at the quadrature points of the current
       * cell.
       * @param[in] compute_strainrate A flag determining whether the strain
       * rate should be computed or not in the output structure.
       * @param[out] material_model_inputs The output structure that contains
       * the solution values evaluated at the quadrature points.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void
      compute_material_model_input_values (const LinearAlgebra::BlockVector                    &input_solution,
                                           const FEValues<dim,dim>                                     &input_finite_element_values,
                                           const bool                                                   compute_strainrate,
                                           typename MaterialModel::Interface<dim>::MaterialModelInputs &material_model_inputs) const;


      /**
       * Return whether the Stokes matrix depends on the values of the
       * solution at the previous time step. This is the case is the
       * coefficients that appear in the matrix (i.e., the viscosity and, in
       * the case of a compressible model, the density) depend on the
       * solution.
       *
       * This function exists to ensure that the Stokes matrix is rebuilt in
       * time steps where it may have changed, while we want to save the
       * effort of rebuilding it whenever we don't need to.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      bool
      stokes_matrix_depends_on_solution () const;

      /**
       * Generate and output some statistics like timing information and
       * memory consumption. Whether this function does anything or not is
       * controlled through the variable aspect::output_parallel_statistics.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      void output_program_stats();

      /**
       * This function is called at the end of each time step and writes the
       * statistics object that contains data like the current time, the
       * number of linear solver iterations, and whatever the postprocessors
       * have generated, to disk.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      void output_statistics();

      /**
       * This routine computes the initial Stokes residual that is
       * needed as a convergence criterion in models with the iterated
       * IMPES solver. We calculate it in the same way as the
       * tolerance for the linear solver, using the norm of the pressure
       * RHS for the pressure part and a residual with zero velocity
       * for the velocity part to get the part of the RHS not balanced
       * by the static pressure.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      double
      compute_initial_stokes_residual();
      /**
       * @}
       */

      /**
       * @name Variables that have to do with input, output, parallel
       * communication and interfacing with other parts of the program.
       * @{
       */
      Parameters                          parameters;
      Introspection<dim>                  introspection;

      MPI_Comm                            mpi_communicator;

      /**
       * This stream will log into the file output/log.txt (used automatically
       * by pcout).
       */
      std::ofstream log_file_stream;

      typedef boost::iostreams::tee_device<std::ostream, std::ofstream> TeeDevice;
      typedef boost::iostreams::stream< TeeDevice > TeeStream;

      TeeDevice iostream_tee_device;
      TeeStream iostream_tee_stream;

      /**
       * Output stream for logging information. Will only output on processor
       * 0.
       */
      ConditionalOStream                  pcout;

      /**
       * An object that stores a bunch of statistics such as the number of
       * linear solver iterations, the time corresponding to each time step,
       * etc, as well as whatever the various postprocessors want to put into
       * it.
       *
       * This variable is written to disk after every time step, by the
       * Simulator::output_statistics() function.
       */
      TableHandler                        statistics;

      Postprocess::Manager<dim>           postprocess_manager;
      TimerOutput                         computing_timer;

      /**
       * In output_statistics(), where we output the statistics object above,
       * we do the actual writing on a separate thread. This variable is the
       * handle we get for this thread so that we can wait for it to finish,
       * either if we want to write the statistics object for the next thread,
       * or if we want to terminate altogether.
       */
      Threads::Thread<>                   output_statistics_thread;
      /**
       * @}
       */

      /**
       * @name Variables that describe the physical setup of the problem
       * @{
       */
      const std::auto_ptr<GeometryModel::Interface<dim> >            geometry_model;
      const IntermediaryConstructorAction                            post_geometry_model_creation_action;
      const std::auto_ptr<MaterialModel::Interface<dim> >            material_model;
      const std::auto_ptr<HeatingModel::Interface<dim> >             heating_model;
      const std::auto_ptr<GravityModel::Interface<dim> >             gravity_model;
      const std::auto_ptr<BoundaryTemperature::Interface<dim> >      boundary_temperature;
      const std::auto_ptr<BoundaryComposition::Interface<dim> >      boundary_composition;
      const std::auto_ptr<InitialConditions::Interface<dim> >        initial_conditions;
      const std::auto_ptr<CompositionalInitialConditions::Interface<dim> > compositional_initial_conditions;
      const std::auto_ptr<AdiabaticConditions::Interface<dim> >      adiabatic_conditions;
      std::map<types::boundary_id,std_cxx1x::shared_ptr<VelocityBoundaryConditions::Interface<dim> > > velocity_boundary_conditions;
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
       * @name Variables related to simulation termination
       * @{
       */
      TerminationCriteria::Manager<dim>                         termination_manager;
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

      MeshRefinement::Manager<dim>                              mesh_refinement_manager;

      const MappingQ<dim>                                       mapping;

      const FESystem<dim>                                       finite_element;

      DoFHandler<dim>                                           dof_handler;

      /**
       * Constraint objects. The first of these describes all constraints that
       * are not time dependent (e.g., hanging nodes, no-normal-flux
       * constraints), whereas the second one is initialized at the top of
       * every time step by copying from the first and then adding to it
       * constraints that are time dependent (e.g., time dependent velocity or
       * temperature boundary conditions).
       *
       * 'constraints' is computed in setup_dofs(), 'current_constraints' is
       * done in compute_current_constraints().
       */
      ConstraintMatrix                                          constraints;
      ConstraintMatrix                                          current_constraints;

      /**
       * The latest correction computed by normalize_pressure(). We store this
       * so we can undo the correction in denormalize_pressure().
       */
      double                                                    pressure_adjustment;

      /**
       * Scaling factor for the pressure as explained in the
       * Kronbichler/Heister/Bangerth paper to ensure that the linear system
       * that results from the Stokes equations is well conditioned.
       */
      double                                                    pressure_scaling;

      /**
       * A variable that determines whether we need to do the correction of
       * the Stokes right hand side vector to ensure that the average
       * divergence is zero. This is necessary for compressible models, but
       * only if there are no in/outflow boundaries.
       */
      bool                           do_pressure_rhs_compatibility_modification;

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

      LinearAlgebra::BlockVector                                current_linearization_point;

      // only used if is_compressible()
      LinearAlgebra::BlockVector                                pressure_shape_function_integrals;



      std_cxx1x::shared_ptr<LinearAlgebra::PreconditionAMG>     Amg_preconditioner;
      std_cxx1x::shared_ptr<LinearAlgebra::PreconditionILU>     Mp_preconditioner;
      std_cxx1x::shared_ptr<LinearAlgebra::PreconditionILU>     T_preconditioner;
//TODO: use n_compositional_field separate preconditioners
      std_cxx1x::shared_ptr<LinearAlgebra::PreconditionILU>     C_preconditioner;

      bool                                                      rebuild_stokes_matrix;
      bool                                                      rebuild_stokes_preconditioner;

      std::vector<LinearAlgebra::Vector> net_rotations_translations;
      /**
       * @}
       */

      /**
       * A member class that isolates the functions and variables that deal
       * with the free surface implementation.  If there are no free surface
       * boundary indicators, then there is no instantiation of this class at
       * all.
       */
      class FreeSurfaceHandler
      {
        public:
          /**
           * Initialize the free surface handler, allowing it to read in
           * relevant parameters as well as giving it a reference to the
           * Simulator that owns it, since it needs to make fairly extensive
           * changes to the internals of the simulator.
           */
          FreeSurfaceHandler(Simulator<dim> &, ParameterHandler &prm);

          /**
           * The main execution step for the free surface implementation. This
           * computes the motion of the free surface, moves the boundary nodes
           * accordingly, redistributes the internal nodes in order to
           * preserve mesh regularity, and calculates the Arbitrary-
           * Lagrangian-Eulerian correction terms for advected quantities.
           */
          void execute();

          /**
           * Allocates and sets up the members of the FreeSurfaceHandler. This
           * is called by Simulator<dim>::setup_dofs()
           */
          void setup_dofs();

          /**
           * Loop over all the mesh vertices and move them so that they are in
           * the positions determined by the free surface implementation.
           * Called in execute(), and also called after redistributing mesh so
           * that the other processes know what has happened to that part of
           * the mesh.
           */
          void displace_mesh();

          /**
           * Apply stabilization to a cell of the system matrix.  The
           * stabilization is only added to cells on a free surface.  The
           * scheme is based on that of Kaus et. al., 2010.  Called during
           * assemly of the system matrix.
           */
          void apply_stabilization (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                    FullMatrix<double> &local_matrix);

          /**
           * Declare parameters for the free surface handling.
           */
          static
          void declare_parameters (ParameterHandler &prm);

          /**
           * Parse parameters for the free surface handling.
           */
          void parse_parameters (ParameterHandler &prm);

        private:
          /**
           * Set the boundary conditions for the solution of the elliptic
           * problem, which computes the displacements of the internal
           * vertices so that the mesh does not become too distored due to
           * motion of the free surface.  Velocities of vertices on the free
           * surface are set to be the normal of the Stokes velocity solution
           * projected onto that surface.  Velocities of vertices on free-slip
           * boundaries are constrained to be tangential to those boundaries.
           * Velocities of vertices on no-slip boundaries are set to be zero.
           */
          void make_constraints ();

          /**
           * Project the normal part of the Stokes velocity solution onto the
           * free surface. Called by make_constraints()
           */
          void project_normal_velocity_onto_boundary (LinearAlgebra::Vector &output);

          /**
           * Actually solve the elliptic problem for the mesh velocitiy.  Just
           * solves a vector Laplacian equation.
           */
          void solve_elliptic_problem ();

          /**
           * From the mesh velocity called in
           * FreeSurfaceHandler::solve_elliptic_problem() we calculate the
           * mesh displacement with mesh_velocity*time_step.  This function
           * also interpolates the mesh velocity onto the finite element space
           * of the Stokes velocity system so that it can be used for ALE
           * corrections.
           */
          void calculate_mesh_displacement ();

          /**
           * Reference to the Simulator object to which a FreeSurfaceHandler
           * instance belongs
           */
          Simulator<dim> &sim;

          /**
           * Finite element for the free surface implementation.  Should be Q1
           */
          const FESystem<dim>                                       free_surface_fe;

          /**
           * DoFHanlder for the free surface implementation
           */
          DoFHandler<dim>                                           free_surface_dof_handler;

          /**
           * Stabilization parameter for the free surface.  Should be between
           * zero and one. A value of zero means no stabilization.  See Kaus
           * et. al. 2010 for more details.
           */
          double free_surface_theta;

          /**
           * BlockVector which stores the mesh velocity interpolated onto the
           * Stokes velocity finite element space.  This is used for ALE
           * corrections.
           */
          LinearAlgebra::BlockVector mesh_velocity;

          /**
           * Vector for storing the positions of the mesh vertices.  This
           * vector is updated by
           * FreeSurfaceHandler::calculate_mesh_displacement(), and is quite
           * important for making sure the mesh stays the same shape upon
           * redistribution of the system.
           */
          LinearAlgebra::Vector mesh_vertices;

          /**
           * The solution of FreeSurfaceHandler::solve_elliptic_problem().
           */
          LinearAlgebra::Vector mesh_vertex_velocity;

          /**
           * The matrix for solving the elliptic problem for moving the
           * internal vertices.
           */
          LinearAlgebra::SparseMatrix mesh_matrix;

          /**
           * IndexSet for the locally owned DoFs for the mesh system
           */
          IndexSet mesh_locally_owned;

          /**
           * IndexSet for the locally relevant DoFs for the mesh system
           */
          IndexSet mesh_locally_relevant;

          /**
           * Storage for the mesh displacement constraints for solving the
           * elliptic problem
           */
          ConstraintMatrix mesh_displacement_constraints;

          /**
           * Storage for the mesh vertex constraints to keep hanging nodes
           * well-behaved
           */
          ConstraintMatrix mesh_vertex_constraints;


          friend class Simulator<dim>;
      };

      /**
       * Shared pointer for an instance of the FreeSurfaceHandler. this way,
       * if we do not need the machinery for doing free surface stuff, we do
       * not even allocate it.
       */
      std_cxx1x::shared_ptr<FreeSurfaceHandler> free_surface;

      friend class boost::serialization::access;
      friend class SimulatorAccess<dim>;
      friend class FreeSurfaceHandler;  //FreeSurfaceHandler needs access to the internals of the Simulator
  };
}


#endif
