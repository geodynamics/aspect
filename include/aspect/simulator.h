//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------
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

#include <aspect/material_model/interface.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/boundary_temperature/interface.h>
#include <aspect/initial_conditions/interface.h>
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
      }

      namespace CopyData
      {
        template <int dim>      struct StokesPreconditioner;
        template <int dim>      struct StokesSystem;
        template <int dim>      struct TemperatureSystem;
      }
    }
  }

  namespace Postprocess
  {
    template <int dim> class SimulatorAccess;
  }


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
        bool                resume_computation;
        double              start_time;
        double              end_time;
        double              CFL_number;
        bool                convert_to_years;
        std::string         output_directory;
        /**
         * @}
         */

        /**
         * @name Parameters that have to do with terms in the model
         * @{
         */
        bool                include_shear_heating;
        double              radiogenic_heating_rate;
        std::vector<int>    fixed_temperature_boundary_indicators;
        std::vector<int>    zero_velocity_boundary_indicators;
        std::vector<int>    tangential_velocity_boundary_indicators;
        std::vector<int>    prescribed_velocity_boundary_indicators;
        /**
         * @}
         */

        /**
         * @name Parameters that have to do with mesh refinement
         * @{
         */
        unsigned int        initial_global_refinement;
        unsigned int        initial_adaptive_refinement;
        double              refinement_fraction;
        double              coarsening_fraction;
        std::vector<double> additional_refinement_times;
        unsigned int        adaptive_refinement_interval;
        /**
         * @}
         */

        /**
         * @name Parameters that have to do with the stabilization of transport equations
         * @{
         */
        double              stabilization_alpha;
        double              stabilization_c_R;
        double              stabilization_beta;
        /**
         * @}
         */

        /**
         * @name Parameters that have to do with spatial discretizations
         * @{
         */
        unsigned int        stokes_velocity_degree;
        bool                use_locally_conservative_discretization;
        unsigned int        temperature_degree;
        /**
         * @}
         */
      };

      /**
       * Constructor
       *
       * @param prm The run-time parameter object from which this class
       * obtains its settings.
       **/
      Simulator (ParameterHandler &prm);


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
       * The function that runs the overall algorithm. It
       * contains the loop over all time steps as well as the logic
       * of what to do when before the loop starts and within the time
       * loop.
       **/
      void run ();

    private:
      /**
       * @name Top-level functions in the overall flow of the numerical algorithm
       * @{
       */

      /**
       * The function that sets up the DoFHandler objects,
       **/
      void setup_dofs ();
      void set_initial_temperature_field ();
      void compute_initial_pressure_field ();
      void start_timestep ();
      void assemble_stokes_preconditioner ();
      void build_stokes_preconditioner ();
      void assemble_stokes_system ();
      void assemble_temperature_system ();
      void solve ();
      void postprocess ();
      void refine_mesh (const unsigned int max_grid_level);
      /**
       * @}
       */

      /**
       * @name Functions used in saving the state of the program and restarting from a saved state
       * @{
       */
      void create_snapshot();
      void resume_from_snapshot();

      template <class Archive>
      void serialize (Archive &ar, const unsigned int version);
      /**
       * @}
       */

      /**
       * @name Functions used in setting up linear systems
       * @{
       */
      void setup_stokes_matrix (const std::vector<IndexSet> &stokes_partitioning);
      void setup_stokes_preconditioner (const std::vector<IndexSet> &stokes_partitioning);
      void setup_temperature_matrix (const IndexSet &temperature_partitioning);
      /**
       * @}
       */

      /**
       * @name Functions used in the assembly of linear systems
       * @{
       */
      void
      local_assemble_stokes_preconditioner (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                            internal::Assembly::Scratch::StokesPreconditioner<dim> &scratch,
                                            internal::Assembly::CopyData::StokesPreconditioner<dim> &data);

      void
      copy_local_to_global_stokes_preconditioner (const internal::Assembly::CopyData::StokesPreconditioner<dim> &data);


      void
      local_assemble_stokes_system (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                    internal::Assembly::Scratch::StokesSystem<dim>  &scratch,
                                    internal::Assembly::CopyData::StokesSystem<dim> &data);

      void
      copy_local_to_global_stokes_system (const internal::Assembly::CopyData::StokesSystem<dim> &data);


      void
      local_assemble_temperature_system (const std::pair<double,double> global_T_range,
                                         const double                   global_max_velocity,
                                         const double                   global_entropy_variation,
                                         const typename DoFHandler<dim>::active_cell_iterator &cell,
                                         internal::Assembly::Scratch::TemperatureSystem<dim>  &scratch,
                                         internal::Assembly::CopyData::TemperatureSystem<dim> &data);

      void
      copy_local_to_global_temperature_system (const internal::Assembly::CopyData::TemperatureSystem<dim> &data);
      /**
       * @}
       */

      /**
       * @name Helper functions
       * @{
       */
      void make_pressure_rhs_compatible(TrilinosWrappers::MPI::BlockVector &vector);

      void normalize_pressure(TrilinosWrappers::MPI::BlockVector &vector);
      double get_maximal_velocity () const;
      double get_entropy_variation (const double average_temperature) const;
      std::pair<double,double> get_extrapolated_temperature_range () const;
      double compute_time_step () const;

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
                        const double                        global_u_infty,
                        const double                        global_T_variation,
                        const double                        average_temperature,
                        const double                        global_entropy_variation,
                        const std::vector<Point<dim> >     &evaluation_points,
                        const double                        cell_diameter) const;
      /**
       * @}
       */

      /**
       * @name Variables that have to do with input and output
       * @{
       */
      Parameters                          parameters;
      ConditionalOStream                  pcout;
      TableHandler                        statistics;
      Postprocess::Manager<dim>           postprocess_manager;
      TimerOutput computing_timer;
      /**
       * @}
       */

      /**
       * @name Variables that describe the physical setup of the problem
       * @{
       */
      const std::auto_ptr<const GeometryModel::Interface<dim> > geometry_model;
      const std::auto_ptr<MaterialModel::Interface<dim> >       material_model;
      const std::auto_ptr<GravityModel::Interface<dim> >        gravity_model;
      const std::auto_ptr<BoundaryTemperature::Interface<dim> > boundary_temperature;
      std::auto_ptr<const InitialConditions::Interface<dim> >   initial_conditions;
      std::auto_ptr<const AdiabaticConditions<dim> >            adiabatic_conditions;
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
       * @name Variables that describe the spatial discretization
       * @{
       */
      parallel::distributed::Triangulation<dim>                 triangulation;
      double                                                    global_Omega_diameter;
      double                                                    global_volume;

      const MappingQ<dim>                                       mapping;

      const FESystem<dim>                                       stokes_fe;

      DoFHandler<dim>                                           stokes_dof_handler;
      ConstraintMatrix                                          stokes_constraints;
      ConstraintMatrix                                          current_stokes_constraints;

      double                                                    pressure_scaling;

      FE_Q<dim>                                                 temperature_fe;
      DoFHandler<dim>                                           temperature_dof_handler;
      ConstraintMatrix                                          temperature_constraints;
//TODO: Have a current_temperature_constraints object, just like we do for the Stokes system
      /**
       * @}
       */


      /**
       * @name Variables that describe the linear systems and solution vectors
       * @{
       */
      TrilinosWrappers::BlockSparseMatrix                       stokes_matrix;
      TrilinosWrappers::BlockSparseMatrix                       stokes_preconditioner_matrix;

      TrilinosWrappers::MPI::BlockVector                        stokes_solution;
      TrilinosWrappers::MPI::BlockVector                        old_stokes_solution;
      TrilinosWrappers::MPI::BlockVector                        stokes_rhs;

      // only used if is_compressible()
      TrilinosWrappers::MPI::BlockVector                        pressure_shape_function_integrals;

      TrilinosWrappers::SparseMatrix                            temperature_matrix;

      TrilinosWrappers::MPI::Vector                             temperature_solution;
      TrilinosWrappers::MPI::Vector                             old_temperature_solution;
      TrilinosWrappers::MPI::Vector                             old_old_temperature_solution;
      TrilinosWrappers::MPI::Vector                             temperature_rhs;



      std_cxx1x::shared_ptr<TrilinosWrappers::PreconditionAMG>  Amg_preconditioner;
      std_cxx1x::shared_ptr<TrilinosWrappers::PreconditionILU>  Mp_preconditioner;
      std_cxx1x::shared_ptr<TrilinosWrappers::PreconditionILU>  T_preconditioner;

      bool                                                      rebuild_stokes_matrix;
      bool                                                      rebuild_stokes_preconditioner;
      /**
       * @}
       */

      friend class boost::serialization::access;
      friend class Postprocess::SimulatorAccess<dim>;
  };
}


#endif
