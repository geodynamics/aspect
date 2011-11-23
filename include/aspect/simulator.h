//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------
#ifndef __aspect__simulator_h
#define __aspect__simulator_h

#include <deal.II/base/timer.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/conditional_ostream.h>

#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>

#include <aspect/material_model/interface.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/gravity_model/interface.h>
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

  template <int dim>
  class Simulator : public Subscriptor
  {
    public:
      struct Parameters;

      Simulator (ParameterHandler &prm);

      static
      void declare_parameters (ParameterHandler &prm);

      void run ();

    private:
      void setup_dofs ();
      void assemble_stokes_preconditioner ();
      void build_stokes_preconditioner ();
      void assemble_stokes_system ();
      void assemble_temperature_system ();
      void set_initial_temperature_field ();
      void compute_initial_pressure_field ();
      double get_maximal_velocity () const;
      double get_cfl_number () const;
      double get_entropy_variation (const double average_temperature) const;
      std::pair<double,double> get_extrapolated_temperature_range () const;
      void solve ();
      void postprocess ();
      void refine_mesh (const unsigned int max_grid_level);

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

    public:
      struct Parameters
      {
        Parameters (ParameterHandler &prm);

        static void declare_parameters (ParameterHandler &prm);
        void parse_parameters (ParameterHandler &prm);

        double end_time;
        double time_step_scaling;

        bool resume_computation;

        unsigned int initial_global_refinement;
        unsigned int initial_adaptive_refinement;
        double       refinement_fraction;
        double       coarsening_fraction;

        std::vector<double> additional_refinement_times;

        unsigned int adaptive_refinement_interval;

        double       stabilization_alpha;
        double       stabilization_c_R;
        double       stabilization_beta;

        unsigned int stokes_velocity_degree;
        bool         use_locally_conservative_discretization;

        unsigned int temperature_degree;

        bool         include_shear_heating;
        double       radiogenic_heating_rate;
      };

    private:
      Parameters                          parameters;
      ConditionalOStream                  pcout;

      const std::auto_ptr<const GeometryModel::Interface<dim> > geometry_model;
      const std::auto_ptr<MaterialModel::Interface<dim> > material_model;
      const std::auto_ptr<GravityModel::Interface<dim> > gravity_model;
      std::auto_ptr<const InitialConditions::Interface<dim> > initial_conditions;


      std::auto_ptr<const AdiabaticConditions<dim> >      adiabatic_conditions;
      double                                              global_Omega_diameter;
      double                                              global_volume;

      Postprocess::Manager<dim>           postprocess_manager;
      TableHandler                        statistics;

      parallel::distributed::Triangulation<dim> triangulation;

      double                              pressure_scaling;

      const MappingQ<dim>                 mapping;

      const FESystem<dim>                 stokes_fe;

      DoFHandler<dim>                     stokes_dof_handler;
      ConstraintMatrix                    stokes_constraints;

      TrilinosWrappers::BlockSparseMatrix stokes_matrix;
      TrilinosWrappers::BlockSparseMatrix stokes_preconditioner_matrix;

      TrilinosWrappers::MPI::BlockVector  stokes_solution;
      TrilinosWrappers::MPI::BlockVector  old_stokes_solution;
      TrilinosWrappers::MPI::BlockVector  stokes_rhs;
      TrilinosWrappers::MPI::BlockVector  stokes_rhs_helper;


      FE_Q<dim>                           temperature_fe;
      DoFHandler<dim>                     temperature_dof_handler;
      ConstraintMatrix                    temperature_constraints;

      TrilinosWrappers::SparseMatrix      temperature_matrix;

      TrilinosWrappers::MPI::Vector       temperature_solution;
      TrilinosWrappers::MPI::Vector       old_temperature_solution;
      TrilinosWrappers::MPI::Vector       old_old_temperature_solution;
      TrilinosWrappers::MPI::Vector       temperature_rhs;


      double time;
      double time_step;
      double old_time_step;
      unsigned int timestep_number;

      std_cxx1x::shared_ptr<TrilinosWrappers::PreconditionAMG> Amg_preconditioner;
      std_cxx1x::shared_ptr<TrilinosWrappers::PreconditionILU> Mp_preconditioner;
      std_cxx1x::shared_ptr<TrilinosWrappers::PreconditionILU> T_preconditioner;

      bool rebuild_stokes_matrix;
      bool rebuild_stokes_preconditioner;

      TimerOutput computing_timer;

      void setup_stokes_matrix (const std::vector<IndexSet> &stokes_partitioning);
      void setup_stokes_preconditioner (const std::vector<IndexSet> &stokes_partitioning);
      void setup_temperature_matrix (const IndexSet &temperature_partitioning);

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

      void normalize_pressure(TrilinosWrappers::MPI::BlockVector &vector);


      void create_snapshot();
      void resume_from_snapshot();

      template <class Archive>
      void serialize (Archive &ar, const unsigned int version);

      void make_pressure_rhs_compatible(TrilinosWrappers::MPI::BlockVector &vector,
                                        const TrilinosWrappers::MPI::BlockVector &helper);

      friend class boost::serialization::access;
      friend class Postprocess::SimulatorAccess<dim>;
  };
}


#endif
