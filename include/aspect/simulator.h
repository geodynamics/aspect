//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011, 2012 by the authors of the ASPECT code
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
        double              surface_pressure;
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
        std::string         refinement_strategy;
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
       * A function that is responsible for initializing the temperature field
       * before the first time step. This temperature field then serves as the
       * temperature from which the velocity is computed during the first time
       * step, and is subsequently overwritten by the temperature field one gets
       * by advancing by one time step.
       *
       * This function is implemented in
       * <code>source/simulator/initial_conditions.cc</code>.
       */
      void set_initial_temperature_field ();

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
       * Initiate the assembly of the Stokes preconditioner matrix via
       * assemble_stokes_preconditoner(), then set up the data structures
       * to actually build a preconditioner from this matrix.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void build_stokes_preconditioner ();

      /**
       * Initiate the assembly of the Stokes matrix and right hand side.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void assemble_stokes_system ();

      /**
       * Initiate the assembly of the temperature matrix and right hand side
       * and build a preconditioner for the matrix.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void assemble_temperature_system ();

      /**
       * Solve the temperature linear system.
       *
       * This function is implemented in
       * <code>source/simulator/solver.cc</code>.
       */
      void solve_temperature ();

      /**
       * Solve the Stokes linear system.
       *
       * This function is implemented in
       * <code>source/simulator/solver.cc</code>.
       */
      void solve_stokes ();

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
      void make_pressure_rhs_compatible(TrilinosWrappers::MPI::BlockVector &vector);

      /**
       * Computes a running average of a vector.
       * http://en.wikipedia.org/wiki/Moving_average
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      void compute_running_average(std::vector<double> &values, const int npoints) const;

      /**
       * Compute a depth average of the current temperature. The function
       * fills a vector that contains average temperatures over slices of the
       * domain of same depth. The function resizes the output vector to match
       * the number of depth slices.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      void compute_depth_average_temperature(std::vector<double> &values) const;

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
      void normalize_pressure(TrilinosWrappers::MPI::BlockVector &vector);

      /**
       * Compute the maximal velocity throughout the domain. This is needed
       * to compute the size of the time step.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      double get_maximal_velocity (const TrilinosWrappers::MPI::BlockVector &solution) const;

      /**
       * Compute the variation (i.e., the difference between maximal and
       * minimal value) of the entropy $(T-\bar T)^2$ where $\bar T$ is the
       * average temperature throughout the domain given as argument to
       * this function.
       *
       * This function is used in computing the artificial diffusion
       * stabilization term.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      double get_entropy_variation (const double average_temperature) const;

      /**
       * Compute the minimal and maximal temperature througout the domain from a
       * solution vector extrapolated from the previous time steps. This is needed
       * to compute the artificial diffusion stabilization terms.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      std::pair<double,double> get_extrapolated_temperature_range () const;

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
                        const double                        global_u_infty,
                        const double                        global_T_variation,
                        const double                        average_temperature,
                        const double                        global_entropy_variation,
                        const std::vector<Point<dim> >     &evaluation_points,
                        const double                        cell_diameter) const;

      /**
        * generate and output some statistics like timing information and memory consumption.
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
      TableHandler                        statistics;
      Postprocess::Manager<dim>           postprocess_manager;
      TimerOutput                         computing_timer;
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
      TrilinosWrappers::BlockSparseMatrix                       system_matrix;
      TrilinosWrappers::BlockSparseMatrix                       system_preconditioner_matrix;

      TrilinosWrappers::MPI::BlockVector                        solution;
      TrilinosWrappers::MPI::BlockVector                        old_solution;
      TrilinosWrappers::MPI::BlockVector                        old_old_solution;
      TrilinosWrappers::MPI::BlockVector                        system_rhs;

      // only used if is_compressible()
      TrilinosWrappers::MPI::BlockVector                        pressure_shape_function_integrals;



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
