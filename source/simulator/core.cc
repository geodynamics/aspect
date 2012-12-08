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


#include <aspect/simulator.h>
#include <aspect/global.h>

#include <deal.II/base/index_set.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/derivative_approximation.h>
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/distributed/grid_refinement.h>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <locale>
#include <string>


using namespace dealii;


namespace aspect
{
  namespace
  {
    /**
     * Return whether t is an element of the given container object.
     */
    template <typename Container>
    bool is_element (const typename Container::value_type &t,
                     const Container                      &container)
    {
      for (typename Container::const_iterator p = container.begin();
           p != container.end();
           ++p)
        if (*p == t)
          return true;

      return false;
    }
  }

  /**
   * Constructor. Initialize all member variables.
   **/
  template <int dim>
  Simulator<dim>::Simulator (const MPI_Comm mpi_communicator_,
                             ParameterHandler &prm)
    :
    parameters (prm),
    mpi_communicator (Utilities::MPI::duplicate_communicator (mpi_communicator_)),
    pcout (std::cout,
           (Utilities::MPI::
            this_mpi_process(mpi_communicator)
            == 0)),

    computing_timer (pcout, TimerOutput::summary,
                     TimerOutput::wall_times),

    geometry_model (GeometryModel::create_geometry_model<dim>(prm)),
    material_model (MaterialModel::create_material_model<dim>(prm)),
    gravity_model (GravityModel::create_gravity_model<dim>(prm)),
    boundary_temperature (BoundaryTemperature::create_boundary_temperature<dim>(prm)),
    initial_conditions (InitialConditions::create_initial_conditions (prm,
                                                                      *geometry_model,
                                                                      *boundary_temperature,
                                                                      *adiabatic_conditions)),
    compositional_initial_conditions (CompositionalInitialConditions::create_initial_conditions (prm,
                                      *geometry_model)),

    time (std::numeric_limits<double>::quiet_NaN()),
    time_step (0),
    old_time_step (0),
    timestep_number (0),

    triangulation (mpi_communicator,
                   typename Triangulation<dim>::MeshSmoothing
                   (Triangulation<dim>::smoothing_on_refinement |
                    Triangulation<dim>::smoothing_on_coarsening),
                   parallel::distributed::Triangulation<dim>::mesh_reconstruction_after_repartitioning),

    mapping (4),

    finite_element(FE_Q<dim>(parameters.stokes_velocity_degree),
                   dim,
                   (parameters.use_locally_conservative_discretization
                    ?
                    static_cast<const FiniteElement<dim> &>
                    (FE_DGP<dim>(parameters.stokes_velocity_degree-1))
                    :
                    static_cast<const FiniteElement<dim> &>
                    (FE_Q<dim>(parameters.stokes_velocity_degree-1))),
                   1,
                   FE_Q<dim>(parameters.temperature_degree),
                   1,
                   FE_Q<dim>(parameters.composition_degree),
                   parameters.n_compositional_fields),

    dof_handler (triangulation),

    rebuild_stokes_matrix (true),
    rebuild_stokes_preconditioner (true)
  {
    computing_timer.enter_section("Initialization");
    // first do some error checking for the parameters we got
    {
      // make sure velocity boundary indicators don't appear in multiple lists
      std::set<types::boundary_id_t> boundary_indicator_lists[3]
        = { parameters.zero_velocity_boundary_indicators,
            parameters.tangential_velocity_boundary_indicators,
            std::set<types::boundary_id_t>()
          };
      for (std::map<types::boundary_id_t,std::string>::const_iterator
           p = parameters.prescribed_velocity_boundary_indicators.begin();
           p != parameters.prescribed_velocity_boundary_indicators.end();
           ++p)
        boundary_indicator_lists[2].insert (p->first);

      // for each combination of boundary indicator lists, make sure that the
      // intersection is empty
      for (unsigned int i=0; i<sizeof(boundary_indicator_lists)/sizeof(boundary_indicator_lists[0]); ++i)
        for (unsigned int j=i+1; j<sizeof(boundary_indicator_lists)/sizeof(boundary_indicator_lists[0]); ++j)
          {
            std::set<types::boundary_id_t> intersection;
            std::set_intersection (boundary_indicator_lists[i].begin(),
                                   boundary_indicator_lists[i].end(),
                                   boundary_indicator_lists[j].begin(),
                                   boundary_indicator_lists[j].end(),
                                   std::inserter(intersection, intersection.end()));
            AssertThrow (intersection.empty(),
                         ExcMessage ("Some velocity boundary indicators are listed as having more "
                                     "than one type in the input file."));
          }

      // next make sure that all listed indicators are actually used by
      // this geometry
      const std::set<types::boundary_id_t> all_boundary_indicators
        = geometry_model->get_used_boundary_indicators();
      for (unsigned int i=0; i<sizeof(boundary_indicator_lists)/sizeof(boundary_indicator_lists[0]); ++i)
        for (typename std::set<types::boundary_id_t>::const_iterator
             p = boundary_indicator_lists[i].begin();
             p != boundary_indicator_lists[i].end(); ++p)
          AssertThrow (all_boundary_indicators.find (*p)
                       != all_boundary_indicators.end(),
                       ExcMessage ("One of the boundary indicators listed in the input file "
                                   "is not used by the geometry model."));

      // now do the same for the fixed temperature indicators
      for (typename std::set<types::boundary_id_t>::const_iterator
           p = parameters.fixed_temperature_boundary_indicators.begin();
           p != parameters.fixed_temperature_boundary_indicators.end(); ++p)
        AssertThrow (all_boundary_indicators.find (*p)
                     != all_boundary_indicators.end(),
                     ExcMessage ("One of the boundary indicators listed in the input file "
                                 "is not used by the geometry model."));
    }

    // continue with initializing members that can't be initialized for one reason
    // or another in the member initializer list above

    // if any plugin wants access to the Simulator by deriving from SimulatorAccess, initialize it:
    if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(&*geometry_model))
      sim->initialize (*this);
    if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(&*material_model))
      sim->initialize (*this);
    if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(&*gravity_model))
      sim->initialize (*this);
    if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(&*boundary_temperature))
      sim->initialize (*this);
    if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(&*initial_conditions))
      sim->initialize (*this);

    postprocess_manager.parse_parameters (prm);
    postprocess_manager.initialize (*this);

    mesh_refinement_manager.parse_parameters (prm);
    mesh_refinement_manager.initialize (*this);

    geometry_model->create_coarse_mesh (triangulation);
    global_Omega_diameter = GridTools::diameter (triangulation);

    adiabatic_conditions.reset (new AdiabaticConditions<dim>(*geometry_model,
                                                             *gravity_model,
                                                             *material_model,
                                                             *compositional_initial_conditions,
                                                             parameters.surface_pressure,
                                                             parameters.adiabatic_surface_temperature,
                                                             parameters.n_compositional_fields));
    for (std::map<types::boundary_id_t,std::string>::const_iterator
         p = parameters.prescribed_velocity_boundary_indicators.begin();
         p != parameters.prescribed_velocity_boundary_indicators.end();
         ++p)
      {
        VelocityBoundaryConditions::Interface<dim> *bv
          = VelocityBoundaryConditions::create_velocity_boundary_conditions
            (p->second,
             prm,
             *geometry_model);
        velocity_boundary_conditions[p->first].reset (bv);
      }

    pressure_scaling = material_model->reference_viscosity() / geometry_model->length_scale();

    // make sure that we don't have to fill every column of the statistics
    // object in each time step.
    statistics.set_auto_fill_mode(true);

    // finally produce a record of the run-time parameters by writing
    // the currently used values into a file
    // Only write the parameter files on the root node to avoid file system conflicts
    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
      {
        std::ofstream prm_out ((parameters.output_directory + "parameters.prm").c_str());
        AssertThrow (prm_out,
                     ExcMessage (std::string("Couldn't open file <") +
                                 parameters.output_directory + "parameters.prm>."));
        prm.print_parameters(prm_out, ParameterHandler::Text);
      }
    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
      {
        std::ofstream prm_out ((parameters.output_directory + "parameters.tex").c_str());
        AssertThrow (prm_out,
                     ExcMessage (std::string("Couldn't open file <") +
                                 parameters.output_directory + "parameters.tex>."));
        prm.print_parameters(prm_out, ParameterHandler::LaTeX);
      }
    computing_timer.exit_section();
  }


  /**
   * Destructor.
   **/
  template <int dim>
  Simulator<dim>::~Simulator ()
  {
    // wait if there is a thread that's still writing the statistics
    // object (set from the output_statistics() function)
    output_statistics_thread.join();
  }


  namespace
  {
    /**
     * Conversion object where one can provide a function that returns
     * a tensor for the velocity at a given point and it returns something
     * that matches the dealii::Function interface with a number of output
     * components equal to the number of components of the finite element
     * in use.
     */
    template <int dim>
    class VectorFunctionFromVelocityFunctionObject : public Function<dim>
    {
      public:
        /**
         * Given a function object that takes a Point and returns a Tensor<1,dim>,
         * convert this into an object that matches the Function@<dim@>
         * interface.
         *
        * @param n_comp_fields The number of compositional fields used in
        *        this vector-valued problem. Needed to determine
        *        the total size of the output vectors.
         * @param function_object The scalar function that will form one component
         *     of the resulting Function object.
         */
        VectorFunctionFromVelocityFunctionObject (const unsigned int n_comp_fields,
                                                  const std_cxx1x::function<Tensor<1,dim> (const Point<dim> &)> &function_object);

        /**
         * Return the value of the
         * function at the given
         * point. Returns the value the
         * function given to the constructor
         * produces for this point.
         */
        virtual double value (const Point<dim>   &p,
                              const unsigned int  component = 0) const;

        /**
         * Return all components of a
         * vector-valued function at a
         * given point.
         *
         * <tt>values</tt> shall have the right
         * size beforehand,
         * i.e. #n_components.
         */
        virtual void vector_value (const Point<dim>   &p,
                                   Vector<double>     &values) const;

      private:
        /**
         * The function object which we call when this class's value() or
         * value_list() functions are called.
         **/
        const std_cxx1x::function<Tensor<1,dim> (const Point<dim> &)> function_object;
    };


    template <int dim>
    VectorFunctionFromVelocityFunctionObject<dim>::
    VectorFunctionFromVelocityFunctionObject
    (const unsigned int n_comp_fields,
     const std_cxx1x::function<Tensor<1,dim> (const Point<dim> &)> &function_object)
      :
      Function<dim>(dim+2+n_comp_fields),
      function_object (function_object)
    {
    }



    template <int dim>
    double
    VectorFunctionFromVelocityFunctionObject<dim>::value (const Point<dim> &p,
                                                          const unsigned int component) const
    {
      Assert (component < this->n_components,
              ExcIndexRange (component, 0, this->n_components));

      if (component < dim)
        {
          const Tensor<1,dim> v = function_object(p);
          return v[component];
        }
      else
        return 0;
    }



    template <int dim>
    void
    VectorFunctionFromVelocityFunctionObject<dim>::
    vector_value (const Point<dim>   &p,
                  Vector<double>     &values) const
    {
      AssertDimension(values.size(), this->n_components);

      // set everything to zero, and then the right components to their correct values
      values = 0;

      const Tensor<1,dim> v = function_object(p);
      for (unsigned int d=0; d<dim; ++d)
        values(d) = v[d];
    }
  }


  template <int dim>
  void
  Simulator<dim>::
  start_timestep ()
  {
    // first produce some output for the screen to show where we are
    if (parameters.convert_to_years == true)
      pcout << "*** Timestep " << timestep_number
            << ":  t=" << time/year_in_seconds
            << " years"
            << std::endl;
    else
      pcout << "*** Timestep " << timestep_number
            << ":  t=" << time
            << " seconds"
            << std::endl;

    // set global statistics about this time step
    statistics.add_value("Time step number", timestep_number);
    if (parameters.convert_to_years == true)
      statistics.add_value("Time (years)", time / year_in_seconds);
    else
      statistics.add_value("Time (seconds)", time);
    statistics.add_value("Number of mesh cells",
                         triangulation.n_global_active_cells());

    std::vector<unsigned int> system_sub_blocks (dim+2+parameters.n_compositional_fields,0);
    system_sub_blocks[dim] = 1;
    system_sub_blocks[dim+1] = 2;
    for (unsigned int i=dim+2; i<dim+2+parameters.n_compositional_fields; ++i)
      system_sub_blocks[i] = 3;
    std::vector<unsigned int> system_dofs_per_block (3+((parameters.n_compositional_fields>0)?1:0));
    DoFTools::count_dofs_per_block (dof_handler, system_dofs_per_block,
                                    system_sub_blocks);

    statistics.add_value("Number of Stokes degrees of freedom",
                         system_dofs_per_block[0]+system_dofs_per_block[1]);
    statistics.add_value("Number of temperature degrees of freedom",
                         system_dofs_per_block[2]);
    if (parameters.n_compositional_fields > 0)
      statistics.add_value("Number of composition degrees of freedom",
                           system_dofs_per_block[3]);


    // then interpolate the current boundary velocities. this adds to
    // the current_constraints object we already have
    {
      IndexSet system_relevant_set;
      DoFTools::extract_locally_relevant_dofs (dof_handler,
                                               system_relevant_set);

      current_constraints.clear ();
      current_constraints.reinit (system_relevant_set);
      current_constraints.merge (constraints);

      // set the current time and do the interpolation
      // for the prescribed velocity fields
      std::vector<bool> velocity_mask (dim+2+parameters.n_compositional_fields, true);
      for (unsigned int i=dim; i<dim+2+parameters.n_compositional_fields; ++i)
        velocity_mask[i] = false;
      for (typename std::map<types::boundary_id_t,std_cxx1x::shared_ptr<VelocityBoundaryConditions::Interface<dim> > >::iterator
           p = velocity_boundary_conditions.begin();
           p != velocity_boundary_conditions.end(); ++p)
        {
          p->second->set_current_time (time);
          VectorFunctionFromVelocityFunctionObject<dim> vel
          (parameters.n_compositional_fields,
           std_cxx1x::bind (&VelocityBoundaryConditions::Interface<dim>::boundary_velocity,
                            p->second,
                            std_cxx1x::_1));
          VectorTools::interpolate_boundary_values (dof_handler,
                                                    p->first,
                                                    vel,
                                                    current_constraints,
                                                    velocity_mask);
        }
      current_constraints.close();
    }

    //TODO: do this in a more efficient way (TH)? we really only need
    // to make sure that the time dependent velocity boundary conditions
    // end up in the right hand side in the righth way; we currently do
    // that by re-assembling the entire system
    if (!velocity_boundary_conditions.empty())
      rebuild_stokes_matrix = rebuild_stokes_preconditioner = true;

    // notify different system components that we started the next time step
    material_model->update();
    gravity_model->update();
  }



  template <int dim>
  void
  Simulator<dim>::
  setup_system_matrix (const std::vector<IndexSet> &system_partitioning)
  {
    system_matrix.clear ();

    TrilinosWrappers::BlockSparsityPattern sp (system_partitioning,
                                               mpi_communicator);

    Table<2,DoFTools::Coupling> coupling (dim+2+parameters.n_compositional_fields, dim+2+parameters.n_compositional_fields);

    // determine which blocks should be fillable in the matrix.
    // note:
    // - all velocities couple with all velocities
    // - pressure couples with all velocities and the other way
    //   around
    // - temperature only couples with itself when using the impes
    //   scheme
    for (unsigned int c=0; c<dim; ++c)
      for (unsigned int d=0; d<dim; ++d)
        coupling[c][d] = DoFTools::always;
    for (unsigned int c=0; c<dim; ++c)
      coupling[c][dim] = coupling[dim][c] = DoFTools::always;
    coupling[dim+1][dim+1] = DoFTools::always;
    for (unsigned int c=dim+2; c<dim+2+parameters.n_compositional_fields; ++c)
      coupling[c][c] = DoFTools::always;

    DoFTools::make_sparsity_pattern (dof_handler,
                                     coupling, sp,
                                     constraints, false,
                                     Utilities::MPI::
                                     this_mpi_process(mpi_communicator));
    sp.compress();

    system_matrix.reinit (sp);
  }



  template <int dim>
  void Simulator<dim>::
  setup_system_preconditioner (const std::vector<IndexSet> &system_partitioning)
  {
    Amg_preconditioner.reset ();
    Mp_preconditioner.reset ();
    T_preconditioner.reset ();
    C_preconditioner.reset ();

    system_preconditioner_matrix.clear ();

    TrilinosWrappers::BlockSparsityPattern sp (system_partitioning,
                                               mpi_communicator);

    Table<2,DoFTools::Coupling> coupling (dim+2+parameters.n_compositional_fields, dim+2+parameters.n_compositional_fields);
    for (unsigned int c=0; c<dim+2+parameters.n_compositional_fields; ++c)
      for (unsigned int d=0; d<dim+2+parameters.n_compositional_fields; ++d)
        if (c == d)
          coupling[c][d] = DoFTools::always;
        else
          coupling[c][d] = DoFTools::none;

    DoFTools::make_sparsity_pattern (dof_handler,
                                     coupling, sp,
                                     constraints, false,
                                     Utilities::MPI::
                                     this_mpi_process(mpi_communicator));
    sp.compress();

    system_preconditioner_matrix.reinit (sp);
  }



  template <int dim>
  void Simulator<dim>::setup_dofs ()
  {
    computing_timer.enter_section("Setup dof systems");

    dof_handler.distribute_dofs(finite_element);

    // Renumber the DoFs hierarchical so that we get the
    // same numbering if we resume the computation. This
    // is because the numbering depends on the order the
    // cells are created.
    DoFRenumbering::hierarchical (dof_handler);
    std::vector<unsigned int> system_sub_blocks (dim+2+parameters.n_compositional_fields,0);
    system_sub_blocks[dim] = 1;
    system_sub_blocks[dim+1] = 2;
    for (unsigned int i=dim+2; i<dim+2+parameters.n_compositional_fields; ++i)
      system_sub_blocks[i] = i-dim+1;

    DoFRenumbering::component_wise (dof_handler, system_sub_blocks);

    // set up the introspection object that stores all sorts of
    // information about components of the finite element, component
    // masks, etc
    setup_introspection();

    std::vector<unsigned int> system_dofs_per_block (3+parameters.n_compositional_fields);
    DoFTools::count_dofs_per_block (dof_handler, system_dofs_per_block,
                                    system_sub_blocks);

    const unsigned int n_u = system_dofs_per_block[0],
                       n_p = system_dofs_per_block[1],
                       n_T = system_dofs_per_block[2];
    unsigned int       n_C_sum = 0;
    std::vector<unsigned int> n_C (parameters.n_compositional_fields+1);
    for (unsigned int c=0; c<parameters.n_compositional_fields; ++c)
      {
        n_C[c] = system_dofs_per_block[c+3];
        n_C_sum += n_C[c];
      }


    // print dof numbers with 1000s
    // separator since they are frequently
    // large
    std::locale s = pcout.get_stream().getloc();
    // Creating std::locale with an empty string causes problems
    // on some platforms, so catch the exception and ignore
    try
      {
        pcout.get_stream().imbue(std::locale(""));
      }
    catch (std::runtime_error e)
      {
        // If the locale doesn't work, just give up
      }
    pcout << "Number of active cells: "
          << triangulation.n_global_active_cells()
          << " (on "
          << triangulation.n_levels()
          << " levels)"
          << std::endl
          << "Number of degrees of freedom: "
          << n_u + n_p + n_T + n_C_sum
          << " (" << n_u << '+' << n_p << '+'<< n_T;

    for (unsigned int c=0; c<parameters.n_compositional_fields; ++c)
      pcout << '+' << n_C[c];

    pcout <<')'
          << std::endl
          << std::endl;
    pcout.get_stream().imbue(s);


    // now also compute the various partitionings between processors and blocks
    // of vectors and matrices

    std::vector<IndexSet> system_partitioning, system_relevant_partitioning;
    IndexSet system_relevant_set;
    {
      IndexSet system_index_set = dof_handler.locally_owned_dofs();
      system_partitioning.push_back(system_index_set.get_view(0,n_u));
      system_partitioning.push_back(system_index_set.get_view(n_u,n_u+n_p));
      system_partitioning.push_back(system_index_set.get_view(n_u+n_p,n_u+n_p+n_T));

      DoFTools::extract_locally_relevant_dofs (dof_handler,
                                               system_relevant_set);
      system_relevant_partitioning.push_back(system_relevant_set.get_view(0,n_u));
      system_relevant_partitioning.push_back(system_relevant_set.get_view(n_u,n_u+n_p));
      system_relevant_partitioning.push_back(system_relevant_set.get_view(n_u+n_p, n_u+n_p+n_T));

      {
        unsigned int n_C_so_far = 0;

        for (unsigned int c=0; c<parameters.n_compositional_fields; ++c)
          {
            system_partitioning.push_back(system_index_set.get_view(n_u+n_p+n_T+n_C_so_far,
                                                                    n_u+n_p+n_T+n_C_so_far+n_C[c]));
            system_relevant_partitioning.push_back(system_relevant_set.get_view(n_u+n_p+n_T+n_C_so_far,
                                                                                n_u+n_p+n_T+n_C_so_far+n_C[c]));
            n_C_so_far += n_C[c];
          }
      }
    }

    // then compute constraints for the velocity. the constraints we compute
    // here are the ones that are the same for all following time steps. in
    // addition, we may be computing constraints from boundary values for the
    // velocity that are different between time steps. these are then put
    // into current_constraints in start_timestep().
    {
      constraints.clear();
      constraints.reinit(system_relevant_set);

      DoFTools::make_hanging_node_constraints (dof_handler,
                                               constraints);

      // do the interpolation for zero velocity
      std::vector<bool> velocity_mask (dim+2+parameters.n_compositional_fields, true);
      velocity_mask[dim] = false;
      velocity_mask[dim+1] = false;
      for (unsigned int i=dim+2; i<dim+2+parameters.n_compositional_fields; ++i)
        velocity_mask[i] = false;
      for (std::set<types::boundary_id_t>::const_iterator
           p = parameters.zero_velocity_boundary_indicators.begin();
           p != parameters.zero_velocity_boundary_indicators.end(); ++p)
        VectorTools::interpolate_boundary_values (dof_handler,
                                                  *p,
                                                  ZeroFunction<dim>(dim+2+parameters.n_compositional_fields),
                                                  constraints,
                                                  velocity_mask);


      // do the same for no-normal-flux boundaries
      VectorTools::compute_no_normal_flux_constraints (dof_handler,
                                                       /* first_vector_component= */ 0,
                                                       parameters.tangential_velocity_boundary_indicators,
                                                       constraints,
                                                       mapping);
    }

    // now do the same for the temperature variable
    {

      // obtain the boundary indicators that belong to Dirichlet-type
      // temperature boundary conditions and interpolate the temperature
      // there
      std::vector<bool> temperature_mask (dim+2+parameters.n_compositional_fields, false);
      temperature_mask[dim+1] = true;

      for (std::set<types::boundary_id_t>::const_iterator
           p = parameters.fixed_temperature_boundary_indicators.begin();
           p != parameters.fixed_temperature_boundary_indicators.end(); ++p)
        {
          Assert (is_element (*p, geometry_model->get_used_boundary_indicators()),
                  ExcInternalError());
          VectorTools::interpolate_boundary_values (dof_handler,
                                                    *p,
                                                    VectorFunctionFromScalarFunctionObject<dim>(std_cxx1x::bind (&BoundaryTemperature::Interface<dim>::temperature,
                                                        std_cxx1x::cref(*boundary_temperature),
                                                        std_cxx1x::cref(*geometry_model),
                                                        *p,
                                                        std_cxx1x::_1),
                                                        dim+1,
                                                        dim+2+parameters.n_compositional_fields),
                                                    constraints,
                                                    temperature_mask);

        }

      // we do nothing with the compositional fields: homogeneous Neumann boundary conditions

      constraints.close();
    }

    // finally initialize vectors, matrices, etc.

    setup_system_matrix (system_partitioning);
    setup_system_preconditioner (system_partitioning);

    system_rhs.reinit(system_partitioning, mpi_communicator);
    solution.reinit(system_relevant_partitioning, mpi_communicator);
    old_solution.reinit(system_relevant_partitioning, mpi_communicator);
    old_old_solution.reinit(system_relevant_partitioning, mpi_communicator);

    current_linearization_point.reinit (system_relevant_partitioning, MPI_COMM_WORLD);

    if (material_model->is_compressible())
      pressure_shape_function_integrals.reinit (system_partitioning, mpi_communicator);

    rebuild_stokes_matrix         = true;
    rebuild_stokes_preconditioner = true;

    computing_timer.exit_section();
  }


  template <int dim>
  void Simulator<dim>::setup_introspection ()
  {
    // the extractors for the default variables have already been initialized.
    // take care of the ones for the compositional variables now
    for (unsigned int c=0; c<parameters.n_compositional_fields; ++c)
      introspection.extractors.compositional_fields
      .push_back (FEValuesExtractors::Scalar(dim+2+c));

    // next initialize the component masks
    introspection.component_masks.velocities
    = finite_element.component_mask (introspection.extractors.velocities);
    introspection.component_masks.pressure
    = finite_element.component_mask (introspection.extractors.pressure);
    introspection.component_masks.temperature
    = finite_element.component_mask (introspection.extractors.temperature);
    for (unsigned int c=0; c<parameters.n_compositional_fields; ++c)
      introspection.component_masks.compositional_fields
      .push_back (finite_element.component_mask(introspection.extractors.compositional_fields[c]));
  }



  template <int dim>
  void Simulator<dim>::postprocess ()
  {
    computing_timer.enter_section ("Postprocessing");
    pcout << "   Postprocessing:" << std::endl;

    // run all the postprocessing routines and then write
    // the current state of the statistics table to a file
    std::list<std::pair<std::string,std::string> >
    output_list = postprocess_manager.execute (statistics);

    // if we are on processor zero, print to screen
    // whatever the postprocessors have generated
    if (Utilities::MPI::this_mpi_process(mpi_communicator)==0)
      {
        // determine the width of the first column of text so that
        // everything gets nicely aligned; then output everything
        {
          unsigned int width = 0;
          for (std::list<std::pair<std::string,std::string> >::const_iterator
               p = output_list.begin();
               p != output_list.end(); ++p)
            width = std::max<unsigned int> (width, p->first.size());

          for (std::list<std::pair<std::string,std::string> >::const_iterator
               p = output_list.begin();
               p != output_list.end(); ++p)
            pcout << "     "
                  << std::left
                  << std::setw(width)
                  << p->first
                  << " "
                  << p->second
                  << std::endl;
        }

        pcout << std::endl;
      }

    // finally, write the entire set of current results to disk
    output_statistics();

    computing_timer.exit_section ();
  }


// Contrary to step-32, we have found that just refining by the temperature
// works well in 2d, but only leads to refinement in the boundary layer at the
// core-mantle boundary in 3d. Consequently, we estimate the error based
// on the temperature, velocity and other criteria, see the second ASPECT paper;
// the vectors with the resulting error indicators are then normalized, and we
// take the maximum or sum of the indicators to decide whether we want to refine or
// not. In case of more complicated materials with jumps in the density
// profile, we also need to refine where the density jumps. This ensures that
// we also refine into plumes where maybe the temperature gradients aren't as
// strong as in the boundary layer but where nevertheless the gradients in the
// velocity are large.
  template <int dim>
  void Simulator<dim>::refine_mesh (const unsigned int max_grid_level)
  {
    computing_timer.enter_section ("Refine mesh structure, part 1");

    Vector<float> estimated_error_per_cell (triangulation.n_active_cells());
    mesh_refinement_manager.execute (estimated_error_per_cell);

    parallel::distributed::GridRefinement::
    refine_and_coarsen_fixed_fraction (triangulation,
                                       estimated_error_per_cell,
                                       parameters.refinement_fraction,
                                       parameters.coarsening_fraction);

    // limit maximum refinement level
    if (triangulation.n_levels() > max_grid_level)
      for (typename Triangulation<dim>::active_cell_iterator
           cell = triangulation.begin_active(max_grid_level);
           cell != triangulation.end(); ++cell)
        cell->clear_refine_flag ();

    std::vector<const LinearAlgebra::BlockVector *> x_system (2);
    x_system[0] = &solution;
    x_system[1] = &old_solution;

    parallel::distributed::SolutionTransfer<dim,LinearAlgebra::BlockVector>
    system_trans(dof_handler);

    triangulation.prepare_coarsening_and_refinement();
    system_trans.prepare_for_coarsening_and_refinement(x_system);

    triangulation.execute_coarsening_and_refinement ();
    global_volume = GridTools::volume (triangulation, mapping);
    computing_timer.exit_section();

    setup_dofs ();

    computing_timer.enter_section ("Refine mesh structure, part 2");

    {
      LinearAlgebra::BlockVector
      distributed_system (system_rhs);
      LinearAlgebra::BlockVector
      old_distributed_system (system_rhs);
      std::vector<LinearAlgebra::BlockVector *> system_tmp (2);
      system_tmp[0] = &(distributed_system);
      system_tmp[1] = &(old_distributed_system);

      system_trans.interpolate (system_tmp);
      solution     = distributed_system;
      old_solution = old_distributed_system;
    }

    computing_timer.exit_section();
  }


  template <int dim>
  void
  Simulator<dim>::
  solve_timestep ()
  {
    // start any scheme with an extrapolated value from the previous
    // two time steps if those are available
    current_linearization_point = old_solution;
    if (timestep_number > 1)
      {
        //TODO: Trilinos sadd does not like ghost vectors even as input. Copy
        //into distributed vectors for now:
        LinearAlgebra::BlockVector distr_solution (system_rhs);
        distr_solution = current_linearization_point;
        LinearAlgebra::BlockVector distr_old_solution (system_rhs);
        distr_old_solution = old_old_solution;
        distr_solution .sadd ((1 + time_step/old_time_step),
                              -time_step/old_time_step,
                              distr_old_solution);
        current_linearization_point = distr_solution;
      }

    switch (parameters.nonlinear_solver)
      {
        case NonlinearSolver::IMPES:
        {
          assemble_advection_system (0);
          build_advection_preconditioner(0,T_preconditioner);
          solve_advection(0);

          current_linearization_point.block(2) = solution.block(2);

          for (unsigned int c=0; c<parameters.n_compositional_fields; ++c)
            {
              assemble_advection_system (1+c);
              build_advection_preconditioner(1+c,C_preconditioner);
              solve_advection(1+c); // this is correct, 0 would be temperature
              current_linearization_point.block(3+c) = solution.block(3+c);
            }


          // the Stokes matrix depends on the viscosity. if the viscosity
          // depends on other solution variables, then after we need to
          // update the Stokes matrix in every time step and so need to set
          // the following flag. if we change the Stokes matrix we also
          // need to update the Stokes preconditioner.
          if (stokes_matrix_depends_on_solution() == true)
            rebuild_stokes_matrix = rebuild_stokes_preconditioner = true;

          assemble_stokes_system();
          build_stokes_preconditioner();
          solve_stokes();

          break;
        }

        case NonlinearSolver::iterated_IMPES:
        {
          double initial_temperature_residual = 0;
          double initial_stokes_residual      = 0;
          std::vector<double> initial_composition_residual (parameters.n_compositional_fields,0);

          unsigned int iteration = 0;

          do
            {
              assemble_advection_system(0);

              if (iteration == 0)
                build_advection_preconditioner(0,T_preconditioner);

              const double temperature_residual = solve_advection(0);

              current_linearization_point.block(2) = solution.block(2);
              rebuild_stokes_matrix = true;
              std::vector<double> composition_residual (parameters.n_compositional_fields,0);

              for (unsigned int c=0; c<parameters.n_compositional_fields; ++c)
                {
                  assemble_advection_system (c+1);
                  build_advection_preconditioner(c+1,C_preconditioner);
                  composition_residual[c]
                    = solve_advection(1+c); // 1+n is correct, because 0 is for temperature
                  current_linearization_point.block(3+c) = solution.block(3+c);
                }

              // the Stokes matrix depends on the viscosity. if the viscosity
              // depends on other solution variables, then after we need to
              // update the Stokes matrix in every time step and so need to set
              // the following flag. if we change the Stokes matrix we also
              // need to update the Stokes preconditioner.
              if (stokes_matrix_depends_on_solution() == true)
                rebuild_stokes_matrix = rebuild_stokes_preconditioner = true;

              assemble_stokes_system();
              if (iteration == 0)
                build_stokes_preconditioner();

              const double stokes_residual = solve_stokes();

              current_linearization_point = solution;

              pcout << "   Nonlinear residuals: " << temperature_residual
                    << ", " << stokes_residual;

              for (unsigned int c=0; c<parameters.n_compositional_fields; ++c)
                pcout << ", " << composition_residual[c];

              pcout << std::endl
                    << std::endl;

              if (iteration == 0)
                {
                  initial_temperature_residual = temperature_residual;
                  for (unsigned int c=0; c<parameters.n_compositional_fields; ++c)
                    initial_composition_residual[c] = composition_residual[c];
                  initial_stokes_residual      = stokes_residual;
                }
              else
                {
//TODO: make this a parameter in the input file
                  double max = 0.0;
                  for (unsigned int c=0; c<parameters.n_compositional_fields; ++c)
                    max = std::max(composition_residual[c]/initial_composition_residual[c],max);
                  if (std::max(std::max(stokes_residual/initial_stokes_residual,
                                        temperature_residual/initial_temperature_residual),
                               max) < 1e-3)
                    break;
                }

              ++iteration;
//TODO: terminate here if the number of iterations is too large and we see no convergence
            }
          while (iteration < 10);

          break;
        }

        case NonlinearSolver::iterated_Stokes:
        {
          // solve the temperature system once...
          assemble_advection_system (0);
          solve_advection(0);

          for (unsigned int c=0; c<parameters.n_compositional_fields; ++c)
            {
              assemble_advection_system (c+1);
              solve_advection(1+c); // 1+n is correct, because 0 is the temperature
            }

          // ...and then iterate the solution
          // of the Stokes system
          for (int i=0; i<10; ++i)
            {
              // rebuild the matrix if it actually depends on the solution
              // of the previous iteration.
              if (stokes_matrix_depends_on_solution() == true)
                rebuild_stokes_matrix = rebuild_stokes_preconditioner = true;

              assemble_stokes_system();
              build_stokes_preconditioner();
              solve_stokes();
              old_solution = solution;

//TODO: don't we need to set the linearization point here somehow?
              Assert (false, ExcNotImplemented());

              pcout << std::endl;
            }

          break;
        }

        default:
          Assert (false, ExcNotImplemented());
      }
  }


  /**
   * This is the main function of the program, containing the overall
   * logic which function is called when.
   */
  template <int dim>
  void Simulator<dim>::run ()
  {
    last_checkpoint_time = std::time(NULL);
    unsigned int max_refinement_level = parameters.initial_global_refinement +
                                        parameters.initial_adaptive_refinement;
    unsigned int pre_refinement_step = 0;

    // if we want to resume a computation from an earlier point
    // then reload it from a snapshot. otherwise do the basic
    // start-up
    if (parameters.resume_computation == true)
      {
        resume_from_snapshot();
        // we need to remove additional_refinement_times that are in the past
        // and adjust max_refinement_level which is not written to file
        while ((parameters.additional_refinement_times.size() > 0)
               &&
               (parameters.additional_refinement_times.front () < time+time_step))
          {
            ++max_refinement_level;
            parameters.additional_refinement_times
            .erase (parameters.additional_refinement_times.begin());
          }
      }
    else
      {
        triangulation.refine_global (parameters.initial_global_refinement);
        global_volume = GridTools::volume (triangulation, mapping);

        setup_dofs();
      }


  start_time_iteration:

    if (parameters.resume_computation == false)
      {
        computing_timer.enter_section ("Initialization");
        set_initial_temperature_and_compositional_fields ();
        compute_initial_pressure_field ();

        time                      = parameters.start_time;
        timestep_number           = 0;
        time_step = old_time_step = 0;
        computing_timer.exit_section();
      }

    // start the principal loop over time steps:
    do
      {
        start_timestep ();

        // then do the core work: assemble systems and solve
        solve_timestep ();

        pcout << std::endl;

        // update the time step size
        // for now the bool (convection/conduction dominated)
        // is unused, will be added to statistics later
        std::pair<double, bool>   time_step_result;
        old_time_step = time_step;
        time_step_result = compute_time_step();
        time_step = time_step_result.first;

        if (parameters.convert_to_years == true)
          statistics.add_value("Time step size (years)", time_step / year_in_seconds);
        else
          statistics.add_value("Time step size (seconds)", time_step);


        // see if we have to start over with a new refinement cycle
        // at the beginning of the simulation
        if ((timestep_number == 0) &&
            (pre_refinement_step < parameters.initial_adaptive_refinement))
          {
            output_statistics();

            if (parameters.run_postprocessors_on_initial_refinement)
              postprocess ();

            refine_mesh (max_refinement_level);
            ++pre_refinement_step;
            goto start_time_iteration;
          }

        postprocess ();

        // see if this is a time step where additional refinement is requested
        // if so, then loop over as many times as this is necessary
        if ((parameters.additional_refinement_times.size() > 0)
            &&
            (parameters.additional_refinement_times.front () < time+time_step))
          {
            while ((parameters.additional_refinement_times.size() > 0)
                   &&
                   (parameters.additional_refinement_times.front () < time+time_step))
              {
                ++max_refinement_level;
                refine_mesh (max_refinement_level);

                parameters.additional_refinement_times
                .erase (parameters.additional_refinement_times.begin());
              }
          }
        else
          // see if this is a time step where regular refinement is necessary, but only
          // if the previous rule wasn't triggered
          if ((timestep_number > 0)
              &&
              (parameters.adaptive_refinement_interval > 0)
              &&
              (timestep_number % parameters.adaptive_refinement_interval == 0))
            refine_mesh (max_refinement_level);


        // every n time steps output a summary of the current
        // timing information
        if ((timestep_number > 0) && (parameters.timing_output_frequency != 0) &&
            (timestep_number % parameters.timing_output_frequency == 0))
          {
            computing_timer.print_summary ();
            output_statistics();
          }

        // increment time step by one. then prepare
        // for the next time step by shifting solution vectors
        // by one time step
        time += time_step;
        ++timestep_number;
        {
          old_old_solution      = old_solution;
          old_solution          = solution;
        }

        // periodically generate snapshots so that we can resume here
        // if the program aborts or is terminated
        bool do_checkpoint = false;

        // If we base checkpoint frequency on timing, measure the time at process 0
        // This prevents race conditions where some processes will checkpoint and others won't
        if (parameters.checkpoint_time_secs > 0)
          {
            int global_do_checkpoint = ((std::time(NULL)-last_checkpoint_time) >= parameters.checkpoint_time_secs);

            MPI_Bcast(&global_do_checkpoint, 1, MPI_INT, 0, MPI_COMM_WORLD);

            do_checkpoint = (global_do_checkpoint == 1);
          }
        // If we base checkpoint frequency on steps, see if it's time for another checkpoint
        if ((parameters.checkpoint_time_secs == 0) &&
            (parameters.checkpoint_steps > 0) &&
            (timestep_number % parameters.checkpoint_steps == 0))
          do_checkpoint = true;
        if (do_checkpoint)
          {
            create_snapshot();
            // matrices will be regenerated after a resume, so do that here too
            // to be consistent. otherwise we would get different results
            // for a restarted computation than for one that ran straight
            // through
            rebuild_stokes_matrix =
              rebuild_stokes_preconditioner = true;
            last_checkpoint_time = std::time(NULL);
          }
      }
    while (time < parameters.end_time);
  }
}



// explicit instantiation of the functions we implement in this file
namespace aspect
{
#define INSTANTIATE(dim) \
  template class Simulator<dim>;

  ASPECT_INSTANTIATE(INSTANTIATE)
}
