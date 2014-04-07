/*
  Copyright (C) 2011, 2012, 2013 by the authors of the ASPECT code.

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
#include <deal.II/lac/sparsity_tools.h>
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
    introspection (parameters.n_compositional_fields),
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
    // create a boundary composition model, but only if we actually need
    // it. otherwise, allow the user to simply specify nothing at all
    boundary_composition (parameters.fixed_composition_boundary_indicators.empty()
			  ?
			  0
			  :
			  BoundaryComposition::create_boundary_composition<dim>(prm)),
    compositional_initial_conditions (CompositionalInitialConditions::create_initial_conditions (prm,
                                      *geometry_model)),
    adiabatic_conditions(),
    initial_conditions (),

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
      std::set<types::boundary_id> boundary_indicator_lists[3]
        = { parameters.zero_velocity_boundary_indicators,
            parameters.tangential_velocity_boundary_indicators,
            std::set<types::boundary_id>()
          };

      for (std::map<types::boundary_id,std::pair<std::string, std::string> >::const_iterator
           p = parameters.prescribed_velocity_boundary_indicators.begin();
           p != parameters.prescribed_velocity_boundary_indicators.end();
           ++p)
        boundary_indicator_lists[2].insert (p->first);

      // for each combination of boundary indicator lists, make sure that the
      // intersection is empty
      for (unsigned int i=0; i<sizeof(boundary_indicator_lists)/sizeof(boundary_indicator_lists[0]); ++i)
        for (unsigned int j=i+1; j<sizeof(boundary_indicator_lists)/sizeof(boundary_indicator_lists[0]); ++j)
          {
            std::set<types::boundary_id> intersection;
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
      const std::set<types::boundary_id> all_boundary_indicators
        = geometry_model->get_used_boundary_indicators();
      for (unsigned int i=0; i<sizeof(boundary_indicator_lists)/sizeof(boundary_indicator_lists[0]); ++i)
        for (typename std::set<types::boundary_id>::const_iterator
             p = boundary_indicator_lists[i].begin();
             p != boundary_indicator_lists[i].end(); ++p)
          AssertThrow (all_boundary_indicators.find (*p)
                       != all_boundary_indicators.end(),
                       ExcMessage ("One of the boundary indicators listed in the input file "
                                   "is not used by the geometry model."));

      // now do the same for the fixed temperature indicators and the
      // compositional indicators
      for (typename std::set<types::boundary_id>::const_iterator
           p = parameters.fixed_temperature_boundary_indicators.begin();
           p != parameters.fixed_temperature_boundary_indicators.end(); ++p)
        AssertThrow (all_boundary_indicators.find (*p)
                     != all_boundary_indicators.end(),
                     ExcMessage ("One of the fixed boundary temperature indicators listed in the input file "
                                 "is not used by the geometry model."));
      for (typename std::set<types::boundary_id>::const_iterator
           p = parameters.fixed_composition_boundary_indicators.begin();
           p != parameters.fixed_composition_boundary_indicators.end(); ++p)
        AssertThrow (all_boundary_indicators.find (*p)
                     != all_boundary_indicators.end(),
                     ExcMessage ("One of the fixed boundary composition indicators listed in the input file "
                                 "is not used by the geometry model."));
    }

    // continue with initializing members that can't be initialized for one reason
    // or another in the member initializer list above

    // if any plugin wants access to the Simulator by deriving from SimulatorAccess, initialize it:
    if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(geometry_model.get()))
      sim->initialize (*this);
    if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(material_model.get()))
      sim->initialize (*this);
    if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(gravity_model.get()))
      sim->initialize (*this);
    if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(boundary_temperature.get()))
      sim->initialize (*this);
    if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(boundary_composition.get()))
      sim->initialize (*this);

    adiabatic_conditions.reset(new AdiabaticConditions<dim> (*geometry_model,
                                                             *gravity_model,
                                                             *material_model,
                                                             *compositional_initial_conditions,
                                                             parameters.surface_pressure,
                                                             parameters.adiabatic_surface_temperature,
                                                             parameters.n_compositional_fields));

    initial_conditions.reset (InitialConditions::create_initial_conditions (prm,
                                                                            *geometry_model,
                                                                            *boundary_temperature,
                                                                            *adiabatic_conditions));
    if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(initial_conditions.get()))
      sim->initialize (*this);

    postprocess_manager.parse_parameters (prm);
    postprocess_manager.initialize (*this);

    mesh_refinement_manager.parse_parameters (prm);
    mesh_refinement_manager.initialize (*this);

    termination_manager.parse_parameters (prm);
    termination_manager.initialize (*this);

    geometry_model->create_coarse_mesh (triangulation);
    global_Omega_diameter = GridTools::diameter (triangulation);

    for (std::map<types::boundary_id,std::pair<std::string,std::string> >::const_iterator
         p = parameters.prescribed_velocity_boundary_indicators.begin();
         p != parameters.prescribed_velocity_boundary_indicators.end();
         ++p)
      {
        VelocityBoundaryConditions::Interface<dim> *bv
          = VelocityBoundaryConditions::create_velocity_boundary_conditions
            (p->second.second,
             prm,
             *geometry_model);
        if (dynamic_cast<SimulatorAccess<dim>*>(bv) != 0)
          dynamic_cast<SimulatorAccess<dim>*>(bv)->initialize(*this);
        velocity_boundary_conditions[p->first].reset (bv);
      }

    // determine how to treat the pressure. we have to scale it for the solver
    // to make velocities and pressures of roughly the same (numerical) size,
    // and we may have to fix up the right hand side vector before solving for
    // compressible models if there are no in-/outflow boundaries
    pressure_scaling = material_model->reference_viscosity() / geometry_model->length_scale();

    std::set<types::boundary_id> open_velocity_boundary_indicators
    = geometry_model->get_used_boundary_indicators();
    for (std::map<types::boundary_id,std::pair<std::string,std::string> >::const_iterator
         p = parameters.prescribed_velocity_boundary_indicators.begin();
         p != parameters.prescribed_velocity_boundary_indicators.end();
         ++p)
      open_velocity_boundary_indicators.erase (p->first);
    for (std::set<types::boundary_id>::const_iterator
         p = parameters.zero_velocity_boundary_indicators.begin();
         p != parameters.zero_velocity_boundary_indicators.end();
         ++p)
      open_velocity_boundary_indicators.erase (*p);
    for (std::set<types::boundary_id>::const_iterator
         p = parameters.tangential_velocity_boundary_indicators.begin();
         p != parameters.tangential_velocity_boundary_indicators.end();
         ++p)
      open_velocity_boundary_indicators.erase (*p);
//TODO: The "correct" condition would be to not do the correction for compressible
    // models if there are either open boundaries, or if there are boundaries
    // where the prescribed velocity is so that in- and outflux do not exactly
    // match
    do_pressure_rhs_compatibility_modification = (material_model->is_compressible()
                                                  &&
                                                  (parameters.prescribed_velocity_boundary_indicators.size() == 0)
                                                  &&
                                                  (open_velocity_boundary_indicators.size() == 0));

    // make sure that we don't have to fill every column of the statistics
    // object in each time step.
    statistics.set_auto_fill_mode(true);

    // finally produce a record of the run-time parameters by writing
    // the currently used values into a file
    // Only write the parameter files on the root node to avoid file system conflicts
    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
      {
        std::ofstream prm_out ((parameters.output_directory + "parameters.prm").c_str());
        AssertThrow (prm_out,
                     ExcMessage (std::string("Couldn't open file <") +
                                 parameters.output_directory + "parameters.prm>."));
        prm.print_parameters(prm_out, ParameterHandler::Text);
      }
    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
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
   */
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

    statistics.add_value("Number of Stokes degrees of freedom",
                         introspection.system_dofs_per_block[0] +
                         introspection.system_dofs_per_block[1]);
    statistics.add_value("Number of temperature degrees of freedom",
                         introspection.system_dofs_per_block[2]);
    if (parameters.n_compositional_fields > 0)
      statistics.add_value("Number of composition degrees of freedom",
                           introspection.system_dofs_per_block[3]);


    // then interpolate the current boundary velocities. this adds to
    // the current_constraints object we already have
    {
      current_constraints.clear ();
      current_constraints.reinit (introspection.index_sets.system_relevant_set);
      current_constraints.merge (constraints);

      // set the current time and do the interpolation
      // for the prescribed velocity fields
      for (typename std::map<types::boundary_id,std_cxx1x::shared_ptr<VelocityBoundaryConditions::Interface<dim> > >::iterator
           p = velocity_boundary_conditions.begin();
           p != velocity_boundary_conditions.end(); ++p)
        {
          p->second->set_current_time (time);
          VectorFunctionFromVelocityFunctionObject<dim> vel
          (parameters.n_compositional_fields,
           std_cxx1x::bind (&VelocityBoundaryConditions::Interface<dim>::boundary_velocity,
                            p->second,
                            std_cxx1x::_1));

          // here we create a mask for interpolate_boundary_values out of the 'selector'
          std::vector<bool> mask(introspection.component_masks.velocities.size(), false);
          Assert(introspection.component_masks.velocities[0]==true, ExcInternalError()); // in case we ever move the velocity around
          const std::string & comp = parameters.prescribed_velocity_boundary_indicators[p->first].first;

          if (comp.length()>0)
            {
              for (std::string::const_iterator direction=comp.begin();direction!=comp.end();++direction)
                {
                  AssertThrow(*direction>='x' && *direction<='z', ExcMessage("Error in selector of prescribed velocity boundary component"));
                  AssertThrow(dim==3 || *direction!='z', ExcMessage("for dim=2, prescribed velocity component z is invalid"))
                  mask[*direction-'x']=true;
                }
              for (unsigned int i=0;i<introspection.component_masks.velocities.size();++i)
                mask[i] = mask[i] & introspection.component_masks.velocities[i];
            }
          else
            {
              for (unsigned int i=0;i<introspection.component_masks.velocities.size();++i)
                  mask[i]=introspection.component_masks.velocities[i];

              Assert(introspection.component_masks.velocities[0]==true, ExcInternalError()); // in case we ever move the velocity down
            }

          VectorTools::interpolate_boundary_values (dof_handler,
                                                    p->first,
                                                    vel,
                                                    current_constraints,
                                                    mask);
        }
      current_constraints.close();
    }

    //TODO: do this in a more efficient way (TH)? we really only need
    // to make sure that the time dependent velocity boundary conditions
    // end up in the right hand side in the right way; we currently do
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

    Table<2,DoFTools::Coupling> coupling (introspection.n_components,
                                          introspection.n_components);

    // determine which blocks should be fillable in the matrix.
    // note:
    // - all velocities couple with all velocities
    // - pressure couples with all velocities and the other way
    //   around
    // - temperature only couples with itself when using the impes
    //   scheme
    // - compositional fields only couple with themselves
    {
      const typename Introspection<dim>::ComponentIndices &x
        = introspection.component_indices;

      for (unsigned int c=0; c<dim; ++c)
        for (unsigned int d=0; d<dim; ++d)
          coupling[x.velocities[c]][x.velocities[d]] = DoFTools::always;
      for (unsigned int d=0; d<dim; ++d)
        {
          coupling[x.velocities[d]][x.pressure] = DoFTools::always;
          coupling[x.pressure][x.velocities[d]] = DoFTools::always;
        }
      coupling[x.temperature][x.temperature] = DoFTools::always;
      for (unsigned int c=0; c<parameters.n_compositional_fields; ++c)
        coupling[x.compositional_fields[c]][x.compositional_fields[c]]
          = DoFTools::always;
    }

#ifdef USE_PETSC
    LinearAlgebra::CompressedBlockSparsityPattern sp(introspection.index_sets.system_relevant_partitioning);

#else
    TrilinosWrappers::BlockSparsityPattern sp (system_partitioning,
                                               mpi_communicator);
#endif

    DoFTools::make_sparsity_pattern (dof_handler,
                                     coupling, sp,
                                     constraints, false,
                                     Utilities::MPI::
                                     this_mpi_process(mpi_communicator));

#ifdef USE_PETSC
    SparsityTools::distribute_sparsity_pattern(sp,
        dof_handler.locally_owned_dofs_per_processor(),
        mpi_communicator, introspection.index_sets.system_relevant_set);

    sp.compress();

    system_matrix.reinit (system_partitioning, system_partitioning, sp, mpi_communicator);
#else
    sp.compress();

    system_matrix.reinit (sp);
#endif
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

    Table<2,DoFTools::Coupling> coupling (introspection.n_components,
                                          introspection.n_components);
    for (unsigned int c=0; c<introspection.n_components; ++c)
      for (unsigned int d=0; d<introspection.n_components; ++d)
        if (c == d)
          coupling[c][d] = DoFTools::always;
        else
          coupling[c][d] = DoFTools::none;


#ifdef USE_PETSC
    LinearAlgebra::CompressedBlockSparsityPattern sp(introspection.index_sets.system_relevant_partitioning);

#else
    TrilinosWrappers::BlockSparsityPattern sp (system_partitioning,
                                               mpi_communicator);
#endif
    DoFTools::make_sparsity_pattern (dof_handler,
                                     coupling, sp,
                                     constraints, false,
                                     Utilities::MPI::
                                     this_mpi_process(mpi_communicator));
#ifdef USE_PETSC
    SparsityTools::distribute_sparsity_pattern(sp,
        dof_handler.locally_owned_dofs_per_processor(),
        mpi_communicator, introspection.index_sets.system_relevant_set);

    sp.compress();

    system_preconditioner_matrix.reinit (system_partitioning, system_partitioning, sp, mpi_communicator);
#else
    sp.compress();

    system_preconditioner_matrix.reinit (sp);
#endif
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
    DoFRenumbering::component_wise (dof_handler,
                                    introspection.components_to_blocks);

    // set up the introspection object that stores all sorts of
    // information about components of the finite element, component
    // masks, etc
    setup_introspection();

    // print dof numbers. Do so with 1000s separator since they are frequently
    // large
    {
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
            << dof_handler.n_dofs()
            << " ("
            << introspection.system_dofs_per_block[0];
      for (unsigned int b=1; b<introspection.system_dofs_per_block.size(); ++b)
        pcout << '+' << introspection.system_dofs_per_block[b];
      pcout <<')'
            << std::endl
            << std::endl;

      pcout.get_stream().imbue(s);
    }

    //reinit the constraints matrix and make hanging node constraints
    constraints.clear();
    constraints.reinit(introspection.index_sets.system_relevant_set);
    DoFTools::make_hanging_node_constraints (dof_handler,
                                             constraints);

    //Now set up the constraints for periodic boundary conditions
    {
      typedef std::set< std::pair< std::pair< types::boundary_id, types::boundary_id>, unsigned int> >
               periodic_boundary_set;
      periodic_boundary_set pbs = geometry_model->get_periodic_boundary_pairs();

      for(periodic_boundary_set::iterator p = pbs.begin(); p != pbs.end(); ++p)
        {
          //Throw error if we are trying to use the same boundary for more than one boundary condition
          Assert( is_element( (*p).first.first, parameters.fixed_temperature_boundary_indicators ) == false &&
                  is_element( (*p).first.second, parameters.fixed_temperature_boundary_indicators ) == false &&
                  is_element( (*p).first.first, parameters.zero_velocity_boundary_indicators ) == false &&
                  is_element( (*p).first.second, parameters.zero_velocity_boundary_indicators ) == false &&
                  is_element( (*p).first.first, parameters.tangential_velocity_boundary_indicators ) == false &&
                  is_element( (*p).first.second, parameters.tangential_velocity_boundary_indicators ) == false &&
                  parameters.prescribed_velocity_boundary_indicators.find( (*p).first.first)
                             == parameters.prescribed_velocity_boundary_indicators.end() &&
                  parameters.prescribed_velocity_boundary_indicators.find( (*p).first.second)
                             == parameters.prescribed_velocity_boundary_indicators.end(),
                  ExcInternalError());

#if (DEAL_II_MAJOR*100 + DEAL_II_MINOR) >= 801
          DoFTools::make_periodicity_constraints(dof_handler,
                                                 (*p).first.first,  //first boundary id
                                                 (*p).first.second, //second boundary id
                                                 (*p).second,       //cartesian direction for translational symmetry
                                                 constraints);
#endif
        }


    }
    // then compute constraints for the velocity. the constraints we compute
    // here are the ones that are the same for all following time steps. in
    // addition, we may be computing constraints from boundary values for the
    // velocity that are different between time steps. these are then put
    // into current_constraints in start_timestep().
    {


      // do the interpolation for zero velocity
      for (std::set<types::boundary_id>::const_iterator
           p = parameters.zero_velocity_boundary_indicators.begin();
           p != parameters.zero_velocity_boundary_indicators.end(); ++p)
        VectorTools::interpolate_boundary_values (dof_handler,
                                                  *p,
                                                  ZeroFunction<dim>(introspection.n_components),
                                                  constraints,
                                                  introspection.component_masks.velocities);


      // do the same for no-normal-flux boundaries
      VectorTools::compute_no_normal_flux_constraints (dof_handler,
                                                       /* first_vector_component= */
                                                       introspection.component_indices.velocities[0],
                                                       parameters.tangential_velocity_boundary_indicators,
                                                       constraints,
                                                       mapping);
    }

    // now do the same for the temperature variable
    {
      // obtain the boundary indicators that belong to Dirichlet-type
      // temperature boundary conditions and interpolate the temperature
      // there
      for (std::set<types::boundary_id>::const_iterator
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
                                                        introspection.component_masks.temperature.first_selected_component(),
                                                        introspection.n_components),
                                                    constraints,
                                                    introspection.component_masks.temperature);

        }

      // now do the same for the composition variable:
      // obtain the boundary indicators that belong to Dirichlet-type
      // composition boundary conditions and interpolate the composition
      // there
      for (unsigned int c=0; c<parameters.n_compositional_fields; ++c)
        for (std::set<types::boundary_id>::const_iterator
             p = parameters.fixed_composition_boundary_indicators.begin();
             p != parameters.fixed_composition_boundary_indicators.end(); ++p)
          {
            Assert (is_element (*p, geometry_model->get_used_boundary_indicators()),
                    ExcInternalError());
            VectorTools::interpolate_boundary_values (dof_handler,
                                                      *p,
                                                      VectorFunctionFromScalarFunctionObject<dim>(std_cxx1x::bind (&BoundaryComposition::Interface<dim>::composition,
                                                          std_cxx1x::cref(*boundary_composition),
                                                          std_cxx1x::cref(*geometry_model),
                                                          *p,
                                                          std_cxx1x::_1,
                                                          c),
                                                          introspection.component_masks.compositional_fields[c].first_selected_component(),
                                                          introspection.n_components),
                                                      constraints,
                                                      introspection.component_masks.compositional_fields[c]);
	  }
    }
    constraints.close();

    // finally initialize vectors, matrices, etc.

    setup_system_matrix (introspection.index_sets.system_partitioning);
    setup_system_preconditioner (introspection.index_sets.system_partitioning);

    system_rhs.reinit(introspection.index_sets.system_partitioning, mpi_communicator);
    solution.reinit(introspection.index_sets.system_partitioning, introspection.index_sets.system_relevant_partitioning, mpi_communicator);
    old_solution.reinit(introspection.index_sets.system_partitioning, introspection.index_sets.system_relevant_partitioning, mpi_communicator);
    old_old_solution.reinit(introspection.index_sets.system_partitioning, introspection.index_sets.system_relevant_partitioning, mpi_communicator);
    current_linearization_point.reinit (introspection.index_sets.system_partitioning, introspection.index_sets.system_relevant_partitioning, mpi_communicator);

    if (do_pressure_rhs_compatibility_modification)
      pressure_shape_function_integrals.reinit (introspection.index_sets.system_partitioning, mpi_communicator);

    rebuild_stokes_matrix         = true;
    rebuild_stokes_preconditioner = true;

    setup_nullspace_removal();

    computing_timer.exit_section();
  }


  template <int dim>
  void Simulator<dim>::setup_introspection ()
  {
    // initialize the component masks, if not already set
    if (introspection.component_masks.velocities == ComponentMask())
      {
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


    // now also compute the various partitionings between processors and blocks
    // of vectors and matrices
    DoFTools::count_dofs_per_block (dof_handler,
                                    introspection.system_dofs_per_block,
                                    introspection.components_to_blocks);
    {
      const types::global_dof_index n_u = introspection.system_dofs_per_block[0],
                                    n_p = introspection.system_dofs_per_block[1],
                                    n_T = introspection.system_dofs_per_block[2];
      std::vector<types::global_dof_index> n_C (parameters.n_compositional_fields+1);
      for (unsigned int c=0; c<parameters.n_compositional_fields; ++c)
        n_C[c] = introspection.system_dofs_per_block
                 [introspection.block_indices.compositional_fields[c]];


      IndexSet system_index_set = dof_handler.locally_owned_dofs();
      introspection.index_sets.system_partitioning.clear ();
      introspection.index_sets.system_partitioning.push_back(system_index_set.get_view(0,n_u));
      introspection.index_sets.system_partitioning.push_back(system_index_set.get_view(n_u,n_u+n_p));
      introspection.index_sets.system_partitioning.push_back(system_index_set.get_view(n_u+n_p,n_u+n_p+n_T));
      introspection.index_sets.stokes_partitioning.clear ();
      introspection.index_sets.stokes_partitioning.push_back(system_index_set.get_view(0,n_u));
      introspection.index_sets.stokes_partitioning.push_back(system_index_set.get_view(n_u,n_u+n_p));

      DoFTools::extract_locally_relevant_dofs (dof_handler,
                                               introspection.index_sets.system_relevant_set);
      introspection.index_sets.system_relevant_partitioning.clear ();
      introspection.index_sets.system_relevant_partitioning
      .push_back(introspection.index_sets.system_relevant_set.get_view(0,n_u));
      introspection.index_sets.system_relevant_partitioning
      .push_back(introspection.index_sets.system_relevant_set.get_view(n_u,n_u+n_p));
      introspection.index_sets.system_relevant_partitioning
      .push_back(introspection.index_sets.system_relevant_set.get_view(n_u+n_p, n_u+n_p+n_T));

      {
        unsigned int n_C_so_far = 0;

        for (unsigned int c=0; c<parameters.n_compositional_fields; ++c)
          {
            introspection.index_sets.system_partitioning.push_back(system_index_set.get_view(n_u+n_p+n_T+n_C_so_far,
                                                                   n_u+n_p+n_T+n_C_so_far+n_C[c]));
            introspection.index_sets.system_relevant_partitioning
            .push_back(introspection.index_sets.system_relevant_set.get_view(n_u+n_p+n_T+n_C_so_far,
                                                                             n_u+n_p+n_T+n_C_so_far+n_C[c]));
            n_C_so_far += n_C[c];
          }
      }
    }

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

    // limit minimum refinement level
    for (typename Triangulation<dim>::active_cell_iterator
         cell = triangulation.begin_active(0);
         cell != triangulation.end_active(parameters.min_grid_level); ++cell)
      cell->clear_coarsen_flag ();

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
        distr_solution = old_solution;
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
          assemble_advection_system (TemperatureOrComposition::temperature());
          build_advection_preconditioner(TemperatureOrComposition::temperature(),
                                         T_preconditioner);
          solve_advection(TemperatureOrComposition::temperature());

          current_linearization_point.block(introspection.block_indices.temperature)
            = solution.block(introspection.block_indices.temperature);

          for (unsigned int c=0; c<parameters.n_compositional_fields; ++c)
            {
              assemble_advection_system (TemperatureOrComposition::composition(c));
              build_advection_preconditioner(TemperatureOrComposition::composition(c),
                                             C_preconditioner);
              solve_advection(TemperatureOrComposition::composition(c)); // this is correct, 0 would be temperature
              current_linearization_point.block(introspection.block_indices.compositional_fields[c])
                = solution.block(introspection.block_indices.compositional_fields[c]);
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
        case NonlinearSolver::Stokes_only:
        {
          unsigned int iteration = 0;

          do
            {
              // the Stokes matrix depends on the viscosity. if the viscosity
              // depends on other solution variables, then we need to
              // update the Stokes matrix in every iteration and so need to set
              // the rebuild_stokes_matrix flag. if we change the Stokes matrix we also
              // need to update the Stokes preconditioner.
              //
              // there is a similar case where this nonlinear solver can be used, namely for
              // compressible models. in that case, the matrix does not depend on
              // the previous solution, but we still need to iterate since the right
              // hand side depends on it. in those cases, the matrix does not change,
              // but if we have to repeat computing the right hand side, we need to
              // also rebuild the matrix if we end up with inhomogenous velocity
              // boundary conditions (i.e., if there are prescribed velocity boundary
              // indicators)
              if ((stokes_matrix_depends_on_solution() == true)
                  ||
                  (parameters.prescribed_velocity_boundary_indicators.size() > 0))
                rebuild_stokes_matrix = true;
              if (stokes_matrix_depends_on_solution() == true)
                rebuild_stokes_preconditioner = true;

              assemble_stokes_system();
              build_stokes_preconditioner();
              const double stokes_residual = solve_stokes();
              current_linearization_point = solution;

              pcout << "      Nonlinear Stokes residual: " << stokes_residual << std::endl;
              if (stokes_residual < 1e-8)
                break;

              ++iteration;
            }
          while (iteration < parameters.max_nonlinear_iterations);
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
              assemble_advection_system(TemperatureOrComposition::temperature());

              if (iteration == 0)
                build_advection_preconditioner(TemperatureOrComposition::temperature(),
                                               T_preconditioner);

              const double temperature_residual = solve_advection(TemperatureOrComposition::temperature());

              current_linearization_point.block(introspection.block_indices.temperature)
                = solution.block(introspection.block_indices.temperature);
              rebuild_stokes_matrix = true;
              std::vector<double> composition_residual (parameters.n_compositional_fields,0);

              for (unsigned int c=0; c<parameters.n_compositional_fields; ++c)
                {
                  assemble_advection_system (TemperatureOrComposition::composition(c));
                  build_advection_preconditioner(TemperatureOrComposition::composition(c),
                                                 C_preconditioner);
                  composition_residual[c]
                    = solve_advection(TemperatureOrComposition::composition(c));
                  current_linearization_point.block(introspection.block_indices.compositional_fields[c])
                    = solution.block(introspection.block_indices.compositional_fields[c]);
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

              pcout << "      Nonlinear residuals: " << temperature_residual
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
                  double max = 0.0;
                  for (unsigned int c=0; c<parameters.n_compositional_fields; ++c)
                    max = std::max(composition_residual[c]/initial_composition_residual[c],max);
                  max = std::max(stokes_residual/initial_stokes_residual, max);
                  max = std::max(temperature_residual/initial_temperature_residual, max);
                  pcout << "      residual: " << max << std::endl;
                  if (max < parameters.nonlinear_tolerance)
                    break;
                }

              ++iteration;
//TODO: terminate here if the number of iterations is too large and we see no convergence
            }
          while (iteration < parameters.max_nonlinear_iterations);

          break;
        }

        case NonlinearSolver::iterated_Stokes:
        {
          // solve the temperature system once...
          assemble_advection_system (TemperatureOrComposition::temperature());
          build_advection_preconditioner (TemperatureOrComposition::temperature (),
                                          T_preconditioner);
          solve_advection(TemperatureOrComposition::temperature());
          current_linearization_point.block(introspection.block_indices.temperature)
            = solution.block(introspection.block_indices.temperature);

          for (unsigned int c=0; c<parameters.n_compositional_fields; ++c)
            {
              assemble_advection_system (TemperatureOrComposition::composition(c));
              build_advection_preconditioner (TemperatureOrComposition::composition (c),
                                              C_preconditioner);
              solve_advection(TemperatureOrComposition::composition(c));
              current_linearization_point.block(introspection.block_indices.compositional_fields[c])
                = solution.block(introspection.block_indices.compositional_fields[c]);
            }


          // residual vector (only for the velocity)
          LinearAlgebra::Vector residual (introspection.index_sets.system_partitioning[0], mpi_communicator);
          LinearAlgebra::Vector tmp (introspection.index_sets.system_partitioning[0], mpi_communicator);

          // ...and then iterate the solution of the Stokes system
          double initial_stokes_residual = 0;
          for (unsigned int i=0; i< parameters.max_nonlinear_iterations; ++i)
            {
              // rebuild the matrix if it actually depends on the solution
              // of the previous iteration.
              if ((stokes_matrix_depends_on_solution() == true)
                  ||
                  (parameters.prescribed_velocity_boundary_indicators.size() > 0))
                rebuild_stokes_matrix = rebuild_stokes_preconditioner = true;

              assemble_stokes_system();
              build_stokes_preconditioner();
              const double stokes_residual = solve_stokes();

              if (i==0)
                  initial_stokes_residual = stokes_residual;
              else
                {
                  pcout << "      residual: " << stokes_residual/initial_stokes_residual << std::endl;
                  if (stokes_residual/initial_stokes_residual < parameters.nonlinear_tolerance)
                    {
                      break; // convergence reached, exist nonlinear iteration.
                    }
                }

              current_linearization_point.block(introspection.block_indices.velocities)
                = solution.block(introspection.block_indices.velocities);
              current_linearization_point.block(introspection.block_indices.pressure)
                = solution.block(introspection.block_indices.pressure);


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

    // start the timer for periodic checkpoints after the setup above
    time_t last_checkpoint_time = std::time(NULL);


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
        // returned by compute_time_step is unused, will be
        // added to statistics later
        old_time_step = time_step;
        time_step = compute_time_step().first;
        time_step = termination_manager.check_for_last_time_step(time_step);

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
        if (((timestep_number > 0) && (parameters.timing_output_frequency != 0) &&
             (timestep_number % parameters.timing_output_frequency == 0))
            ||
            (parameters.timing_output_frequency == 1))
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

        // check whether to terminate the simulation. the
        // first part of the pair indicates whether to terminate
        // the execution; the second indicates whether to do one
        // more checkpoint
        const std::pair<bool,bool> termination = termination_manager.execute();

        // periodically generate snapshots so that we can resume here
        // if the program aborts or is terminated
        bool do_checkpoint = false;

        // If we base checkpoint frequency on timing, measure the time at process 0
        // This prevents race conditions where some processes will checkpoint and others won't
        if (parameters.checkpoint_time_secs > 0)
          {
            int global_do_checkpoint = ((std::time(NULL)-last_checkpoint_time) >=
                                        parameters.checkpoint_time_secs);
            MPI_Bcast(&global_do_checkpoint, 1, MPI_INT, 0, mpi_communicator);

            do_checkpoint = (global_do_checkpoint == 1);
          }

        // If we base checkpoint frequency on steps, see if it's time for another checkpoint
        if ((parameters.checkpoint_time_secs == 0) &&
            (parameters.checkpoint_steps > 0) &&
            (timestep_number % parameters.checkpoint_steps == 0))
          do_checkpoint = true;

        // Do a checkpoint either if indicated by checkpoint parameters, or if this
        // is the end of simulation and the termination criteria say to checkpoint
        if (do_checkpoint || (termination.first && termination.second))
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

        // see if we want to terminate
        if (termination.first)
          break;
      }
    while (true);
  }
}



// explicit instantiation of the functions we implement in this file
namespace aspect
{
#define INSTANTIATE(dim) \
  template class Simulator<dim>;

  ASPECT_INSTANTIATE(INSTANTIATE)
}
