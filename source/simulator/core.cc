/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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
   * Constructor of the IntermediaryConstructorAction class. Since the
   * class has no members, there is nothing to initialize -- all we
   * need to do is execute the 'action' argument.
   */
  template <int dim>
  Simulator<dim>::IntermediaryConstructorAction::
  IntermediaryConstructorAction (std_cxx1x::function<void ()> action)
  {
    action();
  }



  /**
   * Constructor. Initialize all member variables.
   **/
  template <int dim>
  Simulator<dim>::Simulator (const MPI_Comm mpi_communicator_,
                             ParameterHandler &prm)
    :
    parameters (prm, mpi_communicator_),
    introspection (!parameters.use_direct_stokes_solver,
                   parameters.names_of_compositional_fields),
    mpi_communicator (Utilities::MPI::duplicate_communicator (mpi_communicator_)),
    iostream_tee_device(std::cout, log_file_stream),
    iostream_tee_stream(iostream_tee_device),
    pcout (iostream_tee_stream,
           (Utilities::MPI::
            this_mpi_process(mpi_communicator)
            == 0)),

    computing_timer (pcout, TimerOutput::never,
                     TimerOutput::wall_times),

    geometry_model (GeometryModel::create_geometry_model<dim>(prm)),
    // make sure the parameters object gets a chance to
    // parse those parameters that depend on symbolic names
    // for boundary components
    post_geometry_model_creation_action (std_cxx1x::bind (&Parameters::parse_geometry_dependent_parameters,
                                                          std_cxx1x::ref(parameters),
                                                          std_cxx1x::ref(prm),
                                                          std_cxx1x::cref(*geometry_model))),
    material_model (MaterialModel::create_material_model<dim>(prm)),
    heating_model (HeatingModel::create_heating_model<dim>(prm)),
    gravity_model (GravityModel::create_gravity_model<dim>(prm)),
    // create a boundary temperature model, but only if we actually need
    // it. otherwise, allow the user to simply specify nothing at all
    boundary_temperature (parameters.fixed_temperature_boundary_indicators.empty()
                          ?
                          0
                          :
                          BoundaryTemperature::create_boundary_temperature<dim>(prm)),
    // create a boundary composition model, but only if we actually need
    // it. otherwise, allow the user to simply specify nothing at all
    boundary_composition (parameters.fixed_composition_boundary_indicators.empty()
                          ?
                          0
                          :
                          BoundaryComposition::create_boundary_composition<dim>(prm)),
    initial_conditions (InitialConditions::create_initial_conditions<dim>(prm)),
    compositional_initial_conditions (CompositionalInitialConditions::create_initial_conditions<dim>(prm)),
    adiabatic_conditions (AdiabaticConditions::create_adiabatic_conditions<dim>(prm)),

    time (std::numeric_limits<double>::quiet_NaN()),
    time_step (0),
    old_time_step (0),
    timestep_number (0),

    triangulation (mpi_communicator,
                   typename Triangulation<dim>::MeshSmoothing
                   (Triangulation<dim>::smoothing_on_refinement |
                    Triangulation<dim>::smoothing_on_coarsening),
                   parallel::distributed::Triangulation<dim>::mesh_reconstruction_after_repartitioning),

    //Fourth order mapping doesn't really make sense for free surface
    //calculations (we disable curved boundaries) or we we only have straight
    //boundaries. So we either pick MappingQ(4,true) or MappingQ(1,false)
    mapping (
      (parameters.free_surface_enabled
       ||
       geometry_model->has_curved_elements() == false
      )?1:4,
      (parameters.free_surface_enabled
       ||
       geometry_model->has_curved_elements() == false
      )?false:true),

    // define the finite element. obviously, what we do here needs
    // to match the data we provide in the Introspection class
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
    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
      {
        // only open the logfile on processor 0, the other processors won't be
        // writing into the stream anyway
        log_file_stream.open((parameters.output_directory + "log.txt").c_str(),
                             parameters.resume_computation ? std::ios_base::app : std::ios_base::out);

        // we already printed the header to the screen, so here we just dump it
        // into the logfile.
        print_aspect_header(log_file_stream);
      }

    computing_timer.enter_section("Initialization");

    // first do some error checking for the parameters we got
    {
      // make sure velocity boundary indicators don't appear in multiple lists
      std::set<types::boundary_id> boundary_indicator_lists[4]
        = { parameters.zero_velocity_boundary_indicators,
            parameters.tangential_velocity_boundary_indicators,
            parameters.free_surface_boundary_indicators,
            std::set<types::boundary_id>()
          };

      for (std::map<types::boundary_id,std::pair<std::string, std::string> >::const_iterator
           p = parameters.prescribed_velocity_boundary_indicators.begin();
           p != parameters.prescribed_velocity_boundary_indicators.end();
           ++p)
        boundary_indicator_lists[3].insert (p->first);

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
                         ExcMessage ("Boundary indicator <"
                                     +
                                     Utilities::int_to_string(*intersection.begin())
                                     +
                                     "> with symbolic name <"
                                     +
                                     geometry_model->translate_id_to_symbol_name (*intersection.begin())
                                     +
                                     "> is listed as having more "
                                     "than one type of velocity boundary condition in the input file."));
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

    // if any plugin wants access to the Simulator by deriving from SimulatorAccess, initialize it and
    // call the initialize() functions immediately after.
    //
    // we also need to let all models parse their parameters. this is done *after* setting
    // up their SimulatorAccess base class so that they can query, for example, the
    // geometry model's description of symbolic names for boundary parts. note that
    // the geometry model is the only model whose runtime parameters are already read
    // at the time it is created
    if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(geometry_model.get()))
      sim->initialize (*this);
    geometry_model->initialize ();

    if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(material_model.get()))
      sim->initialize (*this);
    material_model->parse_parameters (prm);
    material_model->initialize ();

    if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(heating_model.get()))
      sim->initialize (*this);
    heating_model->parse_parameters (prm);
    heating_model->initialize ();

    if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(gravity_model.get()))
      sim->initialize (*this);
    gravity_model->parse_parameters (prm);
    gravity_model->initialize ();

    if (boundary_temperature.get())
      {
        if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(boundary_temperature.get()))
          sim->initialize (*this);
        boundary_temperature->parse_parameters (prm);
        boundary_temperature->initialize ();
      }

    if (boundary_composition.get())
      {
        if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(boundary_composition.get()))
          sim->initialize (*this);
        boundary_composition->parse_parameters (prm);
        boundary_composition->initialize ();
      }

    if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(initial_conditions.get()))
      sim->initialize (*this);
    initial_conditions->parse_parameters (prm);
    initial_conditions->initialize ();

    if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(compositional_initial_conditions.get()))
      sim->initialize (*this);
    if (compositional_initial_conditions.get())
      {
        compositional_initial_conditions->parse_parameters (prm);
        compositional_initial_conditions->initialize ();
      }

    if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(adiabatic_conditions.get()))
      sim->initialize (*this);
    adiabatic_conditions->parse_parameters (prm);
    adiabatic_conditions->initialize ();


    // Initialize the free surface handler
    if (parameters.free_surface_enabled)
      {
        //It should be possible to make the free surface work with any of a number of nonlinear
        //schemes, but I do not see a way to do it in generality --IR
        AssertThrow( parameters.nonlinear_solver == NonlinearSolver::IMPES,
                     ExcMessage("The free surface scheme is only implemented for the IMPES solver") );
        //Pressure normalization doesn't really make sense with a free surface, and if we do
        //use it, we can run into problems with geometry_model->depth().
        AssertThrow ( parameters.pressure_normalization == "no",
                      ExcMessage("The free surface scheme can only be used with no pressure normalization") );
        free_surface.reset( new FreeSurfaceHandler( *this, prm ) );
      }

    postprocess_manager.initialize (*this);
    postprocess_manager.parse_parameters (prm);

    mesh_refinement_manager.initialize (*this);
    mesh_refinement_manager.parse_parameters (prm);

    termination_manager.initialize (*this);
    termination_manager.parse_parameters (prm);

    geometry_model->create_coarse_mesh (triangulation);
    global_Omega_diameter = GridTools::diameter (triangulation);

    for (std::map<types::boundary_id,std::pair<std::string,std::string> >::const_iterator
         p = parameters.prescribed_velocity_boundary_indicators.begin();
         p != parameters.prescribed_velocity_boundary_indicators.end();
         ++p)
      {
        VelocityBoundaryConditions::Interface<dim> *bv
          = VelocityBoundaryConditions::create_velocity_boundary_conditions<dim>
            (p->second.second);
        velocity_boundary_conditions[p->first].reset (bv);
        if (dynamic_cast<SimulatorAccess<dim>*>(bv) != 0)
          dynamic_cast<SimulatorAccess<dim>*>(bv)->initialize(*this);
        bv->parse_parameters (prm);
        bv->initialize ();
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

    // we need to do the rhs compatibility modification, if the model is
    // compressible, and there is no open boundary to balance the pressure
    do_pressure_rhs_compatibility_modification = (material_model->is_compressible()
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
         * @param n_components total number of components of the finite element system.
         * @param function_object The scalar function that will form one component
         *     of the resulting Function object.
         */
        VectorFunctionFromVelocityFunctionObject (const unsigned int n_components,
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
    (const unsigned int n_components,
     const std_cxx1x::function<Tensor<1,dim> (const Point<dim> &)> &function_object)
      :
      Function<dim>(n_components),
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

    unsigned int n_stokes_dofs = introspection.system_dofs_per_block[0];
    if (introspection.block_indices.velocities != introspection.block_indices.pressure)
      n_stokes_dofs += introspection.system_dofs_per_block[introspection.block_indices.pressure];

    statistics.add_value("Number of Stokes degrees of freedom", n_stokes_dofs);
    statistics.add_value("Number of temperature degrees of freedom",
                         introspection.system_dofs_per_block[introspection.block_indices.temperature]);
    if (parameters.n_compositional_fields > 0)
      statistics.add_value("Number of degrees of freedom for all compositions",
                           parameters.n_compositional_fields
                           * introspection.system_dofs_per_block[introspection.block_indices.compositional_fields[0]]);

    // then interpolate the current boundary velocities. copy constraints
    // into current_constraints and then add to current_constraints
    compute_current_constraints ();


    //TODO: do this in a more efficient way (TH)? we really only need
    // to make sure that the time dependent velocity boundary conditions
    // end up in the right hand side in the right way; we currently do
    // that by re-assembling the entire system
    if (!velocity_boundary_conditions.empty())
      rebuild_stokes_matrix = rebuild_stokes_preconditioner = true;



    // notify different system components that we started the next time step
    // TODO: implement this for all plugins that might need it at one place.
    // Temperature BC are currently updated in compute_current_constraints
    material_model->update();
    gravity_model->update();
    heating_model->update();
    adiabatic_conditions->update();
  }


  template <int dim>
  void
  Simulator<dim>::
  compute_current_constraints ()
  {
    current_constraints.clear ();
    current_constraints.reinit (introspection.index_sets.system_relevant_set);
    current_constraints.merge (constraints);
    {
      // set the current time and do the interpolation
      // for the prescribed velocity fields
      for (typename std::map<types::boundary_id,std_cxx1x::shared_ptr<VelocityBoundaryConditions::Interface<dim> > >::iterator
           p = velocity_boundary_conditions.begin();
           p != velocity_boundary_conditions.end(); ++p)
        {
          p->second->update ();
          VectorFunctionFromVelocityFunctionObject<dim> vel
          (introspection.n_components,
           std_cxx1x::bind (&VelocityBoundaryConditions::Interface<dim>::boundary_velocity,
                            p->second,
                            std_cxx1x::_1));

          // here we create a mask for interpolate_boundary_values out of the 'selector'
          std::vector<bool> mask(introspection.component_masks.velocities.size(), false);
          Assert(introspection.component_masks.velocities[0]==true,
                 ExcInternalError()); // in case we ever move the velocity around
          const std::string &comp = parameters.prescribed_velocity_boundary_indicators[p->first].first;

          if (comp.length()>0)
            {
              for (std::string::const_iterator direction=comp.begin(); direction!=comp.end(); ++direction)
                {
                  switch (*direction)
                    {
                      case 'x':
                        mask[0] = true;
                        break;
                      case 'y':
                        mask[1] = true;
                        break;
                      case 'z':
                        // we must be in 3d, or 'z' should never have gotten through
                        Assert (dim==3, ExcInternalError());
                        mask[2] = true;
                        break;
                      default:
                        Assert (false, ExcInternalError());
                    }
                }
            }
          else
            {
              // no mask given -- take all velocities
              for (unsigned int i=0; i<introspection.component_masks.velocities.size(); ++i)
                mask[i]=introspection.component_masks.velocities[i];

              Assert(introspection.component_masks.velocities[0]==true,
                     ExcInternalError()); // in case we ever move the velocity down
            }

          VectorTools::interpolate_boundary_values (dof_handler,
                                                    p->first,
                                                    vel,
                                                    current_constraints,
                                                    mask);
        }
    }

    // do the same for the temperature variable: evaluate the current boundary temperature
    // and add these constraints as well
    {
      // If there is a fixed boundary temperature,
      // update the temperature boundary condition.
      if (boundary_temperature.get())
        boundary_temperature->update();

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
                                                    current_constraints,
                                                    introspection.component_masks.temperature);
        }
    }

    // now do the same for the composition variable:
    {
      // If there are fixed boundary compositions,
      // update the composition boundary condition.
      if (boundary_composition.get())
        boundary_composition->update();

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
                                                      current_constraints,
                                                      introspection.component_masks.compositional_fields[c]);
          }
    }
    current_constraints.close();
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

    LinearAlgebra::BlockCompressedSparsityPattern sp;
#ifdef ASPECT_USE_PETSC
    sp.reinit (introspection.index_sets.system_relevant_partitioning);
#else
    sp.reinit (system_partitioning,
               system_partitioning,
               introspection.index_sets.system_relevant_partitioning,
               mpi_communicator);
#endif

    DoFTools::make_sparsity_pattern (dof_handler,
                                     coupling, sp,
                                     constraints, false,
                                     Utilities::MPI::
                                     this_mpi_process(mpi_communicator));

#ifdef ASPECT_USE_PETSC
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

    LinearAlgebra::BlockCompressedSparsityPattern sp;
#ifdef ASPECT_USE_PETSC
    sp.reinit (introspection.index_sets.system_relevant_partitioning);
#else
    sp.reinit (system_partitioning,
               system_partitioning,
               introspection.index_sets.system_relevant_partitioning,
               mpi_communicator);
#endif

    DoFTools::make_sparsity_pattern (dof_handler,
                                     coupling, sp,
                                     constraints, false,
                                     Utilities::MPI::
                                     this_mpi_process(mpi_communicator));
#ifdef ASPECT_USE_PETSC
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

      for (periodic_boundary_set::iterator p = pbs.begin(); p != pbs.end(); ++p)
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

          DoFTools::make_periodicity_constraints(dof_handler,
                                                 (*p).first.first,  //first boundary id
                                                 (*p).first.second, //second boundary id
                                                 (*p).second,       //cartesian direction for translational symmetry
                                                 constraints);
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

    if (parameters.free_surface_enabled)
      free_surface->setup_dofs();
    setup_nullspace_removal();

    computing_timer.exit_section();
  }


  /**
   * This is an internal deal.II function stolen from dof_tools.cc
   */
  template <int dim, int spacedim>
  std::vector<unsigned char>
  get_local_component_association (const FiniteElement<dim,spacedim>  &fe,
                                   const ComponentMask        &component_mask)
  {
    std::vector<unsigned char> local_component_association (fe.dofs_per_cell,
                                                            (unsigned char)(-1));

    // compute the component each local dof belongs to.
    // if the shape function is primitive, then this
    // is simple and we can just associate it with
    // what system_to_component_index gives us
    for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
      if (fe.is_primitive(i))
        local_component_association[i] =
          fe.system_to_component_index(i).first;
      else
        // if the shape function is not primitive, then either use the
        // component of the first nonzero component corresponding
        // to this shape function (if the component is not specified
        // in the component_mask), or the first component of this block
        // that is listed in the component_mask (if the block this
        // component corresponds to is indeed specified in the component
        // mask)
        {
          const unsigned int first_comp =
            fe.get_nonzero_components(i).first_selected_component();

          if ((fe.get_nonzero_components(i)
               &
               component_mask).n_selected_components(fe.n_components()) == 0)
            local_component_association[i] = first_comp;
          else
            // pick the component selected. we know from the previous 'if'
            // that within the components that are nonzero for this
            // shape function there must be at least one for which the
            // mask is true, so we will for sure run into the break()
            // at one point
            for (unsigned int c=first_comp; c<fe.n_components(); ++c)
              if (component_mask[c] == true)
                {
                  local_component_association[i] = c;
                  break;
                }
        }

    Assert (std::find (local_component_association.begin(),
                       local_component_association.end(),
                       (unsigned char)(-1))
            ==
            local_component_association.end(),
            ExcInternalError());

    return local_component_association;
  }

  /**
   * Returns an IndexSet that contains all locally active DoFs that belong to
   * the given component_mask.
   *
   * This function should be moved into deal.II at some point.
   */
  template <int dim>
  IndexSet extract_component_subset(DoFHandler<dim> &dof_handler, const ComponentMask &component_mask)
  {
    std::vector<unsigned char> local_asoc =
      get_local_component_association (dof_handler.get_fe(),
                                       ComponentMask(dof_handler.get_fe().n_components(), true));

    IndexSet ret(dof_handler.n_dofs());

    unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;
    std::vector<types::global_dof_index> indices(dofs_per_cell);
    for (typename DoFHandler<dim>::active_cell_iterator cell=dof_handler.begin_active();
         cell!=dof_handler.end(); ++cell)
      if (cell->is_locally_owned())
        {
          cell->get_dof_indices(indices);
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            if (component_mask[local_asoc[i]])
              ret.add_index(indices[i]);
        }

    return ret;
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
      types::global_dof_index n_u = introspection.system_dofs_per_block[0],
                              n_p = introspection.system_dofs_per_block[introspection.block_indices.pressure],
                              n_T = introspection.system_dofs_per_block[introspection.block_indices.temperature];

      Assert(!parameters.use_direct_stokes_solver ||
             (introspection.block_indices.velocities == introspection.block_indices.pressure),
             ExcInternalError());

      // only count pressure once if velocity and pressure are in the same block,
      // i.e., direct solver is used.
      if (introspection.block_indices.velocities == introspection.block_indices.pressure)
        n_p = 0;

      std::vector<types::global_dof_index> n_C (parameters.n_compositional_fields);
      for (unsigned int c=0; c<parameters.n_compositional_fields; ++c)
        n_C[c] = introspection.system_dofs_per_block
                 [introspection.block_indices.compositional_fields[c]];


      IndexSet system_index_set = dof_handler.locally_owned_dofs();
      introspection.index_sets.system_partitioning.clear ();
      introspection.index_sets.system_partitioning.push_back(system_index_set.get_view(0,n_u));
      if (n_p != 0)
        {
          introspection.index_sets.locally_owned_pressure_dofs = system_index_set.get_view(n_u,n_u+n_p);
          introspection.index_sets.system_partitioning.push_back(introspection.index_sets.locally_owned_pressure_dofs);
        }
      else
        {
          introspection.index_sets.locally_owned_pressure_dofs = system_index_set & extract_component_subset(dof_handler, introspection.component_masks.pressure);
        }
      introspection.index_sets.system_partitioning.push_back(system_index_set.get_view(n_u+n_p,n_u+n_p+n_T));

      introspection.index_sets.stokes_partitioning.clear ();
      introspection.index_sets.stokes_partitioning.push_back(system_index_set.get_view(0,n_u));
      if (n_p != 0)
        {
          introspection.index_sets.stokes_partitioning.push_back(system_index_set.get_view(n_u,n_u+n_p));
        }

      DoFTools::extract_locally_relevant_dofs (dof_handler,
                                               introspection.index_sets.system_relevant_set);
      introspection.index_sets.system_relevant_partitioning.clear ();
      introspection.index_sets.system_relevant_partitioning
      .push_back(introspection.index_sets.system_relevant_set.get_view(0,n_u));
      if (n_p != 0)
        {
          introspection.index_sets.system_relevant_partitioning
          .push_back(introspection.index_sets.system_relevant_set.get_view(n_u,n_u+n_p));
        }
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

    mesh_refinement_manager.tag_additional_cells ();


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

    if (parameters.free_surface_enabled)
      x_system.push_back( &free_surface->mesh_velocity );

    parallel::distributed::SolutionTransfer<dim,LinearAlgebra::BlockVector>
    system_trans(dof_handler);

    std::vector<const LinearAlgebra::Vector *> x_fs_system (1);
    std::auto_ptr<parallel::distributed::SolutionTransfer<dim,LinearAlgebra::Vector> >
    freesurface_trans;

    if (parameters.free_surface_enabled)
      {
        x_fs_system[0] = &free_surface->mesh_vertices;
        freesurface_trans.reset (new parallel::distributed::SolutionTransfer<dim,LinearAlgebra::Vector>
                                 (free_surface->free_surface_dof_handler));
      }

    triangulation.prepare_coarsening_and_refinement();
    system_trans.prepare_for_coarsening_and_refinement(x_system);

    if (parameters.free_surface_enabled)
      freesurface_trans->prepare_for_coarsening_and_refinement(x_fs_system);

    triangulation.execute_coarsening_and_refinement ();
    global_volume = GridTools::volume (triangulation, mapping);
    computing_timer.exit_section();

    setup_dofs ();

    computing_timer.enter_section ("Refine mesh structure, part 2");

    {
      LinearAlgebra::BlockVector distributed_system;
      LinearAlgebra::BlockVector old_distributed_system;
      LinearAlgebra::BlockVector distributed_mesh_velocity;

      distributed_system.reinit(introspection.index_sets.system_partitioning, mpi_communicator);
      old_distributed_system.reinit(introspection.index_sets.system_partitioning, mpi_communicator);
      if (parameters.free_surface_enabled)
        distributed_mesh_velocity.reinit(introspection.index_sets.system_partitioning, mpi_communicator);

      std::vector<LinearAlgebra::BlockVector *> system_tmp (2);
      system_tmp[0] = &distributed_system;
      system_tmp[1] = &old_distributed_system;

      if (parameters.free_surface_enabled)
        system_tmp.push_back(&distributed_mesh_velocity);

      // transfer the data previously stored into the vectors indexed by
      // system_tmp. then ensure that the interpolated solution satisfies
      // hanging node constraints
      //
      // note that the 'constraints' variable contains hanging node constraints
      // and constraints from periodic boundary conditions (as well as from
      // zero and tangential velocity boundary conditions), but not from
      // non-homogeneous boundary conditions. the latter are added not in setup_dofs(),
      // which we call above, but added to 'current_constraints' in start_timestep(),
      // which we do not want to call here.
      //
      // however, what we have should be sufficient: we have everything that
      // is necessary to make the solution vectors *conforming* on the current
      // mesh.
      system_trans.interpolate (system_tmp);

      constraints.distribute (distributed_system);
      solution     = distributed_system;

      constraints.distribute (old_distributed_system);
      old_solution = old_distributed_system;

      if (parameters.free_surface_enabled)
        {
          constraints.distribute (distributed_mesh_velocity);
          free_surface->mesh_velocity = distributed_mesh_velocity;
        }
    }

    // do the same as above also for the free surface solution
    if (parameters.free_surface_enabled)
      {
        LinearAlgebra::Vector distributed_mesh_vertices;
        distributed_mesh_vertices.reinit(free_surface->mesh_locally_owned,
                                         mpi_communicator);

        std::vector<LinearAlgebra::Vector *> system_tmp (1);
        system_tmp[0] = &distributed_mesh_vertices;

        freesurface_trans->interpolate (system_tmp);
        free_surface->mesh_vertex_constraints.distribute (distributed_mesh_vertices);
        free_surface->mesh_vertices = distributed_mesh_vertices;
      }

    if (parameters.free_surface_enabled)
      free_surface->displace_mesh ();

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
          //We do the free surface execution at the beginning of the timestep for a specific reason.
          //The time step size is calculated AFTER the whole solve_timestep() function.  If we call
          //free_surface_execute() after the Stokes solve, it will be before we know what the appropriate
          //time step to take is, and we will timestep the boundary incorrectly.
          if (parameters.free_surface_enabled)
            free_surface->execute ();

          assemble_advection_system (AdvectionField::temperature());
          build_advection_preconditioner(AdvectionField::temperature(),
                                         T_preconditioner);
          solve_advection(AdvectionField::temperature());

          current_linearization_point.block(introspection.block_indices.temperature)
            = solution.block(introspection.block_indices.temperature);

          for (unsigned int c=0; c<parameters.n_compositional_fields; ++c)
            {
              assemble_advection_system (AdvectionField::composition(c));
              build_advection_preconditioner(AdvectionField::composition(c),
                                             C_preconditioner);
              solve_advection(AdvectionField::composition(c));
            }

          for (unsigned int c=0; c<parameters.n_compositional_fields; ++c)
            current_linearization_point.block(introspection.block_indices.compositional_fields[c])
              = solution.block(introspection.block_indices.compositional_fields[c]);

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
              assemble_advection_system(AdvectionField::temperature());

              if (iteration == 0)
                {
                  build_advection_preconditioner(AdvectionField::temperature(),
                                                 T_preconditioner);
                  initial_temperature_residual = system_rhs.block(introspection.block_indices.temperature).l2_norm();
                }

              const double temperature_residual = solve_advection(AdvectionField::temperature());

              current_linearization_point.block(introspection.block_indices.temperature)
                = solution.block(introspection.block_indices.temperature);
              rebuild_stokes_matrix = true;
              std::vector<double> composition_residual (parameters.n_compositional_fields,0);

              for (unsigned int c=0; c<parameters.n_compositional_fields; ++c)
                {
                  assemble_advection_system (AdvectionField::composition(c));
                  build_advection_preconditioner(AdvectionField::composition(c),
                                                 C_preconditioner);
                  if (iteration == 0)
                    initial_composition_residual[c] = system_rhs.block(introspection.block_indices.compositional_fields[c]).l2_norm();

                  composition_residual[c]
                    = solve_advection(AdvectionField::composition(c));
                }

              // for consistency we update the current linearization point only after we have solved
              // all fields, so that we use the same point in time for every field when solving
              for (unsigned int c=0; c<parameters.n_compositional_fields; ++c)
                current_linearization_point.block(introspection.block_indices.compositional_fields[c])
                  = solution.block(introspection.block_indices.compositional_fields[c]);

              // the Stokes matrix depends on the viscosity. if the viscosity
              // depends on other solution variables, then after we need to
              // update the Stokes matrix in every time step and so need to set
              // the following flag. if we change the Stokes matrix we also
              // need to update the Stokes preconditioner.
              if (stokes_matrix_depends_on_solution() == true)
                rebuild_stokes_matrix = rebuild_stokes_preconditioner = true;

              assemble_stokes_system();
              if (iteration == 0)
                {
                  build_stokes_preconditioner();
                  initial_stokes_residual = compute_initial_stokes_residual();
                }

              const double stokes_residual = solve_stokes();

              current_linearization_point = solution;

              // write the residual output in the same order as the output when
              // solving the equations
              pcout << "      Nonlinear residuals: " << temperature_residual;

              for (unsigned int c=0; c<parameters.n_compositional_fields; ++c)
                pcout << ", " << composition_residual[c];

              pcout << ", " << stokes_residual;

              pcout << std::endl
                    << std::endl;

              double max = 0.0;
              for (unsigned int c=0; c<parameters.n_compositional_fields; ++c)
                if (initial_composition_residual[c]>0)
                  max = std::max(composition_residual[c]/initial_composition_residual[c],max);
              if (initial_stokes_residual>0)
                max = std::max(stokes_residual/initial_stokes_residual, max);
              if (initial_temperature_residual>0)
                max = std::max(temperature_residual/initial_temperature_residual, max);
              pcout << "      residual: " << max << std::endl;
              if (max < parameters.nonlinear_tolerance)
                break;

              ++iteration;
//TODO: terminate here if the number of iterations is too large and we see no convergence
            }
          while (iteration < parameters.max_nonlinear_iterations);

          break;
        }

        case NonlinearSolver::iterated_Stokes:
        {
          // solve the temperature system once...
          assemble_advection_system (AdvectionField::temperature());
          build_advection_preconditioner (AdvectionField::temperature (),
                                          T_preconditioner);
          solve_advection(AdvectionField::temperature());
          current_linearization_point.block(introspection.block_indices.temperature)
            = solution.block(introspection.block_indices.temperature);

          for (unsigned int c=0; c<parameters.n_compositional_fields; ++c)
            {
              assemble_advection_system (AdvectionField::composition(c));
              build_advection_preconditioner (AdvectionField::composition (c),
                                              C_preconditioner);
              solve_advection(AdvectionField::composition(c));
            }

          for (unsigned int c=0; c<parameters.n_compositional_fields; ++c)
            current_linearization_point.block(introspection.block_indices.compositional_fields[c])
              = solution.block(introspection.block_indices.compositional_fields[c]);

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
                  pcout << "   Residual after nonlinear iteration " << i+1 << ": " << stokes_residual/initial_stokes_residual << std::endl;
                  if (stokes_residual/initial_stokes_residual < parameters.nonlinear_tolerance)
                    {
                      break; // convergence reached, exist nonlinear iteration.
                    }
                }

              current_linearization_point.block(introspection.block_indices.velocities)
                = solution.block(introspection.block_indices.velocities);
              if (introspection.block_indices.velocities != introspection.block_indices.pressure)
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
        // Instead of calling global_refine(n) we flag all cells for
        // refinement and then allow the mesh refinement plugins to unflag
        // the cells if desired. This procedure is repeated n times. If there
        // is no plugin that modifies the flags, it is equivalent to
        // refine_global(n).
        for (unsigned int n=0; n<parameters.initial_global_refinement; ++n)
          {
            for (typename Triangulation<dim>::active_cell_iterator
                 cell = triangulation.begin_active();
                 cell != triangulation.end(); ++cell)
              cell->set_refine_flag ();

            mesh_refinement_manager.tag_additional_cells ();
            triangulation.execute_coarsening_and_refinement();
          }

        global_volume = GridTools::volume (triangulation, mapping);

        time                      = parameters.start_time;
        timestep_number           = 0;
        time_step = old_time_step = 0;

        setup_dofs();
      }

    // start the timer for periodic checkpoints after the setup above
    time_t last_checkpoint_time = std::time(NULL);


  start_time_iteration:

    if (parameters.resume_computation == false)
      {
        computing_timer.enter_section ("Initialization");

        time                      = parameters.start_time;
        timestep_number           = 0;
        time_step = old_time_step = 0;

        set_initial_temperature_and_compositional_fields ();
        compute_initial_pressure_field ();

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
        time_step = std::min (compute_time_step().first,
                              parameters.maximum_time_step);
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
            if (parameters.timing_output_frequency ==0)
              computing_timer.print_summary ();

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
          if (
            (timestep_number > 0
             &&
             (parameters.adaptive_refinement_interval > 0)
             &&
             (timestep_number % parameters.adaptive_refinement_interval == 0))
            ||
            (timestep_number==0 && parameters.adaptive_refinement_interval == 1)
          )
            refine_mesh (max_refinement_level);

        // every n time steps output a summary of the current
        // timing information
        if (((timestep_number > 0) && (parameters.timing_output_frequency != 0) &&
             (timestep_number % parameters.timing_output_frequency == 0))
            ||
            (parameters.timing_output_frequency == 1)||(parameters.timing_output_frequency == 0))
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

    // we disable automatic summary printing so that it won't happen when
    // throwing an exception. Therefore, we have to do this manually here:
    computing_timer.print_summary ();
  }
}



// explicit instantiation of the functions we implement in this file
namespace aspect
{
#define INSTANTIATE(dim) \
  template class Simulator<dim>;

  ASPECT_INSTANTIATE(INSTANTIATE)
}
