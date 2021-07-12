/*
  Copyright (C) 2011 - 2020 by the authors of the ASPECT code.

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


#include <aspect/simulator.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/initial_temperature/interface.h>
#include <aspect/initial_composition/interface.h>
#include <aspect/postprocess/particles.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>


namespace aspect
{

  template <int dim>
  void Simulator<dim>::set_initial_temperature_and_compositional_fields ()
  {
    // create a fully distributed vector since we
    // need to write into it and we can not
    // write into vectors with ghost elements
    LinearAlgebra::BlockVector initial_solution;
    double max_sum_comp = 0.0;

    // we need to track whether we need to normalize the totality of fields
    bool normalize_composition = false;

    initial_solution.reinit(system_rhs, false);

    // below, we would want to call VectorTools::interpolate on the
    // entire FESystem. there currently is no way to restrict the
    // interpolation operations to only a subset of vector
    // components (oversight in deal.II?), specifically to the
    // temperature component. this causes more work than necessary
    // but worse yet, it doesn't work for the DGP(q) pressure element
    // if we use a locally conservative formulation since there the
    // pressure element is non-interpolating (we get an exception
    // even though we are, strictly speaking, not interested in
    // interpolating the pressure; but, as mentioned, there is no way
    // to tell VectorTools::interpolate that)
    //
    // to work around this problem, the following code is essentially
    // a (simplified) copy of the code in VectorTools::interpolate
    // that only works on the temperature component
    //
    // TODO: it would be great if we had a cleaner way than iterating to 1+n_fields.
    // Additionally, the n==1 logic for normalization at the bottom is not pretty.
    for (unsigned int n=0; n<1+introspection.n_compositional_fields; ++n)
      {
        AdvectionField advf = ((n == 0) ? AdvectionField::temperature()
                               : AdvectionField::composition(n-1));

        const unsigned int base_element = advf.base_element(introspection);

        // get the temperature/composition support points
        const std::vector<Point<dim>> support_points
                                   = finite_element.base_element(base_element).get_unit_support_points();
        Assert (support_points.size() != 0,
                ExcInternalError());

        // create an FEValues object with just the temperature/composition element
        FEValues<dim> fe_values (*mapping, finite_element,
                                 support_points,
                                 update_quadrature_points);

        std::vector<types::global_dof_index> local_dof_indices (finite_element.dofs_per_cell);

        const VectorFunctionFromScalarFunctionObject<dim, double> &advf_init_function
          =
            (advf.is_temperature()
             ?
             VectorFunctionFromScalarFunctionObject<dim, double>(
               [&](const Point<dim> &p) -> double
        {
          return initial_temperature_manager.initial_temperature(p);
        },
        introspection.component_indices.temperature,
        introspection.n_components)
        :
        VectorFunctionFromScalarFunctionObject<dim, double>(
          [&](const Point<dim> &p) -> double
        {
          return initial_composition_manager.initial_composition(p, n-1);
        },
        introspection.component_indices.compositional_fields[n-1],
        introspection.n_components));

        const ComponentMask advf_mask =
          (advf.is_temperature()
           ?
           introspection.component_masks.temperature
           :
           introspection.component_masks.compositional_fields[n-1]);

        VectorTools::interpolate(*mapping,
                                 dof_handler,
                                 advf_init_function,
                                 initial_solution,
                                 advf_mask);

        if (parameters.normalized_fields.size()>0 && n==1)
          for (const auto &cell : dof_handler.active_cell_iterators())
            if (cell->is_locally_owned())
              {
                fe_values.reinit (cell);

                // Go through the support points for each dof
                for (unsigned int i=0; i<finite_element.base_element(base_element).dofs_per_cell; ++i)
                  {
                    // if it is specified in the parameter file that the sum of all compositional fields
                    // must not exceed one, this should be checked
                    double sum = 0;
                    for (unsigned int m=0; m<parameters.normalized_fields.size(); ++m)
                      sum += initial_composition_manager.initial_composition(fe_values.quadrature_point(i),
                                                                             parameters.normalized_fields[m]);

                    if (std::abs(sum) > 1.0+std::numeric_limits<double>::epsilon())
                      {
                        max_sum_comp = std::max(sum, max_sum_comp);
                        normalize_composition = true;
                      }
                  }
              }

        initial_solution.compress(VectorOperation::insert);

        // if at least one processor decides that it needs
        // to normalize, do the same on all processors.
        if (Utilities::MPI::max (normalize_composition ? 1 : 0,
                                 mpi_communicator)
            == 1)
          {
            const double global_max
              = Utilities::MPI::max (max_sum_comp, mpi_communicator);

            if (n==1)
              pcout << "Sum of compositional fields is not one, fields will be normalized"
                    << std::endl;

            for (unsigned int m=0; m<parameters.normalized_fields.size(); ++m)
              if (n-1==parameters.normalized_fields[m])
                initial_solution.block(introspection.block_indices.compositional_fields[n-1]) /= global_max;
          }
      }

    // then apply constraints and copy the
    // result into vectors with ghost elements. to do so,
    // we need the current constraints to be correct for
    // the current time
    compute_current_constraints ();
    current_constraints.distribute(initial_solution);

    // Now copy the temperature and initial composition blocks into the solution variables

    for (unsigned int n=0; n<1+introspection.n_compositional_fields; ++n)
      {
        AdvectionField advf = ((n == 0) ? AdvectionField::temperature()
                               : AdvectionField::composition(n-1));

        const unsigned int blockidx = advf.block_index(introspection);

        solution.block(blockidx) = initial_solution.block(blockidx);
        old_solution.block(blockidx) = initial_solution.block(blockidx);
        old_old_solution.block(blockidx) = initial_solution.block(blockidx);
        current_linearization_point.block(blockidx) = initial_solution.block(blockidx);
      }
  }


  template <int dim>
  void Simulator<dim>::interpolate_particle_properties (const AdvectionField &advection_field)
  {
    TimerOutput::Scope timer (computing_timer, "Particles: Interpolate");

    // below, we would want to call VectorTools::interpolate on the
    // entire FESystem. there currently is no way to restrict the
    // interpolation operations to only a subset of vector
    // components (oversight in deal.II?), specifically to the
    // temperature component. this causes more work than necessary
    // but worse yet, it doesn't work for the DGP(q) pressure element
    // if we use a locally conservative formulation since there the
    // pressure element is non-interpolating (we get an exception
    // even though we are, strictly speaking, not interested in
    // interpolating the pressure; but, as mentioned, there is no way
    // to tell VectorTools::interpolate that)
    //
    // to work around this problem, the following code is essentially
    // a (simplified) copy of the code in VectorTools::interpolate
    // that only works on the given component

    // create a fully distributed vector since we
    // need to write into it and we can not
    // write into vectors with ghost elements

    const Postprocess::Particles<dim> &particle_postprocessor =
      postprocess_manager.template get_matching_postprocessor<Postprocess::Particles<dim>>();

    const Particle::Interpolator::Interface<dim> *particle_interpolator = &particle_postprocessor.get_particle_world().get_interpolator();
    const Particle::Property::Manager<dim> *particle_property_manager = &particle_postprocessor.get_particle_world().get_property_manager();

    unsigned int particle_property;

    if (parameters.mapped_particle_properties.size() != 0)
      {
        const std::pair<std::string,unsigned int> particle_property_and_component = parameters.mapped_particle_properties.find(advection_field.compositional_variable)->second;

        particle_property = particle_property_manager->get_data_info().get_position_by_field_name(particle_property_and_component.first)
                            + particle_property_and_component.second;
      }
    else
      {
        particle_property = std::count(introspection.compositional_field_methods.begin(),
                                       introspection.compositional_field_methods.begin() + advection_field.compositional_variable,
                                       Parameters<dim>::AdvectionFieldMethod::particles);
        AssertThrow(particle_property <= particle_property_manager->get_data_info().n_components(),
                    ExcMessage("Can not automatically match particle properties to fields, because there are"
                               "more fields that are marked as particle advected than particle properties"));
      }

    LinearAlgebra::BlockVector particle_solution;

    particle_solution.reinit(system_rhs, false);

    const unsigned int base_element = advection_field.base_element(introspection);

    // get the temperature/composition support points
    const std::vector<Point<dim>> support_points
                               = finite_element.base_element(base_element).get_unit_support_points();
    Assert (support_points.size() != 0,
            ExcInternalError());

    // create an FEValues object with just the temperature/composition element
    FEValues<dim> fe_values (*mapping, finite_element,
                             support_points,
                             update_quadrature_points);

    std::vector<types::global_dof_index> local_dof_indices (finite_element.dofs_per_cell);

    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          fe_values.reinit (cell);
          const std::vector<Point<dim>> quadrature_points = fe_values.get_quadrature_points();

          ComponentMask property_mask  (particle_property_manager->get_data_info().n_components(),false);
          property_mask.set(particle_property,true);

          std::vector<std::vector<double>> particle_properties;

          try
            {
              particle_properties =
                particle_interpolator->properties_at_points(particle_postprocessor.get_particle_world().get_particle_handler(),
                                                            quadrature_points,
                                                            property_mask,
                                                            cell);
            }
          // interpolators that throw exceptions usually do not result in
          // anything good, because they result in an unwinding of the stack
          // and, if only one processor triggers an exception, the
          // destruction of objects often causes a deadlock or completely
          // unrelated MPI error messages. Thus, if an exception is
          // generated, catch it, print an error message, and abort the program.
          catch (std::exception &exc)
            {
              std::cerr << std::endl << std::endl
                        << "----------------------------------------------------"
                        << std::endl;
              std::cerr << "Exception on MPI process <"
                        << Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)
                        << "> while interpolating particle properties: "
                        << std::endl
                        << exc.what() << std::endl
                        << "Aborting!" << std::endl
                        << "----------------------------------------------------"
                        << std::endl;

              // terminate the program!
              MPI_Abort (MPI_COMM_WORLD, 1);
            }

          // go through the composition dofs and set their global values
          // to the particle field interpolated at these points
          cell->get_dof_indices (local_dof_indices);
          for (unsigned int i=0; i<finite_element.base_element(base_element).dofs_per_cell; ++i)
            {
              const unsigned int system_local_dof
                = finite_element.component_to_system_index(advection_field.component_index(introspection),
                                                           /*dof index within component=*/i);

              particle_solution(local_dof_indices[system_local_dof]) = particle_properties[i][particle_property];
            }
        }

    particle_solution.compress(VectorOperation::insert);

    // we should not have written at all into any of the blocks with
    // the exception of the current composition block
    for (unsigned int b=0; b<particle_solution.n_blocks(); ++b)
      if (b != advection_field.block_index(introspection))
        Assert (particle_solution.block(b).l2_norm() == 0,
                ExcInternalError());

    // overwrite the relevant composition block only
    const unsigned int blockidx = advection_field.block_index(introspection);
    solution.block(blockidx) = particle_solution.block(blockidx);

    // In the first timestep, and for iterative Advection schemes only
    // in the first nonlinear iteration, initialize all solution vectors with the initial
    // particle solution, identical to the end of the
    // Simulator<dim>::set_initial_temperature_and_compositional_fields ()
    // function.
    if (timestep_number == 0 && nonlinear_iteration == 0)
      {
        old_solution.block(blockidx) = particle_solution.block(blockidx);
        old_old_solution.block(blockidx) = particle_solution.block(blockidx);
      }
  }


  template <int dim>
  void Simulator<dim>::compute_initial_pressure_field ()
  {
    // Note that this code will overwrite the velocity solution with 0 if
    // velocity and pressure are in the same block (i.e., direct solver is
    // used). As the velocity is all zero anyway, this is currently not a
    // problem.

    // we'd like to interpolate the initial pressure onto the pressure
    // variable but that's a bit involved because the pressure may either
    // be an FE_Q (for which we can interpolate) or an FE_DGP (for which
    // we can't since the element has no nodal basis.
    //
    // fortunately, in the latter case, the element is discontinuous and
    // we can compute a local projection onto the pressure space
    if (parameters.use_locally_conservative_discretization == false)
      {
        // allocate a vector that is distributed but doesn't have
        // ghost elements (vectors with ghost elements are not
        // writable); the stokes_rhs vector is a valid template for
        // this kind of thing. interpolate into it and later copy it into the
        // solution vector that does have the necessary ghost elements
        LinearAlgebra::BlockVector system_tmp;
        system_tmp.reinit (system_rhs);

        // First grab the correct pressure to work on:
        const FEVariable<dim> &pressure_variable
          = parameters.include_melt_transport ?
            introspection.variable("fluid pressure")
            : introspection.variable("pressure");
        const unsigned int pressure_comp = pressure_variable.first_component_index;
        const ComponentMask pressure_component_mask = pressure_variable.component_mask;

        // interpolate the pressure given by the adiabatic conditions
        // object onto the solution space. note that interpolate
        // wants a function that represents all components of the
        // solution vector, so create such a function object
        // that is simply zero for all velocity components
        auto lambda = [&](const Point<dim> &p) -> double
        {
          return adiabatic_conditions->pressure(p);
        };

        VectorFunctionFromScalarFunctionObject<dim> vector_function_object(
          lambda,
          pressure_comp,
          introspection.n_components);

        VectorTools::interpolate (*mapping, dof_handler,
                                  vector_function_object,
                                  system_tmp,
                                  pressure_component_mask);

        // we may have hanging nodes, so apply constraints
        constraints.distribute (system_tmp);

        const unsigned int pressure_block = pressure_variable.block_index;
        old_solution.block(pressure_block) = system_tmp.block(pressure_block);
      }
    else
      {
        // Find the local projection for the discontinuous pressure
        // element. This is only going to work if, indeed, the element
        // is discontinuous.
        Assert (finite_element.base_element(introspection.base_elements.pressure).dofs_per_face == 0,
                ExcNotImplemented());

        LinearAlgebra::BlockVector system_tmp;
        system_tmp.reinit (system_rhs);

        QGauss<dim> quadrature(parameters.stokes_velocity_degree+1);

        Utilities::project_cellwise<dim,LinearAlgebra::BlockVector>(*mapping,
                                                                    dof_handler,
                                                                    introspection.component_indices.pressure,
                                                                    quadrature,
                                                                    [&](const typename DoFHandler<dim>::active_cell_iterator & /*cell*/,
                                                                        const std::vector<Point<dim>> &q_points,
                                                                        std::vector<double> &values) -> void
        {
          for (unsigned int i=0; i<values.size(); ++i)
            values[i] = adiabatic_conditions->pressure(q_points[i]);
          return;
        },
        system_tmp);

        old_solution.block(introspection.block_indices.pressure) = system_tmp.block(introspection.block_indices.pressure);
      }

    // normalize the pressure in such a way that the surface pressure
    // equals a known and desired value
    this->last_pressure_normalization_adjustment = normalize_pressure(old_solution);

    // set all solution vectors to the same value as the previous solution
    solution = old_solution;
    old_old_solution = old_solution;
  }
}



// explicit instantiation of the functions we implement in this file
namespace aspect
{
#define INSTANTIATE(dim) \
  template void Simulator<dim>::set_initial_temperature_and_compositional_fields(); \
  template void Simulator<dim>::compute_initial_pressure_field(); \
  template void Simulator<dim>::interpolate_particle_properties(const AdvectionField &);


  ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}
