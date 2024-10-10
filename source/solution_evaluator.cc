/*
  Copyright (C) 2024 - 2024 by the authors of the ASPECT code.

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

#include <aspect/global.h>
#include <aspect/utilities.h>
#include <aspect/solution_evaluator.h>
#include <aspect/simulator.h>
#include <aspect/melt.h>

namespace aspect
{
  namespace internal
  {
    /**
     * The functions in this namespace allow us to use scalar and vector-valued FEEvaluation objects in the same
     * way, as the return type of FEEvaluation is double for one component, but a Tensor for more than one
     * component. These function converts scalar to Tensor return values.
     */
    namespace convert
    {
      template <int dim>
      Tensor<1,dim> to_tensor(const Tensor<1,dim> &in)
      {
        return in;
      }



      template <int dim>
      Tensor<1,1> to_tensor(const double in)
      {
        Tensor<1,1> result;
        result[0] = in;
        return result;
      }



      template <int dim, int n_components>
      Tensor<1,n_components,Tensor<1,dim>> to_tensor2(const Tensor<1,n_components,Tensor<1,dim>> &in)
      {
        return in;
      }



      template <int dim, int n_components>
      Tensor<1,1,Tensor<1,dim>> to_tensor2(const Tensor<1,dim> &in)
      {
        Tensor<1,1,Tensor<1,dim>> result;
        result[0] = in;
        return result;
      }
    }



    /**
     * Implementation of the base class DynamicFEPointEvaluation that wraps
     * an FEPointEvaluation object.
     */
    template <int dim, int n_components>
    class DynamicFEPointEvaluationImpl: public DynamicFEPointEvaluation<dim>
    {
      public:
        DynamicFEPointEvaluationImpl(NonMatching::MappingInfo<dim> &mapping,
                                     const FiniteElement<dim> &fe,
                                     const unsigned int        first_selected_component)
          : DynamicFEPointEvaluation<dim>(first_selected_component, n_components),
            evaluation(mapping, fe, first_selected_component)
        {}

        void evaluate(const ArrayView<double> &solution_values,
                      const EvaluationFlags::EvaluationFlags flags) override
        {
          evaluation.evaluate(solution_values, flags);
        }

        small_vector<double>
        get_value(const unsigned int evaluation_point) const override
        {
          const Tensor<1,n_components> x = convert::to_tensor<n_components>(evaluation.get_value(evaluation_point));
          small_vector<double> result (n_components);
          for (int c=0; c<n_components; ++c)
            result[c] = x[c];
          return result;
        }

        void
        get_value(const unsigned int evaluation_point,
                  const ArrayView<double> &solution) const override
        {
          Assert(solution.size() == n_components, ExcMessage("The size of the solution vector does not match the number of components."));

          const Tensor<1,n_components> x = convert::to_tensor<n_components>(evaluation.get_value(evaluation_point));
          for (int c=0; c<n_components; ++c)
            solution[c] = x[c];
        }

        small_vector<Tensor<1,dim>>
        get_gradient(const unsigned int evaluation_point) const override
        {
          const Tensor<1,n_components,Tensor<1,dim>> x = convert::to_tensor2<dim,n_components>(evaluation.get_gradient(evaluation_point));
          small_vector<Tensor<1,dim>> result (n_components);
          for (int c=0; c<n_components; ++c)
            result[c] = x[c];
          return result;
        }

        void
        get_gradient(const unsigned int evaluation_point,
                     const ArrayView<Tensor<1,dim>> &gradients) const override
        {
          Assert(gradients.size() == n_components, ExcMessage("The size of the gradient vector does not match the number of components."));

          const Tensor<1,n_components,Tensor<1,dim>> x = convert::to_tensor2<dim,n_components>(evaluation.get_gradient(evaluation_point));
          for (int c=0; c<n_components; ++c)
            gradients[c] = x[c];
        }

      private:
        FEPointEvaluation<n_components, dim, dim, double> evaluation;
    };




    template <int dim>
    static std::unique_ptr<DynamicFEPointEvaluation<dim>> make(NonMatching::MappingInfo<dim> &mapping,
                                                                const FiniteElement<dim> &fe,
                                                                const unsigned int        first_selected_component,
                                                                int n_fields)
    {
      switch (n_fields)
        {
          case 1:
            return std::make_unique<DynamicFEPointEvaluationImpl<dim,1>>(mapping, fe, first_selected_component);
          case 2:
            return std::make_unique<DynamicFEPointEvaluationImpl<dim,2>>(mapping, fe, first_selected_component);
          case 3:
            return std::make_unique<DynamicFEPointEvaluationImpl<dim,3>>(mapping, fe, first_selected_component);
          case 4:
            return std::make_unique<DynamicFEPointEvaluationImpl<dim,4>>(mapping, fe, first_selected_component);
          case 5:
            return std::make_unique<DynamicFEPointEvaluationImpl<dim,5>>(mapping, fe, first_selected_component);
          case 6:
            return std::make_unique<DynamicFEPointEvaluationImpl<dim,6>>(mapping, fe, first_selected_component);
          case 7:
            return std::make_unique<DynamicFEPointEvaluationImpl<dim,7>>(mapping, fe, first_selected_component);
          case 8:
            return std::make_unique<DynamicFEPointEvaluationImpl<dim,8>>(mapping, fe, first_selected_component);
          case 9:
            return std::make_unique<DynamicFEPointEvaluationImpl<dim,9>>(mapping, fe, first_selected_component);
          case 10:
            return std::make_unique<DynamicFEPointEvaluationImpl<dim,10>>(mapping, fe, first_selected_component);

          default:
            Assert(false, ExcNotImplemented());
            return std::unique_ptr<DynamicFEPointEvaluation<dim>>();
        }

    }
  }



  template <int dim>
  SolutionEvaluator<dim>::SolutionEvaluator(const SimulatorAccess<dim> &simulator,
                                            const UpdateFlags update_flags)
    :
    mapping_info(simulator.get_mapping(),
                 update_flags),
    velocity(mapping_info,
             simulator.get_fe(),
             simulator.introspection().component_indices.velocities[0]),
    pressure(std::make_unique<FEPointEvaluation<1, dim>>(mapping_info,
                                                          simulator.get_fe(),
                                                          simulator.introspection().component_indices.pressure)),
    temperature(mapping_info,
                simulator.get_fe(),
                simulator.introspection().component_indices.temperature),
    melt_component_indices(),
    simulator_access(simulator)
  {
    // Create the evaluators for all compositional fields
    {
      const auto &component_indices = simulator_access.introspection().component_indices.compositional_fields;

      // We consider each group of consecutive compositions of the same type together. This is because it is more efficient
      // than evaluating each one individually.
      for (const unsigned int base_element_index : simulator.introspection().get_composition_base_element_indices())
        {
          std::vector<unsigned int> indices = simulator.introspection().get_compositional_field_indices_with_base_element(base_element_index);

          // We can evaluate at most N at a time, if we have more than that of the same type, we tackle
          // them in groups of N:
          const unsigned int N = 10;
          while (indices.size()>N)
            {
              compositions.emplace_back(internal::make<dim>(mapping_info,
                                                            simulator_access.get_fe(),
                                                            component_indices[indices[0]],
                                                            N
                                                           ));

              indices.erase(indices.begin(), indices.begin() + N);
            }

          compositions.emplace_back(internal::make<dim>(mapping_info,
                                                        simulator_access.get_fe(),
                                                        component_indices[indices[0]],
                                                        indices.size()
                                                       ));
        }
    }

    // The FE_DGP pressure element used in locally conservative discretization is not
    // supported by the fast path of FEPointEvaluation. Replace with slow path.
    if (simulator_access.get_parameters().use_locally_conservative_discretization == true)
      pressure = std::make_unique<FEPointEvaluation<1, dim>>(simulator_access.get_mapping(),
                                                              simulator_access.get_fe(),
                                                              update_flags,
                                                              simulator.introspection().component_indices.pressure);

    // Create the melt evaluators, but only if we use melt transport in the model
    if (simulator_access.include_melt_transport())
      {
        // Store the melt component indices to avoid repeated string lookups later on
        melt_component_indices[0] = simulator_access.introspection().variable("fluid velocity").first_component_index;
        melt_component_indices[1] = simulator_access.introspection().variable("fluid pressure").first_component_index;
        melt_component_indices[2] = simulator_access.introspection().variable("compaction pressure").first_component_index;

        fluid_velocity = std::make_unique<FEPointEvaluation<dim, dim>>(mapping_info,
                                                                        simulator_access.get_fe(),
                                                                        melt_component_indices[0]);
        if (simulator_access.get_parameters().use_locally_conservative_discretization == false)
          fluid_pressure = std::make_unique<FEPointEvaluation<1, dim>>(mapping_info,
                                                                        simulator_access.get_fe(),
                                                                        melt_component_indices[1]);
        else
          {
            fluid_pressure = std::make_unique<FEPointEvaluation<1, dim>>(simulator_access.get_mapping(),
                                                                          simulator_access.get_fe(),
                                                                          update_flags,
                                                                          melt_component_indices[1]);
          }

        if (simulator_access.get_melt_handler().melt_parameters.use_discontinuous_p_c == false)
          compaction_pressure = std::make_unique<FEPointEvaluation<1, dim>>(mapping_info,
                                                                             simulator_access.get_fe(),
                                                                             melt_component_indices[2]);
        else
          compaction_pressure = std::make_unique<FEPointEvaluation<1, dim>>(simulator_access.get_mapping(),
                                                                             simulator_access.get_fe(),
                                                                             update_flags,
                                                                             melt_component_indices[2]);


      }
  }



  template <int dim>
  void
  SolutionEvaluator<dim>::reinit(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                 const ArrayView<Point<dim>> &positions)
  {
    mapping_info.reinit(cell,positions);

    if (simulator_access.get_parameters().use_locally_conservative_discretization == true)
      {
        pressure->reinit(cell, positions);

        if (simulator_access.include_melt_transport())
          {
            fluid_pressure->reinit (cell, positions);
          }
      }

    if (simulator_access.include_melt_transport()
        && simulator_access.get_melt_handler().melt_parameters.use_discontinuous_p_c == true)
      {
        compaction_pressure->reinit (cell, positions);
      }
  }



  template <int dim>
  void
  SolutionEvaluator<dim>::evaluate(const ArrayView<double> &solution_values,
                                   const std::vector<EvaluationFlags::EvaluationFlags> &evaluation_flags)
  {
    const auto &component_indices = simulator_access.introspection().component_indices;

    for (unsigned int i=1; i<dim; ++i)
      Assert(evaluation_flags[component_indices.velocities[i]] == evaluation_flags[component_indices.velocities[0]],
             ExcMessage("The evaluation flags for the velocity components are not the same."));

    for (unsigned int i=1; i<simulator_access.n_compositional_fields(); ++i)
      Assert(evaluation_flags[component_indices.compositional_fields[i]] == evaluation_flags[component_indices.compositional_fields[0]],
             ExcMessage("The evaluation flags for the compositional fields are not the same."));

    if (evaluation_flags[component_indices.velocities[0]] & EvaluationFlags::values ||
        evaluation_flags[component_indices.velocities[0]] & EvaluationFlags::gradients)
      velocity.evaluate (solution_values, evaluation_flags[component_indices.velocities[0]]);

    if (evaluation_flags[component_indices.pressure] & EvaluationFlags::values ||
        evaluation_flags[component_indices.pressure] & EvaluationFlags::gradients)
      pressure->evaluate (solution_values, evaluation_flags[component_indices.pressure]);

    if (evaluation_flags[component_indices.temperature] & EvaluationFlags::values ||
        evaluation_flags[component_indices.temperature] & EvaluationFlags::gradients)
      temperature.evaluate (solution_values, evaluation_flags[component_indices.temperature]);

    if (component_indices.compositional_fields.size() > 0)
      {
        if (evaluation_flags[component_indices.compositional_fields[0]] & EvaluationFlags::values ||
            evaluation_flags[component_indices.compositional_fields[0]] & EvaluationFlags::gradients)
          for (auto &eval: compositions)
            eval->evaluate (solution_values, evaluation_flags[eval->get_first_component()]);
      }

    if (simulator_access.include_melt_transport())
      {
        if (evaluation_flags[melt_component_indices[0]] & EvaluationFlags::values ||
            evaluation_flags[melt_component_indices[0]] & EvaluationFlags::gradients)
          fluid_velocity->evaluate (solution_values, evaluation_flags[melt_component_indices[0]]);
        if (evaluation_flags[melt_component_indices[1]] & EvaluationFlags::values ||
            evaluation_flags[melt_component_indices[1]] & EvaluationFlags::gradients)
          fluid_pressure->evaluate (solution_values, evaluation_flags[melt_component_indices[1]]);
        if (evaluation_flags[melt_component_indices[2]] & EvaluationFlags::values ||
            evaluation_flags[melt_component_indices[2]] & EvaluationFlags::gradients)
          compaction_pressure->evaluate (solution_values, evaluation_flags[melt_component_indices[2]]);
      }
  }



  namespace
  {
    template <int n_compositional_fields>
    double
    get_value(const Tensor<1,n_compositional_fields> &solution,
              const unsigned int component_index)
    {
      AssertIndexRange(component_index, n_compositional_fields);
      return solution[component_index];
    }

    template <int n_compositional_fields>
    double
    get_value(const double &solution,
              const unsigned int component_index)
    {
      (void) component_index;
      AssertIndexRange(component_index, 1);
      return solution;
    }

    template <int dim, int n_compositional_fields>
    Tensor<1,dim>
    get_gradient(const Tensor<1,n_compositional_fields,Tensor<1,dim>> &gradient,
                 const unsigned int component_index)
    {
      AssertIndexRange(component_index, n_compositional_fields);
      return gradient[component_index];
    }


    template <int dim, int n_compositional_fields>
    Tensor<1,dim>
    get_gradient(const Tensor<1,dim> &gradient,
                 const unsigned int component_index)
    {
      (void) component_index;
      AssertIndexRange(component_index, 1);
      return gradient;
    }
  }


  template <int dim>
  void
  SolutionEvaluator<dim>::get_solution(const unsigned int evaluation_point,
                                       const ArrayView<double> &solution,
                                       const std::vector<EvaluationFlags::EvaluationFlags> &evaluation_flags) const
  {
    Assert(solution.size() == simulator_access.introspection().n_components,
           ExcDimensionMismatch(solution.size(), simulator_access.introspection().n_components));

    const auto &component_indices = simulator_access.introspection().component_indices;

    if (evaluation_flags[component_indices.velocities[0]] & EvaluationFlags::values)
      {
        const Tensor<1,dim> velocity_value = velocity.get_value(evaluation_point);
        for (unsigned int j=0; j<dim; ++j)
          solution[component_indices.velocities[j]] = velocity_value[j];
      }


    if (evaluation_flags[component_indices.pressure] & EvaluationFlags::values)
      solution[component_indices.pressure] = pressure->get_value(evaluation_point);

    if (evaluation_flags[component_indices.temperature] & EvaluationFlags::values)
      solution[component_indices.temperature] = temperature.get_value(evaluation_point);

    for (const auto &eval : compositions)
      {
        const unsigned int start_index = eval->get_first_component();

        if (evaluation_flags[start_index] & EvaluationFlags::values)
          {
            const unsigned int n_components = eval->get_n_components();
            eval->get_value(evaluation_point,
            {&solution[start_index],n_components});
          }
      }

    if (simulator_access.include_melt_transport())
      {
        if (evaluation_flags[melt_component_indices[0]] & EvaluationFlags::values)
          {
            const Tensor<1,dim> fluid_velocity_value = fluid_velocity->get_value(evaluation_point);
            for (unsigned int j=0; j<dim; ++j)
              solution[melt_component_indices[0]+j] = fluid_velocity_value[j];
          }

        if (evaluation_flags[melt_component_indices[1]] & EvaluationFlags::values)
          solution[melt_component_indices[1]] = fluid_pressure->get_value(evaluation_point);

        if (evaluation_flags[melt_component_indices[2]] & EvaluationFlags::values)
          solution[melt_component_indices[2]] = compaction_pressure->get_value(evaluation_point);
      }
  }



  template <int dim>
  void
  SolutionEvaluator<dim>::get_gradients(const unsigned int evaluation_point,
                                        const ArrayView<Tensor<1,dim>> &gradients,
                                        const std::vector<EvaluationFlags::EvaluationFlags> &evaluation_flags) const
  {
    Assert(gradients.size() == simulator_access.introspection().n_components,
           ExcDimensionMismatch(gradients.size(), simulator_access.introspection().n_components));

    const auto &component_indices = simulator_access.introspection().component_indices;

    if (evaluation_flags[component_indices.velocities[0]] & EvaluationFlags::gradients)
      {
        const Tensor<2,dim> velocity_gradient = velocity.get_gradient(evaluation_point);
        for (unsigned int j=0; j<dim; ++j)
          gradients[component_indices.velocities[j]] = velocity_gradient[j];
      }

    if (evaluation_flags[component_indices.pressure] & EvaluationFlags::gradients)
      gradients[component_indices.pressure] = pressure->get_gradient(evaluation_point);

    if (evaluation_flags[component_indices.temperature] & EvaluationFlags::gradients)
      gradients[component_indices.temperature] = temperature.get_gradient(evaluation_point);

    for (const auto &eval : compositions)
      {
        const unsigned int start_index = eval->get_first_component();
        const unsigned int n_components = eval->get_n_components();

        if (evaluation_flags[start_index] & EvaluationFlags::gradients)
          eval->get_gradient(evaluation_point,
          {&gradients[start_index],n_components});
      }

    if (simulator_access.include_melt_transport())
      {
        if (evaluation_flags[melt_component_indices[0]] & EvaluationFlags::gradients)
          {
            const Tensor<2,dim> fluid_velocity_gradient = fluid_velocity->get_gradient(evaluation_point);
            for (unsigned int j=0; j<dim; ++j)
              gradients[melt_component_indices[0]+j] = fluid_velocity_gradient[j];
          }

        if (evaluation_flags[melt_component_indices[1]] & EvaluationFlags::gradients)
          gradients[melt_component_indices[1]] = fluid_pressure->get_gradient(evaluation_point);

        if (evaluation_flags[melt_component_indices[2]] & EvaluationFlags::gradients)
          gradients[melt_component_indices[2]] = compaction_pressure->get_gradient(evaluation_point);
      }
  }



  template <int dim>
  FEPointEvaluation<dim, dim> &
  SolutionEvaluator<dim>::get_velocity_or_fluid_velocity_evaluator(const bool use_fluid_velocity)
  {
    if (use_fluid_velocity)
      return *fluid_velocity;
    else
      return velocity;

    return velocity;
  }



  template <int dim>
  NonMatching::MappingInfo<dim> &
  SolutionEvaluator<dim>::get_mapping_info()
  {
    return mapping_info;
  }



  template <int dim>
  unsigned int
  SolutionEvaluator<dim>::n_components() const
  {
    return simulator_access.introspection().n_components;
  }



  template <int dim>
  std::unique_ptr<SolutionEvaluator<dim>>
  construct_solution_evaluator (const SimulatorAccess<dim> &simulator_access,
                                const UpdateFlags update_flags)
  {
    return std::make_unique<SolutionEvaluator<dim>>(simulator_access, update_flags);
  }
}


// explicit instantiation of the functions we implement in this file
namespace aspect
{
#define INSTANTIATE(dim) \
  template class SolutionEvaluator<dim>; \
  template std::unique_ptr<SolutionEvaluator<dim>> construct_solution_evaluator (const SimulatorAccess<dim> &simulator_access, \
                                                                                  const UpdateFlags update_flags);

  ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}
