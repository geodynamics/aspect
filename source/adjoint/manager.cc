/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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

#include <aspect/adjoint/manager.h>
#include <aspect/postprocess/dynamic_topography.h>
#include <aspect/simulator.h>

#include <deal.II/base/table_handler.h>
#include <deal.II/base/utilities.h>

#include <algorithm>
#include <set>
#include <string>

namespace aspect
{
  namespace Adjoint
  {
    template <int dim>
    Manager<dim>::Manager()
      : simulator(nullptr)
    {}



    template <int dim>
    void
    Manager<dim>::declare_parameters (ParameterHandler &prm)
    {
      ObjectiveManager<dim>::declare_parameters (prm);
    }



    template <int dim>
    void
    Manager<dim>::initialize_simulator(Simulator<dim> &simulator_object)
    {
      simulator = &simulator_object;

      // Initialize the read-only SimulatorAccess interface explicitly. Calling
      // this->initialize_simulator() here would recurse into this overload.
      this->SimulatorAccess<dim>::initialize_simulator(simulator_object);
      objective_manager.initialize_simulator(simulator_object);
    }



    template <int dim>
    void
    Manager<dim>::initialize_simulator(const Simulator<dim> &)
    {
      // This override exists only to reject initialization through the generic
      // SimulatorAccess interface, which cannot provide the mutable simulator
      // access required by the adjoint workflow.
      AssertThrow(false,
                  ExcMessage("The adjoint manager requires mutable simulator access. "
                             "Call initialize_simulator(Simulator<dim> &) instead."));
    }



    template <int dim>
    void
    Manager<dim>::parse_parameters(ParameterHandler &prm)
    {
      Assert(simulator != nullptr, ExcInternalError());

      objective_manager.parse_parameters(prm);
    }



    template <int dim>
    void
    Manager<dim>::initialize()
    {}



    template <int dim>
    const KernelRepository<dim> &
    Manager<dim>::get_kernels() const
    {
      return kernels;
    }



    template <int dim>
    std::map<std::string, double>
    Manager<dim>::get_objective_values() const
    {
      std::map<std::string, double> values;

      for (const auto &objective_result : objective_results)
        values[objective_result->objective_name] = objective_result->value;

      return values;
    }



    template <int dim>
    void
    Manager<dim>::validate_instantaneous_stokes_setup() const
    {
      using MassConservation = typename Parameters<dim>::Formulation::MassConservation;

      const std::vector<std::string> &objective_names =
        objective_manager.get_active_plugin_names();
      const bool uses_dynamic_topography_objective =
        std::find(objective_names.begin(),
                  objective_names.end(),
                  "dynamic topography") != objective_names.end();

      if (uses_dynamic_topography_objective)
        AssertThrow(this->get_postprocess_manager().template has_matching_active_plugin<Postprocess::DynamicTopography<dim>>(),
                    ExcMessage("The dynamic topography adjoint objective requires the dynamic topography postprocessor. "
                               "Add dynamic topography to Postprocess/List of postprocessors."));

      AssertThrow(this->get_parameters().formulation_mass_conservation == MassConservation::incompressible,
                  ExcMessage("The dynamic-topography adjoint currently supports only the incompressible "
                             "mass conservation formulation. Compressible, reference-density, and projected-density "
                             "formulations need explicit adjoint physics terms before they can be used."));

      AssertThrow(this->get_material_model().is_compressible() == false,
                  ExcMessage("The dynamic-topography adjoint currently supports only incompressible material models."));
    }



    template <int dim>
    void
    Manager<dim>::prepare_required_forward_postprocessors()
    {
      const std::set<std::string> required_postprocessors =
        objective_manager.required_forward_postprocessors();

      if (required_postprocessors.empty())
        return;

      TableHandler temporary_statistics;

      // These postprocessors are executed early to prepare forward quantities
      // required by the objective functionals. They may be executed again during
      // the normal postprocessing phase after the adjoint kernels are available.
      simulator->postprocess_manager.execute_selected_postprocessors(required_postprocessors,
                                                                     temporary_statistics);
    }



    template <int dim>
    void
    Manager<dim>::solve_adjoint_states()
    {
      AssertThrow(objective_results.empty() == false,
                  ExcMessage("Adjoint right hand sides must be assembled before solving adjoint states."));

      adjoint_states.clear();

      const LinearAlgebra::BlockVector saved_solution = simulator->solution;
      const LinearAlgebra::BlockVector saved_linearization_point = simulator->current_linearization_point;
      const LinearAlgebra::BlockVector saved_system_rhs = simulator->system_rhs;
      const double saved_pressure_normalization_adjustment = simulator->last_pressure_normalization_adjustment;

      const unsigned int velocity_block_index = this->introspection().block_indices.velocities;
      const unsigned int pressure_block_index = this->introspection().block_indices.pressure;

      for (const auto &objective_result : objective_results)
        {
          auto adjoint_state = std::make_unique<AdjointState<dim>>();
          adjoint_state->objective_name = objective_result->objective_name;
          adjoint_state->solution.reinit(this->introspection().index_sets.system_partitioning,
                                         this->introspection().index_sets.system_relevant_partitioning,
                                         this->get_mpi_communicator());
          adjoint_state->solution = 0.0;

          simulator->system_rhs = objective_result->rhs;

          simulator->current_linearization_point.block(velocity_block_index) = 0.0;
          simulator->current_linearization_point.block(pressure_block_index) = 0.0;
          simulator->solution.block(velocity_block_index) = 0.0;
          simulator->solution.block(pressure_block_index) = 0.0;

          simulator->solve_stokes(adjoint_state->solution);
          adjoint_states.push_back(std::move(adjoint_state));
        }

      simulator->solution = saved_solution;
      simulator->current_linearization_point = saved_linearization_point;
      simulator->system_rhs = saved_system_rhs;
      simulator->last_pressure_normalization_adjustment = saved_pressure_normalization_adjustment;
    }



    template <int dim>
    void
    Manager<dim>::calculate_kernels()
    {
      KernelCalculator<dim> kernel_calculator(*simulator);
      kernels = kernel_calculator.calculate(objective_results, adjoint_states);

      objective_manager.add_direct_kernel_contributions(kernels);
    }



    template <int dim>
    void
    Manager<dim>::solve_instantaneous_stokes()
    {
      Assert(simulator != nullptr, ExcInternalError());

      validate_instantaneous_stokes_setup();

      simulator->assemble_and_solve_stokes();

      prepare_required_forward_postprocessors();
      objective_results = objective_manager.evaluate_and_assemble_right_hand_sides();
      solve_adjoint_states();
      calculate_kernels();
    }
  }
}

namespace aspect
{
#define INSTANTIATE(dim) \
  template class Adjoint::Manager<dim>;

  ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}
