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

#include <aspect/adjoint/objective_functional.h>
#include <aspect/adjoint/kernel_calculator.h>
#include <aspect/utilities.h>

#include <deal.II/base/utilities.h>

#include <set>
#include <tuple>

namespace aspect
{
  namespace Adjoint
  {
    namespace
    {
      std::tuple<aspect::internal::Plugins::UnusablePluginList,
          aspect::internal::Plugins::UnusablePluginList,
          aspect::internal::Plugins::PluginList<ObjectiveFunctional<2>>,
          aspect::internal::Plugins::PluginList<ObjectiveFunctional<3>>> registered_objectives;
    }



    template <int dim>
    std::set<std::string>
    ObjectiveFunctional<dim>::required_forward_postprocessors() const
    {
      return {};
    }



    template <int dim>
    void
    ObjectiveManager<dim>::register_objective_functional(const std::string &name,
                                                         const std::string &description,
                                                         void (*declare_parameters_function) (ParameterHandler &),
                                                         std::unique_ptr<ObjectiveFunctional<dim>> (*factory_function) ())
    {
      std::get<dim>(registered_objectives).register_plugin(name,
                                                           description,
                                                           declare_parameters_function,
                                                           factory_function);
    }



    template <int dim>
    void
    ObjectiveManager<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Adjoint");
      {
        prm.declare_entry("List of objectives", "",
                          Patterns::MultipleSelection(std::get<dim>(registered_objectives).get_pattern_of_names()),
                          "Comma separated list of adjoint objective functionals.\n\n"
                          "The following objective functionals are available:\n\n"
                          +
                          std::get<dim>(registered_objectives).get_description_string());
      }
      prm.leave_subsection();

      std::get<dim>(registered_objectives).declare_parameters(prm);
    }



    template <int dim>
    void
    ObjectiveManager<dim>::parse_parameters(ParameterHandler &prm)
    {
      this->plugin_objects.clear();
      this->plugin_names.clear();

      prm.enter_subsection("Adjoint");
      {
        this->plugin_names =
          Utilities::split_string_list(prm.get("List of objectives"));
      }
      prm.leave_subsection();

      AssertThrow(Utilities::has_unique_entries(this->plugin_names),
                  ExcMessage("The list of strings for the parameter "
                             "'Adjoint/List of objectives' contains entries more than once. "
                             "This is not allowed. Please check your parameter file."));

      for (const std::string &objective_name : this->plugin_names)
        {
          this->plugin_objects.emplace_back(std::get<dim>(registered_objectives)
                                            .create_plugin(objective_name,
                                                           "Adjoint::List of objectives"));
          this->plugin_objects.back()->initialize_simulator(this->get_simulator());
          this->plugin_objects.back()->parse_parameters(prm);
          this->plugin_objects.back()->initialize();
        }
    }



    template <int dim>
    std::set<std::string>
    ObjectiveManager<dim>::required_forward_postprocessors() const
    {
      std::set<std::string> required_postprocessors;

      for (const auto &objective : this->get_active_plugins())
        {
          const std::set<std::string> objective_required_postprocessors =
            objective->required_forward_postprocessors();

          required_postprocessors.insert(objective_required_postprocessors.begin(),
                                         objective_required_postprocessors.end());
        }

      return required_postprocessors;
    }



    template <int dim>
    std::vector<std::unique_ptr<ObjectiveResult<dim>>>
    ObjectiveManager<dim>::evaluate_and_assemble_right_hand_sides() const
    {
      const auto &objective_functionals = this->get_active_plugins();
      const auto &objective_names = this->get_active_plugin_names();

      AssertThrow(objective_functionals.size() == objective_names.size(),
                  ExcInternalError());

      std::vector<std::unique_ptr<ObjectiveResult<dim>>> objective_results;
      objective_results.reserve(objective_names.size());

      auto objective = objective_functionals.begin();
      for (unsigned int objective_index = 0;
           objective_index < objective_names.size();
           ++objective_index, ++objective)
        {
          auto result = std::make_unique<ObjectiveResult<dim>>();
          result->objective_name = objective_names[objective_index];
          result->value = (*objective)->evaluate();
          result->rhs.reinit(this->introspection().index_sets.system_partitioning,
                             this->get_mpi_communicator());
          result->rhs = 0.0;

          (*objective)->assemble_adjoint_rhs(result->rhs);
          objective_results.push_back(std::move(result));
        }

      return objective_results;
    }



    template <int dim>
    void
    ObjectiveManager<dim>::add_direct_kernel_contributions(KernelRepository<dim> &kernels) const
    {
      const auto &objective_functionals = this->get_active_plugins();
      const auto &objective_names = this->get_active_plugin_names();

      AssertThrow(objective_functionals.size() == objective_names.size(),
                  ExcInternalError());

      auto objective = objective_functionals.begin();
      for (unsigned int objective_index = 0;
           objective_index < objective_names.size();
           ++objective_index, ++objective)
        (*objective)->add_direct_kernel_contributions(objective_names[objective_index],
                                                      kernels);
    }
  }
}

namespace aspect
{
#define INSTANTIATE(dim) \
  template std::set<std::string> Adjoint::ObjectiveFunctional<dim>::required_forward_postprocessors() const; \
  template class Adjoint::ObjectiveManager<dim>;

  ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}
