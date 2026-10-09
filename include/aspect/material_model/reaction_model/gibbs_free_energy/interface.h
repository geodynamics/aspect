/*
  Copyright (C) 2011 - 2026 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_reaction_model_gibbs_free_energy_interface_h
#define _aspect_material_model_reaction_model_gibbs_free_energy_interface_h

#include <aspect/plugins.h>
#include <deal.II/base/parameter_handler.h>

#include <memory>
#include <string>

namespace aspect
{
  namespace MaterialModel
  {
    namespace ReactionModel
    {
      struct GibbsFreeEnergyInputs
      {
        /** Reference temperature at which the reference Gibbs free-energy difference is evaluated. Units: K. */
        double reference_temperature = 0.0;
        /** Reference pressure at which the reference Gibbs free-energy difference is evaluated. Units: Pa. */
        double reference_pressure = 0.0;
        /** Reference value of G_B-G_A. Units: J/mol. */
        double reference_delta_gibbs_free_energy = 0.0;
        /** Difference S_B-S_A in molar entropy. Units: J/(mol K). */
        double delta_entropy = 0.0;
        /** Difference V_B-V_A in molar volume. Units: m^3/mol. */
        double delta_volume = 0.0;
      };

      namespace GibbsFreeEnergyModel
      {
        /**
         * Interface for models that compute the raw Gibbs free-energy
         * difference delta_G = G_B-G_A between two phases A and B.
         * The phase names define the ordering of the difference, but not the
         * direction in which the reaction currently proceeds. A negative
         * delta_G favors A -> B, while a positive delta_G favors B -> A.
         */
        template <int dim>
        class Interface
        {
          public:
            virtual ~Interface() = default;

            /**
             * Compute the raw Gibbs free-energy difference G_B-G_A. Phases A
             * and B specify the ordering of the difference, not which phase is
             * currently consumed. The sign of the returned value determines
             * the favored direction: a negative value favors A -> B, a
             * positive value favors B -> A, and zero represents equilibrium.
             */
            virtual double
            compute_delta_gibbs_free_energy(const double temperature,
                                            const double pressure,
                                            const GibbsFreeEnergyInputs &reference_state,
                                            const unsigned int reaction_index) const = 0;

            static void declare_parameters(ParameterHandler &prm);
            virtual void parse_parameters(ParameterHandler &prm, const unsigned int n_reactions) = 0;
        };

        template <int dim>
        std::unique_ptr<Interface<dim>> create_gibbs_free_energy_model(const std::string &model_name);

        template <int dim>
        using PluginList = aspect::internal::Plugins::PluginList<Interface<dim>>;

#define ASPECT_REGISTER_GIBBS_FREE_ENERGY_MODEL(classname, name, description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_GIBBS_FREE_ENERGY_MODEL_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::MaterialModel::ReactionModel::GibbsFreeEnergyModel::Interface<2>, classname<2>> \
    dummy_ ## classname ## _2d(&aspect::MaterialModel::ReactionModel::GibbsFreeEnergyModel::PluginList<2>::register_plugin, name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::MaterialModel::ReactionModel::GibbsFreeEnergyModel::Interface<3>, classname<3>> \
    dummy_ ## classname ## _3d(&aspect::MaterialModel::ReactionModel::GibbsFreeEnergyModel::PluginList<3>::register_plugin, name, description); \
  }
      }
    }
  }
}

#endif
