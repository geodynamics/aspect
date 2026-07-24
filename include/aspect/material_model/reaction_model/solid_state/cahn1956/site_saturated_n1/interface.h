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

#ifndef _aspect_material_model_reaction_model_solid_state_cahn1956_site_saturated_n1_interface_h
#define _aspect_material_model_reaction_model_solid_state_cahn1956_site_saturated_n1_interface_h

#include <aspect/plugins.h>
#include <aspect/simulator_access.h>
#include <deal.II/base/parameter_handler.h>

#include <memory>
#include <string>

namespace aspect
{
  namespace MaterialModel
  {
    namespace ReactionModel
    {
      namespace SolidState
      {
        /**
         * Kinetic laws derived from Cahn (1956) nucleation-and-growth theory in the limit of site saturation (all nucleation sites consumed early)
         * with Avrami exponent n = 1 (interface-controlled growth on grain boundaries). In this limit the transformation kinetics reduce to a rate
         * law that is first-order in the untransformed phase fraction and depends only on the instantaneous local thermodynamic driving force and
         * temperature/pressure (no explicit time or nucleation-rate dependence is needed).
         *
         * This is one specific family of solid-state reaction kinetics among several that ASPECT supports (see the sibling namespaces under
         * aspect::MaterialModel::ReactionModel::SolidState). Other families (e.g. a future JMAK-type formulation), or the general (non-site-saturated)
         * Cahn1956 model, have different physics and therefore different call signatures (each gets their own Interface class and plugin list).
         * Do not add virtual functions here that only make sense for those other families.
         */
        namespace Cahn1956
        {
          namespace SiteSaturatedN1
          {
            /**
             * Base class for objects describing the kinetics of a single two-phase transformation A -> B under site-saturated, n = 1 Cahn (1956) kinetics.
             * Implementations are agnostic to how dG was computed and to how the resulting rate is used downstream.
             *
             * @ingroup ReactionModel
             */
            template <int dim>
            class Interface : public ::aspect::SimulatorAccess<dim>
            {
              public:
                virtual ~Interface() = default;

                /**
                 * Reaction rate dX/dt for the forward transformation A -> B. See InterfaceControlledGrowth docs for the sign convention.
                 */
                virtual double
                reaction_rate(const double temperature, const double pressure, const double delta_gibbs_energy, const double reactant_phase_fraction) const = 0;

                /**
                 * Compute and return the dimensionless (or model-specific) Arrhenius/thermal activation factor isolated from the overall rate law.
                 */
                virtual double arrhenius_factor(const double temperature, const double pressure) const = 0;

                /**
                 * Compute and return the dimensionless (or model-specific) thermodynamic driving force factor isolated from the overall rate law.
                 */
                virtual double thermodynamic_factor(const double temperature, const double delta_gibbs_energy) const = 0;

                /**
                 * Declare this model's own parameters at the current subsection level.
                 */
                static void declare_parameters(ParameterHandler &prm);

                /**
                 * Read this model's own parameters at the current subsection level.
                 */
                virtual void parse_parameters(ParameterHandler &prm) = 0;
            };

            /**
             * Instantiate the registered ReactionModel plugin named `model_name`.
             */
            template <int dim>
            std::unique_ptr<Interface<dim>> create_reaction_model(const std::string &model_name);

            /**
             * Declare the parameters of every registered plugin in this family, so that a "Kinetics model" selector entry can pick among them later.
             * Call once per reaction-step subsection, before reading which model was actually selected.
             */
            template <int dim>
            void declare_parameters(ParameterHandler &prm);

            template <int dim>
            using ReactionModelPluginList = internal::Plugins::PluginList<Interface<dim>>;

#define ASPECT_REGISTER_CAHN1956_SITE_SATURATED_N1_REACTION_MODEL(classname, name, description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_CAHN1956_SITE_SATURATED_N1_REACTION_MODEL_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::MaterialModel::ReactionModel::SolidState::Cahn1956::SiteSaturatedN1::Interface<2>, classname<2>> \
    dummy_ ## classname ## _2d(&aspect::MaterialModel::ReactionModel::SolidState::Cahn1956::SiteSaturatedN1::ReactionModelPluginList<2>::register_plugin, name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::MaterialModel::ReactionModel::SolidState::Cahn1956::SiteSaturatedN1::Interface<3>, classname<3>> \
    dummy_ ## classname ## _3d(&aspect::MaterialModel::ReactionModel::SolidState::Cahn1956::SiteSaturatedN1::ReactionModelPluginList<3>::register_plugin, name, description); \
  }
          }
        }
      }
    }
  }
}

#endif
