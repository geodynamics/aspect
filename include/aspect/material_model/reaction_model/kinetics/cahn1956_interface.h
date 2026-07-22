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

#ifndef _aspect_material_model_reaction_model_kinetics_cahn1956_interface_h
#define _aspect_material_model_reaction_model_kinetics_cahn1956_interface_h

#include <aspect/plugins.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>
#include <deal.II/base/parameter_handler.h>

#include <memory>
#include <string>
#include <vector>

namespace aspect
{
  namespace MaterialModel
  {
    namespace ReactionModel
    {
      /**
       * Kinetic laws derived from Cahn (1956; https://doi.org/10.1016/0001-6160(56)90041-4) nucleation-and-growth theory in the limit of site
       * saturation (all nucleation sites consumed early) with Avrami exponent n = 1 (interface-controlled growth on grain boundaries). In this
       * limit the transformation kinetics reduce to a rate law that is first-order in the untransformed phase fraction and depends only on the
       * instantaneous local thermodynamic driving force and temperature/pressure (no explicit time or nucleation-rate dependence is needed).
       *
       * @ingroup ReactionModel
       */
      template <int dim>
      class Cahn1956Interface : public ::aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * Destructor.
           */
          virtual ~Cahn1956Interface() = default;

          /**
           * Net forward reaction rate dX_B/dt for the transformation A -> B at @p reaction_index.
           */
          virtual double
          net_forward_reaction_rate(const double temperature,
                                    const double pressure,
                                    const double delta_forward_gibbs_energy,
                                    const double cumulative_forward_reaction_progress,
                                    const unsigned int reaction_index) const = 0;

          /**
           * Declare kinetic parameters (e.g., activation energy, kinetic prefactor).
           */
          static void declare_parameters(ParameterHandler &prm);

          /**
           * Read model-specific kinetic parameters for the @p n_reactions.
           */
          virtual void parse_parameters(ParameterHandler &prm, const unsigned int n_reactions) = 0;
      };

      /**
       * Instantiate the reaction kinetics plugin named `model_name`.
       */
      template <int dim>
      std::unique_ptr<Cahn1956Interface<dim>> create_reaction_model(const std::string &model_name);

      template <int dim>
      using ReactionModelPluginList = internal::Plugins::PluginList<Cahn1956Interface<dim>>;

#define ASPECT_REGISTER_REACTION_MODEL(classname, name, description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_REACTION_MODEL_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::MaterialModel::ReactionModel::Cahn1956Interface<2>, classname<2>> \
    dummy_ ## classname ## _2d(&aspect::MaterialModel::ReactionModel::ReactionModelPluginList<2>::register_plugin, name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::MaterialModel::ReactionModel::Cahn1956Interface<3>, classname<3>> \
    dummy_ ## classname ## _3d(&aspect::MaterialModel::ReactionModel::ReactionModelPluginList<3>::register_plugin, name, description); \
  }
    }
  }
}

#endif
