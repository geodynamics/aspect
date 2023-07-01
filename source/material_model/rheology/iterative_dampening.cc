/*
  Copyright (C) 2019 - 2023 by the authors of the ASPECT code.

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


#include <aspect/material_model/rheology/iterative_dampening.h>

#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/parameter_handler.h>
#include <aspect/utilities.h>
#include <aspect/postprocess/particles.h>
#include <aspect/particle/property/interface.h>
#include <aspect/simulator.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace Rheology
    {
      template <int dim>
      IterativeDampening<dim>::IterativeDampening ()
        = default;


      template <int dim>
      void
      IterativeDampening<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Iterative viscosity dampening factor", "1.0",
                           Patterns::Double (0.),
                           "A dampening factor for the viscosity that controls the rate of change "
                           "between the viscosity calculated on the previous and current nonlinear "
                           "iteration. "
                           "Units: none.");
      }

      template <int dim>
      void
      IterativeDampening<dim>::parse_parameters (ParameterHandler &prm)
      {
        iterative_viscosity_dampening_factor = prm.get_double("Iterative viscosity dampening factor");
      }


      template <int dim>
      void
      IterativeDampening<dim>::
      fill_reaction_outputs (const MaterialModel::MaterialModelInputs<dim> &in,
                             const int i,
                             const double old_viscosity,
                             MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        if  (this->simulator_is_past_initialization() && in.requests_property(MaterialProperties::reaction_terms))
          {
            Assert(std::isfinite(in.strain_rate[i].norm()),
                   ExcMessage("Invalid strain_rate in the MaterialModelInputs. This is likely because it was "
                              "not filled by the caller."));

            // Index for viscosity compositional field
            const int vis_index = this->introspection().compositional_index_for_name("viscosity_field");

            out.reaction_terms[i][vis_index] = -old_viscosity + out.viscosities[i];
          }
      }

      template <int dim>
      double
      IterativeDampening<dim>::
      calculate_viscosity (const double old_viscosity,
                           const double new_viscosity) const
      {
        const double dampened_viscosity = std::pow(old_viscosity, iterative_viscosity_dampening_factor) *
                                          std::pow(new_viscosity, 1. - iterative_viscosity_dampening_factor);

        return dampened_viscosity;
      }
    }
  }
}
// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
#define INSTANTIATE(dim) \
  namespace Rheology \
  { \
    template class IterativeDampening<dim>; \
  }

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
