/*
  Copyright (C) 2011 - 2021 by the authors of the ASPECT code.

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


#include <aspect/material_model/equation_of_state/linearized_incompressible.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace EquationOfState
    {
      template <int dim>
      void
      LinearizedIncompressible<dim>::
      evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
               const unsigned int q,
               MaterialModel::EquationOfStateOutputs<dim> &out) const
      {

        Assert(maximum_number_of_compositions+1 >= out.densities.size(),
               ExcMessage("Error: You are trying to evaluate the equation of state with "
                          + Utilities::to_string(out.densities.size()-1) +
                          " compositional fields, which is larger than "
                          + Utilities::to_string(maximum_number_of_compositions) +
                          ", the number of fields the equation of state was set up with."));


        for (unsigned int c=0; c < out.densities.size(); ++c)
          {
            out.densities[c] = reference_rho * (1 - thermal_alpha * (in.temperature[q] - reference_T));
            if (c>0)
              out.densities[c] += compositional_delta_rhos[c-1];

            out.thermal_expansion_coefficients[c] = thermal_alpha;
            out.specific_heat_capacities[c] = reference_specific_heat;
            out.compressibilities[c] = 0.0;
            out.entropy_derivative_pressure[c] = 0.0;
            out.entropy_derivative_temperature[c] = 0.0;
          }
      }



      template <int dim>
      bool
      LinearizedIncompressible<dim>::
      is_compressible () const
      {
        return false;
      }



      template <int dim>
      void
      LinearizedIncompressible<dim>::declare_parameters (ParameterHandler &prm,
                                                         const unsigned int n_compositions)
      {
        prm.declare_entry ("Reference density", "3300.",
                           Patterns::Double (0.),
                           "Reference density $\\rho_0$. "
                           "Units: \\si{\\kilogram\\per\\meter\\cubed}.");
        prm.declare_entry ("Reference temperature", "293.",
                           Patterns::Double (0.),
                           "The reference temperature $T_0$. The reference temperature is used "
                           "in both the density and viscosity formulas. Units: \\si{\\kelvin}.");
        prm.declare_entry ("Reference specific heat", "1250.",
                           Patterns::Double (0.),
                           "The value of the specific heat $C_p$. "
                           "Units: \\si{\\joule\\per\\kelvin\\per\\kilogram}.");
        prm.declare_entry ("Thermal expansion coefficient", "2e-5",
                           Patterns::Double (0.),
                           "The value of the thermal expansion coefficient $\\alpha$. "
                           "Units: \\si{\\per\\kelvin}.");
        if (n_compositions > 0)
          prm.declare_entry ("Density differential for compositional field 1", "0.",
                             Patterns::Double(),
                             "If compositional fields are used, then one would frequently want "
                             "to make the density depend on these fields. In this simple material "
                             "model, we make the following assumptions: if no compositional fields "
                             "are used in the current simulation, then the density is simply the usual "
                             "one with its linear dependence on the temperature. If there are compositional "
                             "fields, then the material model determines how many of them influence the density. "
                             "The composition-dependence adds a term of the kind $+\\Delta \\rho \\; c_1(\\mathbf x)$. "
                             "This parameter describes the value of $\\Delta \\rho$. "
                             "Units: \\si{\\kilogram\\per\\meter\\cubed}/unit change in composition.");
        if (n_compositions > 1)
          prm.declare_entry ("Density differential for compositional field 2", "0.",
                             Patterns::Double(),
                             "If compositional fields are used, then one would frequently want "
                             "to make the density depend on these fields. In this simple material "
                             "model, we make the following assumptions: if no compositional fields "
                             "are used in the current simulation, then the density is simply the usual "
                             "one with its linear dependence on the temperature. If there are compositional "
                             "fields, then the material model determines how many of them influence the density. "
                             "The composition-dependence adds a term of the kind $+\\Delta \\rho \\; c_2(\\mathbf x)$. "
                             "This parameter describes the value of $\\Delta \\rho$. "
                             "Units: \\si{\\kilogram\\per\\meter\\cubed}/unit change in composition.");
      }



      template <int dim>
      void
      LinearizedIncompressible<dim>::parse_parameters (ParameterHandler &prm,
                                                       const unsigned int n_compositions)
      {
        reference_rho              = prm.get_double ("Reference density");
        reference_T                = prm.get_double ("Reference temperature");
        reference_specific_heat    = prm.get_double ("Reference specific heat");
        thermal_alpha              = prm.get_double ("Thermal expansion coefficient");

        maximum_number_of_compositions = n_compositions;
        compositional_delta_rhos.resize(maximum_number_of_compositions);
        if (maximum_number_of_compositions > 0)
          compositional_delta_rhos[0]    = prm.get_double ("Density differential for compositional field 1");
        if (maximum_number_of_compositions > 1)
          compositional_delta_rhos[1]    = prm.get_double ("Density differential for compositional field 2");
      }

    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    namespace EquationOfState
    {
#define INSTANTIATE(dim) \
  template class LinearizedIncompressible<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
    }
  }
}
