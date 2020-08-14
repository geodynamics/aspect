/*
  Copyright (C) 2020 by the authors of the ASPECT code.

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


#include <aspect/material_model/rheology/constant_viscosity_prefactors.h>
#include <aspect/utilities.h>

#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/parameter_handler.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace Rheology
    {
      template <int dim>
      ConstantViscosityPrefactors<dim>::ConstantViscosityPrefactors ()
      {}



      template <int dim>
      double
      ConstantViscosityPrefactors<dim>::compute_viscosity (const double base_viscosity,
                                                           const unsigned int composition_index) const
      {
        return base_viscosity * constant_viscosity_prefactors[composition_index];
      }



      template <int dim>
      void
      ConstantViscosityPrefactors<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Constant viscosity prefactors", "1.0",
                           Patterns::List(Patterns::Double(0)),
                           "List of constant viscosity prefactors (i.e., multiplicative factors) "
                           "for background material and compositional fields, for a total of N+1 "
                           "where N is the number of compositional fields. Units: none.");
      }



      template <int dim>
      void
      ConstantViscosityPrefactors<dim>::parse_parameters (ParameterHandler &prm)
      {
        // increment by one for background:
        const unsigned int n_fields = this->n_compositional_fields() + 1;

        constant_viscosity_prefactors = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Constant viscosity prefactors"))),
                                                                                n_fields,
                                                                                "Constant viscosity prefactors");
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
    template class ConstantViscosityPrefactors<dim>; \
  }

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
