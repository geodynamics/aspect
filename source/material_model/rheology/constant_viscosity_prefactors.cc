/*
  Copyright (C) 2020 - 2024 by the authors of the ASPECT code.

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
        = default;



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
                           "where N is the number of all compositional fields or only those "
                           "corresponding to chemical compositions. Units: none.");
      }



      template <int dim>
      void
      ConstantViscosityPrefactors<dim>::parse_parameters (ParameterHandler &prm)
      {
        // Retrieve the list of composition names
        std::vector<std::string> compositional_field_names = this->introspection().get_composition_names();

        // Retrieve the list of names of fields that represent chemical compositions, and not, e.g.,
        // plastic strain
        std::vector<std::string> chemical_field_names = this->introspection().chemical_composition_field_names();

        // Establish that a background field is required here
        compositional_field_names.insert(compositional_field_names.begin(), "background");
        chemical_field_names.insert(chemical_field_names.begin(),"background");

        Utilities::MapParsing::Options options(chemical_field_names, "Constant viscosity prefactors");
        options.list_of_allowed_keys = compositional_field_names;

        constant_viscosity_prefactors = Utilities::MapParsing::parse_map_to_double_array (prm.get("Constant viscosity prefactors"),
                                        options);
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
