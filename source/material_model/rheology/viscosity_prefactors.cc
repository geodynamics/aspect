/*
  Copyright (C) 2020 - 2023 by the authors of the ASPECT code.

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


#include <aspect/material_model/rheology/viscosity_prefactors.h>
#include <aspect/utilities.h>
#include <aspect/global.h>
#include <aspect/geometry_model/interface.h>

#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/parameter_handler.h>
#include <aspect/simulator_signals.h>
#include <aspect/parameters.h>


// TODO: Add constant_viscosity_prefactors to this file at a point when ASPECT is ready to be broken
namespace aspect
{
  namespace MaterialModel
  {
    namespace Rheology
    {
      template <int dim>
      ViscosityPrefactors<dim>::ViscosityPrefactors ()
        = default;



      template <int dim>
      double
      ViscosityPrefactors<dim>::compute_viscosity (const double base_viscosity,
                                                   const unsigned int composition_index) const
      {
        double time = this->get_time();
        if (this->convert_output_to_years())
          time /= year_in_seconds;

        const double factor = weakening_times[composition_index] <= time ? variable_viscosity_prefactors[composition_index] : 1.0;
        return base_viscosity * factor;
      }



      template <int dim>
      void
      ViscosityPrefactors<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Viscosity prefactors", "1.0",
                           Patterns::List(Patterns::Double(0)),
                           "List of viscosity prefactors (i.e., multiplicative factors) "
                           "for background material and compositional fields, for a total of N+1 "
                           "where N is the number of all compositional fields or only those "
                           "corresponding to chemical compositions. Units: none.");

        prm.declare_entry ("Weakening times", "1.0",
                           Patterns::List(Patterns::Double(0)),
                           "List of times to be specify the onset of weakening in years, or seconds if "
                           "the parameter 'Use years in output instead of seconds = false'. Before the user "
                           "specified times, no weakening is applied, and values provided in Viscosity prefactors "
                           "are applied after the specified time. A value can be specified for the background material "
                           "and compositional fields, for a total of N+1 where N is the number of "
                           "all compositional fields or only those corresponding to chemical compositions. Units: none.");
      }



      template <int dim>
      void
      ViscosityPrefactors<dim>::parse_parameters (ParameterHandler &prm)
      {
        // Retrieve the list of composition names
        std::vector<std::string> compositional_field_names = this->introspection().get_composition_names();

        // Retrieve the list of names of fields that represent chemical compositions, and not, e.g.,
        // plastic strain
        std::vector<std::string> chemical_field_names = this->introspection().chemical_composition_field_names();

        // Establish that a background field is required here
        compositional_field_names.insert(compositional_field_names.begin(), "background");
        chemical_field_names.insert(chemical_field_names.begin(),"background");


        Utilities::MapParsing::Options options(chemical_field_names, "Viscosity prefactors");
        options.list_of_allowed_keys = compositional_field_names;

        variable_viscosity_prefactors = Utilities::MapParsing::parse_map_to_double_array (prm.get("Viscosity prefactors"),
                                        options);
        weakening_times = Utilities::MapParsing::parse_map_to_double_array (prm.get("Weakening times"), options);
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
    template class ViscosityPrefactors<dim>; \
  }

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
