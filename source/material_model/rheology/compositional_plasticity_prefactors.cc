/*
  Copyright (C) 2024 by the authors of the ASPECT code.

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


#include <aspect/material_model/rheology/compositional_plasticity_prefactors.h>
#include <aspect/utilities.h>
#include <aspect/global.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/adiabatic_conditions/interface.h>


#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/parameter_handler.h>
#include <aspect/simulator_signals.h>
#include <aspect/parameters.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace Rheology
    {
      template <int dim>
      CompositionalPlasticityPrefactors<dim>::CompositionalPlasticityPrefactors ()
        = default;

      template <int dim>
      std::pair<double, double>
      CompositionalPlasticityPrefactors<dim>::compute_weakening_factor (const MaterialModel::MaterialModelInputs<dim> &in,
                                                                        const unsigned int volume_fraction_index,
                                                                        const unsigned int i,
                                                                        const double input_cohesion,
                                                                        const double input_friction_angle,
                                                                        const Point<dim> &position) const
      {
        std::pair<double, double> plasticity_parameters (input_cohesion, input_friction_angle);
        switch (prefactor_mechanism)
          {
            case none:
            {
              return plasticity_parameters;
            }
            case porosity:
            {
              // Use the porosity to compute the cohesion and friction angle
              // The porosity is computed in the MaterialModel::MaterialModelInputs
              const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
              plasticity_parameters.first *= std::max((max_porosities_for_cohesion_prefactors[volume_fraction_index] - in.composition[i][porosity_idx]) / max_porosities_for_cohesion_prefactors[volume_fraction_index], 0.0);
              plasticity_parameters.second *= std::max((max_porosities_for_friction_angle_prefactors[volume_fraction_index] - in.composition[i][porosity_idx]) / max_porosities_for_friction_angle_prefactors[volume_fraction_index], 0.0);

              plasticity_parameters.first = std::max(plasticity_parameters.first, minimum_cohesions[volume_fraction_index]);
              plasticity_parameters.second = std::max(plasticity_parameters.second, minimum_friction_angles[volume_fraction_index] * constants::degree_to_radians);
              return plasticity_parameters;
            }
            case function:
            {
              // Use a given function input per composition to get the prefactors
              Utilities::NaturalCoordinate<dim> point =
                this->get_geometry_model().cartesian_to_other_coordinates(position, coordinate_system_prefactor_function);

              // we get time passed as seconds (always) but may want
              // to reinterpret it in years
              if (this->convert_output_to_years())
                prefactor_function->set_time (this->get_time() / year_in_seconds);
              else
                prefactor_function->set_time (this->get_time());

              // determine the prefactors based on position and composition
              // The volume_fraction_index is based on the number of chemical compositional fields.
              // However, this plugin reads a function for every compositional field, regardless of
              // its type. Therefore we have to get the correct index.
              // If no fields or no chemical fields are present, but only background material, the index is zero.
              // If chemical fields are present, volume_fractions will be of size 1+n_chemical_composition_fields.
              // The size of chemical_composition_field_indices will be one less.
              unsigned int index = 0;
              if (this->introspection().composition_type_exists(CompositionalFieldDescription::chemical_composition))
                index = this->introspection().chemical_composition_field_indices()[volume_fraction_index-1];
              double prefactor_from_function =
                prefactor_function->value(Utilities::convert_array_to_point<dim>(point.get_coordinates()),index);

              // Multiply the cohesion and the friction angle by the prefactor.
              plasticity_parameters.first *= prefactor_from_function;
              plasticity_parameters.second *= prefactor_from_function * constants::degree_to_radians;

              return plasticity_parameters;
            }
          }
      }


      template <int dim>
      void
      CompositionalPlasticityPrefactors<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Maximum porosity for cohesion prefactor", "1.0",
                           Patterns::List(Patterns::Double(std::numeric_limits<double>::epsilon(), 1.0)),
                           "List of the maximum porosity value for calculating the prefactor for the cohesion. "
                           "entered as a volume fraction. If chosen to be 0.1 (10 percent porosity), the "
                           "cohesion will be linearly reduced from its default value when the model porosity is 0"
                           "to a maximum reduction when the model porosity is 0.1. Units: none.");

        prm.declare_entry ("Maximum porosity for friction angle prefactor", "1.0",
                           Patterns::List(Patterns::Double(std::numeric_limits<double>::epsilon(), 1.0)),
                           "List of the maximum porosity value for calculating the prefactor for the friction "
                           "angle, entered as a volume fraction. If chosen to be 0.1 (10 percent porosity), the "
                           "friction angle will be linearly reduced from its default value when the model porosity is 0"
                           "to a maximum reduction when the model porosity is 0.1. Units: none.");

        prm.declare_entry ("Plasticity prefactor scheme", "none",
                           Patterns::Selection("none|function|porosity"),
                           "Select what type of plasticity multiplicative prefactor scheme to apply. "
                           "Allowed entries are 'none', 'function', and 'porosity'. 'function' allows "
                           "the cohesion and friction angle prefactors to be specified using a function. "
                           "'porosity' determines a prefactor based on the modeled porosity, this scheme "
                           "requires that the user defines a compositional field called 'porosity'. 'none' "
                           "does not modify the cohesion or the friction angle. Units: none.");

        prm.declare_entry ("Minimum cohesions", "1e6",
                           Patterns::List(Patterns::Double(0.)),
                           "List of the minimum cohesions for limiting the reduction to the cohesion. Units: MPa.");

        prm.declare_entry ("Minimum friction angles", "1",
                           Patterns::List(Patterns::Double(0.)),
                           "List of the minimum friction angles for limiting the reduction to the friction angle. Units: none.");

      }



      template <int dim>
      void
      CompositionalPlasticityPrefactors<dim>::parse_parameters (ParameterHandler &prm)
      {

        const unsigned int n_fields = this->n_compositional_fields() + 1;
        if (prm.get ("Plasticity prefactor scheme") == "none")
          prefactor_mechanism = none;
        else if (prm.get ("Plasticity prefactor scheme") == "porosity")
          {
            AssertThrow(this->introspection().compositional_name_exists("porosity"),
                        ExcMessage("The 'porosity' Plasticity prefactor scheme work only if "
                                   "is a compositional field called porosity."));
            prefactor_mechanism = porosity;

            // Retrieve the list of chemical names
            std::vector<std::string> chemical_field_names = this->introspection().chemical_composition_field_names();
            std::vector<std::string> compositional_field_names = this->introspection().get_composition_names();

            // Establish that a background field is required here
            compositional_field_names.insert(compositional_field_names.begin(), "background");
            chemical_field_names.insert(chemical_field_names.begin(),"background");

            Utilities::MapParsing::Options options(chemical_field_names, "Maximum porosity for cohesion prefactor");

            options.list_of_allowed_keys = compositional_field_names;
            max_porosities_for_cohesion_prefactors = Utilities::MapParsing::parse_map_to_double_array (prm.get("Maximum porosity for cohesion prefactor"),
                                                     options);
            max_porosities_for_friction_angle_prefactors = Utilities::MapParsing::parse_map_to_double_array (prm.get("Maximum porosity for friction angle prefactor"),
                                                           options);

            minimum_cohesions = Utilities::MapParsing::parse_map_to_double_array (prm.get("Minimum cohesions"),
                                                                                  options);
            minimum_friction_angles = Utilities::MapParsing::parse_map_to_double_array (prm.get("Minimum friction angles"),
                                                                                        options);
          }
        else if (prm.get ("Plasticity prefactor scheme") == "function")
          {
            prefactor_mechanism = function;
            prm.enter_subsection("Prefactor function");
            {
              coordinate_system_prefactor_function = Utilities::Coordinates::string_to_coordinate_system(prm.get("Coordinate system"));
              try
                {
                  prefactor_function
                    = std::make_unique<Functions::ParsedFunction<dim>>(n_fields);
                  prefactor_function->parse_parameters (prm);
                }
              catch (...)
                {
                  std::cerr << "FunctionParser failed to parse\n"
                            << "\t friction function\n"
                            << "with expression \n"
                            << "\t' " << prm.get("Function expression") << "'";
                  throw;
                }
            }
            prm.leave_subsection();
          }
        else
          AssertThrow(false, ExcMessage("Not a valid plasticity prefactor scheme"));
      }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    namespace Rheology
    {
#define INSTANTIATE(dim) \
  template class CompositionalPlasticityPrefactors<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
    }
  }
}
