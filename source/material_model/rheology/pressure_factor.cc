/*
  Copyright (C) 2019 - 2026 by the authors of the ASPECT code.

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

#include <aspect/material_model/rheology/pressure_factor.h>
#include <aspect/utilities.h>
#include <aspect/global.h>
#include <aspect/geometry_model/interface.h>

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
      double
      PressureFactor<dim>::
      compute_pressure_factor(const MaterialModel::MaterialModelInputs<dim> &in,
                              const unsigned int q,
                              const unsigned int volume_fraction_index) const
      {
        double current_pressure_factor = 0.;

        switch (pressure_factor_model)
          {
            case PressureFactorModel::none:
            {
              break;
            }
            case PressureFactorModel::function:
            {
              // In this model, there will be a function that defines the pressure
              // factor for each composition. We need to get this and average them.
              double compositional_pressure_factor = pressure_factor_function(in.position[q],volume_fraction_index);

              current_pressure_factor = compositional_pressure_factor;

              break;
            }
            default:
            {
              AssertThrow(false, ExcNotImplemented());
              break;
            }
          }

        return current_pressure_factor;
      }

      template <int dim>
      double
      PressureFactor<dim>::
      compute_scaled_pressure(const MaterialModel::MaterialModelInputs<dim> &in,
                              const double &pressure_to_scale,
                              const unsigned int q,
                              const unsigned int volume_fraction_index) const
      {
        // First find the pressure scaling factor.
        const double lambda = compute_pressure_factor(in, q, volume_fraction_index);

        // Return either the solid pressure scaled by lamdba, as would be expected
        // by something like pore fluid pressure, or return the fluid pressure.
        const double scaled_pressure = return_solid_scaling
                                       ? pressure_to_scale * (1 - lambda)
                                       : pressure_to_scale * lambda;

        return scaled_pressure;
      }

      template <int dim>
      double
      PressureFactor<dim>::
      pressure_factor_function(const Point<dim> &position, const unsigned int n_comp) const
      {
        const Utilities::NaturalCoordinate<dim> point =
          this->get_geometry_model().cartesian_to_other_coordinates(position, coordinate_system);

        // Convert time for the function based on prm file.
        if (this->convert_output_to_years())
          function->set_time (this->get_time() / year_in_seconds);
        else
          function->set_time (this->get_time());

        // TODO: At some point it needs to be checked how this works with a deforming surface mesh
        // if we use a depth-dependent function.
        return function->value(Utilities::convert_array_to_point<dim>(point.get_coordinates()),n_comp);
      }

      template <int dim>
      PressureFactorModel
      PressureFactor<dim>::
      get_pressure_factor_model() const
      {
        return pressure_factor_model;
      }

      template <int dim>
      void
      PressureFactor<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Pressure factor model");
        {
          prm.declare_entry ("Pressure factor model", "none",
                             Patterns::Selection("none|function"),
                             "Model for pressure factor. At the moment, there is only a"
                             "function implementation.");
          prm.declare_entry ("Use solid pressure scaling", "true",
                             Patterns::Bool (),
                             "Whether the user wants to scale the solid pressure as "
                             "P * (1 - lambda), or return the fluid pressure as "
                             "P * lambda.");

          prm.enter_subsection("Function");
          {
            /**
             * Choose the coordinates to evaluate the pressure factor
             * function. The function can be declared in dependence of depth,
             * cartesian coordinates or spherical coordinates. Note that the order
             * of spherical coordinates is r,phi,theta and not r,theta,phi, since
             * this allows for dimension independent expressions.
             */
            prm.declare_entry ("Coordinate system", "cartesian",
                               Patterns::Selection ("cartesian|spherical|depth"),
                               "A selection that determines the assumed coordinate "
                               "system for the function variables. Allowed values "
                               "are `cartesian', `spherical', and `depth'. `spherical' coordinates "
                               "are interpreted as r,phi or r,phi,theta in 2d/3d "
                               "respectively with theta being the polar angle. `depth' "
                               "will create a function, in which only the first "
                               "parameter is non-zero, which is interpreted to "
                               "be the depth of the point.");

            Functions::ParsedFunction<dim>::declare_parameters (prm, 1);
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }

      template <int dim>
      void
      PressureFactor<dim>::parse_parameters (ParameterHandler &prm,
                                             const std::unique_ptr<std::vector<unsigned int>> &expected_n_phases_per_composition)
      {
        prm.enter_subsection("Pressure factor model");
        {
          // Find which pressure factor model has been selected.
          if (prm.get ("Pressure factor model") == "none")
            pressure_factor_model = PressureFactorModel::none;
          else if (prm.get ("Pressure factor model") == "function")
            pressure_factor_model = PressureFactorModel::function;
          else
            AssertThrow(false, ExcMessage("Not a valid pressure factor model!"));

          return_solid_scaling = prm.get_bool ("Use solid pressure scaling");

          // Retrieve the list of composition names
          std::vector<std::string> compositional_field_names = this->introspection().get_composition_names();

          // Retrieve the list of names of fields that represent chemical compositions, and not, e.g.,
          // plastic strain
          std::vector<std::string> chemical_field_names = this->introspection().chemical_composition_field_names();

          // Establish that a background field is required here
          compositional_field_names.insert(compositional_field_names.begin(), "background");
          chemical_field_names.insert(chemical_field_names.begin(), "background");

          // We use Angles of internal friction here because it should have the same keys.
          Utilities::MapParsing::Options options(chemical_field_names, "Angles of internal friction");
          options.list_of_allowed_keys = compositional_field_names;

          if (expected_n_phases_per_composition)
            {
              options.allow_multiple_values_per_key = true;
              options.n_values_per_key = *expected_n_phases_per_composition;

              // check_values_per_key is required to be true to duplicate single values
              // if they are to be used for all phases associated with a given key.
              options.check_values_per_key = true;
            }

          if (pressure_factor_model == PressureFactorModel::function)
            {
              prm.enter_subsection("Function");
              {
                coordinate_system = Utilities::Coordinates::string_to_coordinate_system(prm.get("Coordinate system"));
              }

              try
                {
                  function
                    = std::make_unique<Functions::ParsedFunction<dim>>(chemical_field_names.size());
                  function->parse_parameters (prm);
                }
              catch (...)
                {
                  std::cerr << "ERROR: FunctionParser failed to parse\n"
                            << "\t'Pressure factor function.Function'\n"
                            << "with expression\n"
                            << "\t'" << prm.get("Function expression") << "'\n"
                            << "More information about the cause of the parse error \n"
                            << "is shown below.\n";
                  throw;
                }
              prm.leave_subsection();
            }
        }
        prm.leave_subsection();

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
    template class PressureFactor<dim>; \
  }

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
