/*
  Copyright (C) 2019 - 2024 by the authors of the ASPECT code.

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

#include <aspect/material_model/rheology/friction_models.h>

#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/parameter_handler.h>
#include <aspect/utilities.h>
#include <aspect/geometry_model/interface.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace Rheology
    {
      template <int dim>
      double
      FrictionModels<dim>::
      compute_friction_angle(const double current_edot_ii,
                             const unsigned int volume_fraction_index,
                             const double static_friction_angle,
                             const Point<dim> &position) const
      {

        switch (friction_mechanism)
          {
            case static_friction:
            {
              return static_friction_angle;
            }
            case dynamic_friction:
            {
              // Calculate effective steady-state friction coefficient.
              // This is based on the former material model dynamic friction.
              // The formula below is equivalent to the equation 13 in \cite{van_dinther_seismic_2013}.
              // Although here the dynamic friction coefficient is directly specified. In addition,
              // we also use a reference strain rate in place of a characteristic
              // velocity divided by local element size. This reference strain rate is called
              // the dynamic characteristic strain rate and is used to compute what value between
              // dynamic and static angle of internal friction should be used.
              // Furthermore a smoothness exponent X is added, which determines whether the
              // friction vs strain rate curve is rather step-like or more gradual.
              // mu  = mu_d + (mu_s - mu_d) / ( (1 + strain_rate_dev_inv2/dynamic_characteristic_strain_rate)^X );
              // Angles of friction are used in radians within ASPECT. The coefficient
              // of friction (mu) is the tangent of the internal angle of friction, hence convergence is needed.
              // The incoming variable static_friction_angle holds the static friction angle. It's value is
              // then updated with the strain rate dependent friction angle which is returned by the function.
              const double mu = (std::tan(dynamic_angles_of_internal_friction[volume_fraction_index])
                                 + (std::tan(static_friction_angle) - std::tan(dynamic_angles_of_internal_friction[volume_fraction_index]))
                                 / (1. + std::pow((current_edot_ii / dynamic_characteristic_strain_rate),
                                                  dynamic_friction_smoothness_exponent)));
              const double dynamic_friction_angle = std::atan (mu);
              Assert((mu < 1) && (0 < dynamic_friction_angle) && (dynamic_friction_angle <= 1.6), ExcMessage(
                       "The friction coefficient should be larger than zero and smaller than 1. "
                       "The friction angle should be smaller than 1.6 rad."));
              return dynamic_friction_angle;
            }
            case function:
            {
              // Use a given function input per composition to get the friction angle
              Utilities::NaturalCoordinate<dim> point =
                this->get_geometry_model().cartesian_to_other_coordinates(position, coordinate_system_friction_function);

              // we get time passed as seconds (always) but may want
              // to reinterpret it in years
              if (this->convert_output_to_years())
                friction_function->set_time (this->get_time() / year_in_seconds);
              else
                friction_function->set_time (this->get_time());

              // determine the friction angle based on position and composition
              // The volume_fraction_index is based on the number of chemical compositional fields.
              // However, this plugin reads a function for every compositional field, regardless of
              // its type. Therefore we have to get the correct index.
              // If no fields or no chemical fields are present, but only background material, the index is zero.
              // If chemical fields are present, volume_fractions will be of size 1+n_chemical_composition_fields.
              // The size of chemical_composition_field_indices will be one less.
              unsigned int index = 0;
              if (this->introspection().composition_type_exists(CompositionalFieldDescription::chemical_composition))
                index = this->introspection().chemical_composition_field_indices()[volume_fraction_index-1];
              double friction_from_function =
                friction_function->value(Utilities::convert_array_to_point<dim>(point.get_coordinates()),index);

              // Convert angles from degrees to radians
              friction_from_function *= constants::degree_to_radians;

              return friction_from_function;
            }
          }
        // we should never get here, return something anyway, so the compiler does not complain...
        AssertThrow (false, ExcMessage("Unknown friction model."));
        return static_friction_angle;
      }



      template <int dim>
      FrictionMechanism
      FrictionModels<dim>::
      get_friction_mechanism() const
      {
        return friction_mechanism;
      }



      template <int dim>
      void
      FrictionModels<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Friction mechanism", "none",
                           Patterns::Selection("none|dynamic friction|function"),
                           "Whether to make the friction angle dependent on strain rate or not. This rheology "
                           "is intended to be used together with the visco-plastic rheology model."
                           "\n\n"
                           "\\item ``none'': No dependence of the friction angle is applied. "
                           "\n\n"
                           "\\item ``dynamic friction'': The friction angle is rate dependent."
                           "When 'dynamic angles of internal friction' are specified, "
                           "the friction angle will be weakened for high strain rates with: "
                           "$\\mu = \\mu_d + \\frac{\\mu_s-\\mu_d}{1+\\frac{\\dot{\\epsilon}_{ii}}{\\dot{\\epsilon}_C}}^x$  "
                           "where $\\mu_s$ and $\\mu_d$ are the friction angles at low and high strain rates, "
                           "respectively. $\\dot{\\epsilon}_{ii}$ is the second invariant of the strain rate and "
                           "$\\dot{\\epsilon}_C$ is the 'dynamic characteristic strain rate' where $\\mu = (\\mu_s+\\mu_d)/2$. "
                           "The 'dynamic friction smoothness exponent' x controls how "
                           "smooth or step-like the change from $\\mu_s$ to $\\mu_d$ is. "
                           "The equation is modified after Equation (13) in \\cite{van_dinther_seismic_2013}. "
                           "$\\mu_s$ and $\\mu_d$ can be specified by setting 'Angles of internal friction' and "
                           "'Dynamic angles of internal friction', respectively. "
                           "This relationship is similar to rate-and-state friction constitutive relationships, which "
                           "are applicable to the strength of rocks during earthquakes."
                           "\n\n"
                           "\\item ``function'': Specify the friction angle as a function of space and time "
                           "for each compositional field.");

        // Dynamic friction parameters
        prm.declare_entry ("Dynamic characteristic strain rate", "1e-12",
                           Patterns::Double (0),
                           "The characteristic strain rate value at which the angle of friction is "
                           "equal to $\\mu = (\\mu_s+\\mu_d)/2$. When the effective strain rate "
                           "is very high, the dynamic angle of friction is taken, when it is very low, "
                           "the static angle of internal friction is used. Around the dynamic characteristic "
                           "strain rate, there is a smooth gradient from the static to the dynamic angle "
                           "of internal friction. "
                           "Units: \\si{\\per\\second}.");

        prm.declare_entry ("Dynamic angles of internal friction", "2",
                           Patterns::List(Patterns::Double(0)),
                           "List of dynamic angles of internal friction, $\\phi$, for background material and compositional "
                           "fields, for a total of N$+$1 values, where N is the number of all compositional fields or only "
                           "those corresponding to chemical compositions. "
                           "Dynamic angles of friction are used as the current friction angle when the effective "
                           "strain rate is well above the 'dynamic characteristic strain rate'. "
                           "Units: \\si{\\degree}.");

        prm.declare_entry ("Dynamic friction smoothness exponent", "1",
                           Patterns::Double (0),
                           "An exponential factor in the equation for the calculation of the friction angle when a "
                           "static and a dynamic angle of internal friction are specified. A factor of 1 returns the equation "
                           "to Equation (13) in \\cite{van_dinther_seismic_2013}. A factor between 0 and 1 makes the "
                           "curve of the friction angle vs. the strain rate smoother, while a factor $>$ 1 makes "
                           "the change between static and dynamic friction angle more steplike. "
                           "Units: none.");

        /**
         * If friction is specified as a function input.
         */
        prm.enter_subsection("Friction function");
        {
          /**
           * The function to specify the friction angle per composition can be declared in dependence
           * of depth, cartesian coordinates or spherical coordinates. Note that the order
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

          Functions::ParsedFunction<dim>::declare_parameters(prm,1);
        }
        prm.leave_subsection();
      }



      template <int dim>
      void
      FrictionModels<dim>::parse_parameters (ParameterHandler &prm)
      {
        // Friction dependence parameters
        if (prm.get ("Friction mechanism") == "none")
          friction_mechanism = static_friction;
        else if (prm.get ("Friction mechanism") == "dynamic friction")
          friction_mechanism = dynamic_friction;
        else if (prm.get ("Friction mechanism") == "function")
          friction_mechanism = function;
        else
          AssertThrow(false, ExcMessage("Not a valid friction mechanism option!"));

        // Dynamic friction parameters
        dynamic_characteristic_strain_rate = prm.get_double("Dynamic characteristic strain rate");

        // Retrieve the list of composition names
        std::vector<std::string> compositional_field_names = this->introspection().get_composition_names();

        // Retrieve the list of names of fields that represent chemical compositions, and not, e.g.,
        // plastic strain
        std::vector<std::string> chemical_field_names = this->introspection().chemical_composition_field_names();

        // Establish that a background field is required here
        compositional_field_names.insert(compositional_field_names.begin(), "background");
        chemical_field_names.insert(chemical_field_names.begin(), "background");

        Utilities::MapParsing::Options options(chemical_field_names, "Dynamic angles of internal friction");
        options.list_of_allowed_keys = compositional_field_names;

        dynamic_angles_of_internal_friction = Utilities::MapParsing::parse_map_to_double_array (prm.get("Dynamic angles of internal friction"),
                                              options);

        // Convert angles from degrees to radians
        for (double &angle : dynamic_angles_of_internal_friction)
          {
            AssertThrow(angle <= 90,
                        ExcMessage("Dynamic angles of friction must be <= 90 degrees"));
            angle *= constants::degree_to_radians;
          }

        dynamic_friction_smoothness_exponent = prm.get_double("Dynamic friction smoothness exponent");


        // Get the number of fields for composition-dependent material properties
        // including the background field.
        // TODO Make sure functions only have to be specified per chemical composition,
        // but can still be specified for all fields for backwards compatibility.
        const unsigned int n_fields = this->n_compositional_fields() + 1;

        // if friction is specified as a function
        if (friction_mechanism == function)
          {
            prm.enter_subsection("Friction function");
            {
              coordinate_system_friction_function = Utilities::Coordinates::string_to_coordinate_system(prm.get("Coordinate system"));
              try
                {
                  friction_function
                    = std::make_unique<Functions::ParsedFunction<dim>>(n_fields);
                  friction_function->parse_parameters (prm);
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

      }
    }
  }
}



// Explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
#define INSTANTIATE(dim) \
  namespace Rheology \
  { \
    template class FrictionModels<dim>; \
  }

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
