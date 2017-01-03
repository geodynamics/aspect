/*
  Copyright (C) 2016 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/material_model/ascii_reference_profile.h>

using namespace dealii;

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    void
    AsciiReferenceProfile<dim>::initialize ()
    {
      profile.initialize(7,this->get_mpi_communicator());
    }

    template <int dim>
    void
    AsciiReferenceProfile<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      for (unsigned int i=0; i < in.temperature.size(); ++i)
        {
          const Point<dim> position = in.position[i];
          const double temperature_deviation = in.temperature[i] - this->get_adiabatic_conditions().temperature(position);
          const double pressure_deviation = in.pressure[i] - this->get_adiabatic_conditions().pressure(position);

          const double depth = this->get_geometry_model().depth(position);
          const Point<1> profile_position(depth);

          double visc_temperature_dependence = std::max(std::min(std::exp(-thermal_viscosity_exponent*temperature_deviation/this->get_adiabatic_conditions().temperature(position)),1e3),1e-3);
          if (std::isnan(visc_temperature_dependence))
            visc_temperature_dependence = 1.0;

          double visc_depth_dependence = viscosity_prefactors[0];
          for (unsigned int j=0; j < transition_depths.size(); ++j)
            {
              if(depth>transition_depths[j])
                visc_depth_dependence = viscosity_prefactors[j+1];
            }

          out.viscosities[i] = viscosity * visc_temperature_dependence * visc_depth_dependence;

          out.thermal_conductivities[i] = thermal_conductivity;

          out.thermal_expansion_coefficients[i] = profile.get_data_component(profile_position,4);
          out.specific_heat[i] = profile.get_data_component(profile_position,5);
          out.compressibilities[i] = profile.get_data_component(profile_position,6);

          out.densities[i] = profile.get_data_component(profile_position,2)
                             * (1.0 - out.thermal_expansion_coefficients[i] * temperature_deviation)
                             * (1.0 + out.compressibilities[i] * pressure_deviation);

          out.entropy_derivative_pressure[i] = 0.0;
          out.entropy_derivative_temperature[i] = 0.0;

          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            out.reaction_terms[i][c] = 0.0;
        }
    }


    template <int dim>
    double
    AsciiReferenceProfile<dim>::
    reference_viscosity () const
    {
      return viscosity;
    }

    template <int dim>
    double
    AsciiReferenceProfile<dim>::
    reference_density () const
    {
      return 3340;
    }

    template <int dim>
    double
    AsciiReferenceProfile<dim>::
    reference_cp () const
    {
      return 1200;
    }

    template <int dim>
    bool
    AsciiReferenceProfile<dim>::
    is_compressible () const
    {
      return true;
    }



    template <int dim>
    void
    AsciiReferenceProfile<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Ascii reference profile");
        {
          prm.declare_entry ("Thermal conductivity", "4.0",
                             Patterns::Double (0),
                             "Reference conductivity");
          prm.declare_entry ("Viscosity", "1e21",
                             Patterns::Double (0),
                             "Viscosity");
          prm.declare_entry ("Use TALA", "false",
                             Patterns::Bool (),
                             "Whether to use the TALA instead of the ALA "
                             "approximation.");
          prm.declare_entry ("Thermal viscosity exponent", "0.0",
                             Patterns::Double (0),
                             "The temperature dependence of viscosity. Dimensionless exponent.");
          prm.declare_entry ("Transition depths", "1.5e5, 4.1e5, 6.6e5",
                             Patterns::List (Patterns::Double(0)),
                             "A list of depths where the viscosity changes. Values must "
                             "monotonically increase. "
                             "Units: $m$.");
          prm.declare_entry ("Viscosity prefactors", "10, 0.1, 1, 10",
                             Patterns::List (Patterns::Double(0)),
                             "A list of prefactors for the viscosity for each phase. The reference "
                             "viscosity will be multiplied by this factor to get the corresponding "
                             "viscosity for each phase. "
                             "List must have one more entry than Phase transition depths. "
                             "Units: non-dimensional.");
        }
        prm.leave_subsection();

        aspect::Utilities::AsciiDataProfile<dim>::declare_parameters(prm,
            "$ASPECT_SOURCE_DIR/data/adiabatic-conditions/ascii-data/",
            "simple_test.txt");
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    AsciiReferenceProfile<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Ascii reference profile");
        {
          tala                 = prm.get_bool ("Use TALA");
          thermal_conductivity = prm.get_double ("Thermal conductivity");
          viscosity            = prm.get_double ("Viscosity");
          thermal_viscosity_exponent = prm.get_double ("Thermal viscosity exponent");
          transition_depths    = Utilities::string_to_double
                                 (Utilities::split_string_list(prm.get ("Transition depths")));
          viscosity_prefactors = Utilities::string_to_double
                                 (Utilities::split_string_list(prm.get ("Viscosity prefactors")));

          // make sure to check against the depth lists for size errors, since using depth
          if (viscosity_prefactors.size() != transition_depths.size()+1)
            AssertThrow(false, ExcMessage("Error: The list of Viscosity prefactors needs to have exactly "
                                          "one more entry than the list of Transition depths. "));
        }
        prm.leave_subsection();

        profile.parse_parameters(prm);
      }
      prm.leave_subsection();

      // Declare dependencies on solution variables
      this->model_dependence.viscosity = NonlinearDependence::none;
      this->model_dependence.density = NonlinearDependence::pressure | NonlinearDependence::temperature;
      this->model_dependence.compressibility = NonlinearDependence::none;
      this->model_dependence.specific_heat = NonlinearDependence::none;
      this->model_dependence.thermal_conductivity = NonlinearDependence::none;

    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(AsciiReferenceProfile,
                                   "ascii reference profile",
                                   "A material model that reads in a reference "
                                   "state from a text file.")
  }
}
