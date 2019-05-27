/*
  Copyright (C) 2011 - 2018 by the authors of the ASPECT code.

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


#include <aspect/material_model/multicomponent.h>
#include <aspect/simulator.h>
#include <aspect/utilities.h>

#include <numeric>


namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    void
    Multicomponent<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      for (unsigned int i=0; i < in.temperature.size(); ++i)
        {
          const double temperature = in.temperature[i];
          const std::vector<double> composition = in.composition[i];
          const std::vector<double> volume_fractions = MaterialUtilities::compute_volume_fractions(composition);

          out.viscosities[i] = MaterialUtilities::average_value (volume_fractions, viscosities, viscosity_averaging);
          out.specific_heat[i] = MaterialUtilities::average_value (volume_fractions, specific_heats, MaterialUtilities::arithmetic);


          // Arithmetic averaging of thermal conductivities
          // This may not be strictly the most reasonable thing, but for most Earth materials we hope
          // that they do not vary so much that it is a big problem.
          out.thermal_conductivities[i] = MaterialUtilities::average_value (volume_fractions, thermal_conductivities, MaterialUtilities::arithmetic);

          double density = 0.0;
          for (unsigned int j=0; j < volume_fractions.size(); ++j)
            {
              // not strictly correct if thermal expansivities are different, since we are interpreting
              // these compositions as volume fractions, but the error introduced should not be too bad.
              const double temperature_factor= (1.0 - thermal_expansivities[j] * (temperature - reference_T));
              density += volume_fractions[j] * densities[j] * temperature_factor;
            }
          out.densities[i] = density;


          out.thermal_expansion_coefficients[i] = MaterialUtilities::average_value (volume_fractions, thermal_expansivities, MaterialUtilities::arithmetic);


          // Compressibility at the given positions.
          // The compressibility is given as
          // $\frac 1\rho \frac{\partial\rho}{\partial p}$.
          // (here we use an incompressible medium)
          out.compressibilities[i] = 0.0;
          // Pressure derivative of entropy at the given positions.
          out.entropy_derivative_pressure[i] = 0.0;
          // Temperature derivative of entropy at the given positions.
          out.entropy_derivative_temperature[i] = 0.0;
          // Change in composition due to chemical reactions at the
          // given positions. The term reaction_terms[i][c] is the
          // change in compositional field c at point i.
          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            out.reaction_terms[i][c] = 0.0;

        }
    }

    template <int dim>
    double
    Multicomponent<dim>::
    reference_viscosity () const
    {
      return viscosities[0]; // background
    }

    template <int dim>
    bool
    Multicomponent<dim>::
    is_compressible () const
    {
      return false;
    }

    template <int dim>
    void
    Multicomponent<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Multicomponent");
        {
          prm.declare_entry ("Reference temperature", "293",
                             Patterns::Double (0),
                             "The reference temperature $T_0$. Units: $\\text{K}$.");
          prm.declare_entry ("Densities", "3300.",
                             Patterns::Anything(),
                             "List of densities for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one value is given, then all use the same value.  Units: $kg / m^3$");
          prm.declare_entry ("Viscosities", "1.e21",
                             Patterns::Anything(),
                             "List of viscosities for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one value is given, then all use the same value. Units: $Pa \\, s$");
          prm.declare_entry ("Thermal expansivities", "4.e-5",
                             Patterns::Anything(),
                             "List of thermal expansivities for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one value is given, then all use the same value. Units: $1/K$");
          prm.declare_entry ("Specific heats", "1250.",
                             Patterns::Anything(),
                             "List of specific heats $C_p$ for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one value is given, then all use the same value. Units: $J /kg /K$");
          prm.declare_entry ("Thermal conductivities", "4.7",
                             Patterns::Anything(),
                             "List of thermal conductivities for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one value is given, then all use the same value. Units: $W/m/K$.");
          prm.declare_entry("Viscosity averaging scheme", "harmonic",
                            Patterns::Selection("arithmetic|harmonic|geometric|maximum composition"),
                            "When more than one compositional field is present at a point "
                            "with different viscosities, we need to come up with an average "
                            "viscosity at that point.  Select a weighted harmonic, arithmetic, "
                            "geometric, or maximum composition.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    Multicomponent<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Multicomponent");
        {
          reference_T = prm.get_double ("Reference temperature");

          viscosity_averaging = MaterialUtilities::parse_compositional_averaging_operation ("Viscosity averaging scheme",
                                prm);

          // Establish that a background field is required here
          const bool has_background_field = true;

          // Retrieve the list of composition names
          const std::vector<std::string> list_of_composition_names = this->introspection().get_composition_names();

          // Parse multicomponent properties
          densities = Utilities::parse_map_to_double_array (prm.get("Densities"),
                                                            list_of_composition_names,
                                                            has_background_field,
                                                            "Densities");

          viscosities = Utilities::parse_map_to_double_array (prm.get("Viscosities"),
                                                              list_of_composition_names,
                                                              has_background_field,
                                                              "Viscosities");

          thermal_conductivities = Utilities::parse_map_to_double_array (prm.get("Thermal conductivities"),
                                                                         list_of_composition_names,
                                                                         has_background_field,
                                                                         "Thermal conductivities");

          thermal_expansivities = Utilities::parse_map_to_double_array (prm.get("Thermal expansivities"),
                                                                        list_of_composition_names,
                                                                        has_background_field,
                                                                        "Thermal expansivities");

          specific_heats = Utilities::parse_map_to_double_array (prm.get("Specific heats"),
                                                                 list_of_composition_names,
                                                                 has_background_field,
                                                                 "Specific heats");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Declare dependencies on solution variables
      this->model_dependence.viscosity = NonlinearDependence::compositional_fields;
      this->model_dependence.density = NonlinearDependence::temperature | NonlinearDependence::compositional_fields;
      this->model_dependence.compressibility = NonlinearDependence::none;
      this->model_dependence.specific_heat = NonlinearDependence::compositional_fields;
      this->model_dependence.thermal_conductivity = NonlinearDependence::compositional_fields;
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(Multicomponent,
                                   "multicomponent",
                                   "This model is for use with an arbitrary number of compositional fields, where each field"
                                   " represents a rock type which can have completely different properties from the others."
                                   " However, each rock type itself has constant material properties.  The value of the "
                                   " compositional field is interpreted as a volume fraction. If the sum of the fields is"
                                   " greater than one, they are renormalized.  If it is less than one, material properties "
                                   " for ``background mantle'' make up the rest. When more than one field is present, the"
                                   " material properties are averaged arithmetically.  An exception is the viscosity,"
                                   " where the averaging should make more of a difference.  For this, the user selects"
                                   " between arithmetic, harmonic, geometric, or maximum composition averaging.")
  }
}
