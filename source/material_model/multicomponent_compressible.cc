/*
  Copyright (C) 2011 - 2017 by the authors of the ASPECT code.

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


#include <aspect/material_model/multicomponent_compressible.h>
#include <aspect/simulator.h>
#include <aspect/utilities.h>

#include <numeric>


namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    double
    MulticomponentCompressible<dim>::
    average_value ( const std::vector<double> &volume_fractions,
                    const std::vector<double> &parameter_values,
                    const enum AveragingScheme &average_type) const
    {
      double averaged_parameter = 0.0;

      switch (average_type)
        {
          case arithmetic:
          {
            for (unsigned int i=0; i< volume_fractions.size(); ++i)
              averaged_parameter += volume_fractions[i]*parameter_values[i];
            break;
          }
          case harmonic:
          {
            for (unsigned int i=0; i< volume_fractions.size(); ++i)
              averaged_parameter += volume_fractions[i]/(parameter_values[i]);
            averaged_parameter = 1.0/averaged_parameter;
            break;
          }
          case geometric:
          {
            for (unsigned int i=0; i < volume_fractions.size(); ++i)
              averaged_parameter += volume_fractions[i]*std::log(parameter_values[i]);
            averaged_parameter = std::exp(averaged_parameter);
            break;
          }
          case maximum_composition:
          {
            const unsigned int i = (unsigned int)(std::max_element( volume_fractions.begin(),
                                                                    volume_fractions.end() )
                                                  - volume_fractions.begin());
            averaged_parameter = parameter_values[i];
            break;
          }
          default:
          {
            AssertThrow( false, ExcNotImplemented() );
            break;
          }
        }
      return averaged_parameter;
    }


    template <int dim>
    void
    MulticomponentCompressible<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      for (unsigned int i=0; i < in.temperature.size(); ++i)
        {
          const double pressure = in.pressure[i];
          const double temperature = in.temperature[i];
          const std::vector<double> composition = in.composition[i];
          const std::vector<double> reference_volume_fractions = compute_volume_fractions(composition);

          // Get the mass fractions from the declared reference volume fractions and densities
          unsigned int n_fields = reference_volume_fractions.size();
          std::vector<double> mass_fractions(n_fields);
          std::vector<double> volume_fractions(n_fields);

          double sum_mass_fractions = 0.0;
          for (unsigned int j=0; j < n_fields; ++j)
            {
              mass_fractions[j] = reference_volume_fractions[j] * reference_densities[j];
              sum_mass_fractions += mass_fractions[j];
            }

          // Averaging of the material properties
          // These are correct when the compositional fields correspond to distinct phases, and all are at the same temperature and pressure
          double invdensity = 0.;
          double Cp = 0.;
          double aoverrho = 0.;
          double boverrho = 0.;

          for (unsigned int j=0; j < n_fields; ++j)
            {
              mass_fractions[j] /= sum_mass_fractions;
              const double ak = reference_thermal_expansivities[j]/reference_isothermal_compressibilities[j];
              const double f = (1. + (pressure - ak*(temperature - reference_temperatures[j])) \
                                * reference_Kprimes[j]*reference_isothermal_compressibilities[j]);
              const double density = reference_densities[j]*std::pow(f, 1./reference_Kprimes[j]);
              const double isothermal_compressibility = reference_isothermal_compressibilities[j]/f;
              const double thermal_expansivity = reference_thermal_expansivities[j]/f;
              const double specific_heat = reference_specific_heats[j] - (temperature*reference_thermal_expansivities[j] \
                                                                          * ak * std::pow(f, -1.-(1./reference_Kprimes[j])) \
                                                                          / reference_densities[j]);

              volume_fractions[j] = mass_fractions[j] / density;
              invdensity += mass_fractions[j] / density;
              Cp += mass_fractions[j] * specific_heat;
              aoverrho += mass_fractions[j] * thermal_expansivity / density;
              boverrho += mass_fractions[j] * isothermal_compressibility / density;
            }

          for (unsigned int j=0; j < n_fields; ++j)
            {
              volume_fractions[j] /= invdensity;
            }

          out.densities[i] = 1./invdensity;
          out.specific_heat[i] = Cp;
          out.thermal_expansion_coefficients[i] = aoverrho*out.densities[i];
          out.compressibilities[i] = boverrho*out.densities[i];

          // Arithmetic averaging of thermal conductivities based on volume fractions
          // This may not be strictly the most reasonable thing, but for most Earth materials we hope
          // that they do not vary so much that it is a big problem.
          out.thermal_conductivities[i] = average_value ( volume_fractions, thermal_conductivities, arithmetic);

          // User-defined averaging of the viscosities
          out.viscosities[i] = average_value ( volume_fractions, viscosities, viscosity_averaging);

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
    MulticomponentCompressible<dim>::
    reference_viscosity () const
    {
      return viscosities[0]; // background
    }

    template <int dim>
    bool
    MulticomponentCompressible<dim>::
    is_compressible () const
    {
      return true;
    }

    template <int dim>
    void
    MulticomponentCompressible<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Multicomponent compressible");
        {
          prm.declare_entry ("Reference temperatures", "298.15",
                             Patterns::List(Patterns::Double(0)),
                             "List of reference temperatures $T_0$ for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one value is given, then all use the same value. Units: $K$.");
          prm.declare_entry ("Reference densities", "3300.",
                             Patterns::List(Patterns::Double(0)),
                             "List of densities for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one value is given, then all use the same value.  Units: $kg / m^3$");
          prm.declare_entry ("Reference isothermal compressibilities", "4e-12",
                             Patterns::List(Patterns::Double(0)),
                             "List of isothermal compressibilities for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one value is given, then all use the same value. "
                             "Units: $1/Pa$.");
          prm.declare_entry ("Reference Kprimes", "4.",
                             Patterns::List(Patterns::Double(0)),
                             "List of pressure derivatives of the bulk moduli for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one value is given, then all use the same value. "
                             "Units: [].");
          prm.declare_entry ("Reference thermal expansivities", "4.e-5",
                             Patterns::List(Patterns::Double(0)),
                             "List of thermal expansivities for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one value is given, then all use the same value. Units: $1/K$");
          prm.declare_entry ("Reference specific heats", "1250.",
                             Patterns::List(Patterns::Double(0)),
                             "List of specific heats $C_p$ for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one value is given, then all use the same value. Units: $J /kg /K$");
          prm.declare_entry ("Thermal conductivities", "4.7",
                             Patterns::List(Patterns::Double(0)),
                             "List of thermal conductivities for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one value is given, then all use the same value. Units: $W/m/K$.");
          prm.declare_entry ("Viscosities", "1.e21",
                             Patterns::List(Patterns::Double(0)),
                             "List of viscosities for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one value is given, then all use the same value. Units: $Pa \\, s$");
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
    MulticomponentCompressible<dim>::parse_parameters (ParameterHandler &prm)
    {
      // not pretty, but we need to get the number of compositional fields before
      // simulator access has been initialized here...
      unsigned int n_foreground_fields;
      prm.enter_subsection ("Compositional fields");
      {
        n_foreground_fields = prm.get_integer ("Number of fields");
      }
      prm.leave_subsection();

      const unsigned int n_fields= n_foreground_fields + 1;


      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Multicomponent compressible");
        {

          if (prm.get ("Viscosity averaging scheme") == "harmonic")
            viscosity_averaging = harmonic;
          else if (prm.get ("Viscosity averaging scheme") == "arithmetic")
            viscosity_averaging = arithmetic;
          else if (prm.get ("Viscosity averaging scheme") == "geometric")
            viscosity_averaging = geometric;
          else if (prm.get ("Viscosity averaging scheme") == "maximum composition")
            viscosity_averaging = maximum_composition;
          else
            AssertThrow(false, ExcMessage("Not a valid viscosity averaging scheme"));

          // Parse MulticomponentCompressible properties
          reference_temperatures = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Reference temperatures"))),
                                                                           n_fields,
                                                                           "Reference temperatures");
          reference_densities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Reference densities"))),
                                                                        n_fields,
                                                                        "Reference densities");
          reference_isothermal_compressibilities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Reference isothermal compressibilities"))),
                                                   n_fields,
                                                   "Reference isothermal compressibilities");
          reference_Kprimes = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Reference Kprimes"))),
                                                                      n_fields,
                                                                      "Reference Kprimes");
          reference_thermal_expansivities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Reference thermal expansivities"))),
                                                                                    n_fields,
                                                                                    "Reference thermal expansivities");
          reference_specific_heats = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Reference specific heats"))),
                                                                             n_fields,
                                                                             "Reference specific heats");
          viscosities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Viscosities"))),
                                                                n_fields,
                                                                "Viscosities");
          thermal_conductivities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Thermal conductivities"))),
                                                                           n_fields,
                                                                           "Thermal conductivities");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Declare dependencies on solution variables
      // Declare dependencies on solution variables
      this->model_dependence.thermal_conductivity = NonlinearDependence::compositional_fields;

      this->model_dependence.viscosity = NonlinearDependence::compositional_fields;

      this->model_dependence.density |= NonlinearDependence::temperature
                                        | NonlinearDependence::pressure
                                        | NonlinearDependence::compositional_fields;
      this->model_dependence.compressibility = NonlinearDependence::temperature
                                               | NonlinearDependence::pressure
                                               | NonlinearDependence::compositional_fields;
      this->model_dependence.specific_heat = NonlinearDependence::temperature
                                             | NonlinearDependence::pressure
                                             | NonlinearDependence::compositional_fields;
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(MulticomponentCompressible,
                                   "multicomponent compressible",
                                   "This compressible model can accept an arbitrary number of compositional fields,"
                                   " where each field represents a phase or rock type which can have completely"
                                   " different properties from the others. Each field has a"
                                   " pressure- and temperature-dependent density of the form:"
                                   "\\begin{align}"
                                   "  \\rho(p,T) = \\rho_0"
                                   "              \\left(1 + \\left(P - \\frac{\\alpha_0}{\\beta_{T0}} (T-T_0)"
                                   "              \\right) \\frac{K\'\\beta_{T0}}\\right) "
                                   "              ^{\\frac{1}{K\'}}"
                                   "\\end{align}"
                                   " This density parameterization is based on the Murnaghan equation of state"
                                   " (a constant value of $K\'$. The thermal term assumes a constant value of "
                                   " $\\alpha / \\beta_T$ (giving a thermal pressure that is linear with temperature)."
                                   " In the case where pressure is constant and $K\'=1$, this parameterization reduces"
                                   " to the one in the Simpler model."
                                   " The value of each compositional field is interpreted as a volume fraction at the"
                                   " reference conditions. If the sum of the fields is greater than one, they are renormalized."
                                   " If it is less than one, material properties for ``background mantle'' make up the rest."
                                   " When more than one field is present, the"
                                   " material properties are averaged arithmetically.  An exception is the viscosity,"
                                   " where the averaging should make more of a difference.  For this, the user selects"
                                   " between arithmetic, harmonic, geometric, or maximum composition averaging.")
  }
}
