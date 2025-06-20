/*
  Copyright (C) 2019 - 2025 by the authors of the ASPECT code.

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


#include <aspect/material_model/rheology/grain_boundary_sliding.h>
#include <aspect/material_model/utilities.h>
#include <aspect/utilities.h>

#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/parameter_handler.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace Rheology
    {
      GrainBoundarySlidingParameters::GrainBoundarySlidingParameters()
        : prefactor (numbers::signaling_nan<double>()),
          activation_energy (numbers::signaling_nan<double>()),
          activation_volume (numbers::signaling_nan<double>()),
          stress_exponent (numbers::signaling_nan<double>())
      {}



      template <int dim>
      GrainBoundarySliding<dim>::GrainBoundarySliding()
        = default;


      template <int dim>
      const GrainBoundarySlidingParameters
      GrainBoundarySliding<dim>::compute_slide_parameters (const unsigned int composition,
                                                           const std::vector<double> &phase_function_values,
                                                           const std::vector<unsigned int> &n_phase_transitions_per_composition) const
      {
        GrainBoundarySlidingParameters slide_parameters;
        if (phase_function_values == std::vector<double>())
          {
            // no phases
            slide_parameters.prefactor = prefactors[composition];
            slide_parameters.activation_energy = activation_energies[composition];
            slide_parameters.activation_volume = activation_volumes[composition];
            slide_parameters.stress_exponent = stress_exponents[composition];
            slide_parameters.grain_size_exponent = grain_size_exponents[composition];
          }
        else
          {
            // Average among phases
            slide_parameters.prefactor = MaterialModel::MaterialUtilities::phase_average_value(phase_function_values, n_phase_transitions_per_composition,
                                         prefactors, composition,  MaterialModel::MaterialUtilities::PhaseUtilities::logarithmic);
            slide_parameters.activation_energy = MaterialModel::MaterialUtilities::phase_average_value(phase_function_values, n_phase_transitions_per_composition,
                                                 activation_energies, composition);
            slide_parameters.activation_volume = MaterialModel::MaterialUtilities::phase_average_value(phase_function_values, n_phase_transitions_per_composition,
                                                 activation_volumes, composition);
            slide_parameters.stress_exponent = MaterialModel::MaterialUtilities::phase_average_value(phase_function_values, n_phase_transitions_per_composition,
                                               stress_exponents, composition);
            slide_parameters.grain_size_exponent = MaterialModel::MaterialUtilities::phase_average_value(phase_function_values, n_phase_transitions_per_composition,
                                                   grain_size_exponents, composition);
          }
        return slide_parameters;
      }



      template <int dim>
      double
      GrainBoundarySliding<dim>::compute_viscosity (const double strain_rate,
                                                    const double pressure,
                                                    const double temperature,
                                                    const unsigned int composition,
                                                    const std::vector<double> &phase_function_values,
                                                    const std::vector<unsigned int> &n_phase_transitions_per_composition) const
      {
        return compute_viscosity(strain_rate, pressure, temperature, fixed_grain_size, composition, phase_function_values, n_phase_transitions_per_composition);
      }


      template <int dim>
      double
      GrainBoundarySliding<dim>::compute_viscosity (const double strain_rate,
                                                    const double pressure,
                                                    const double temperature,
                                                    const double grain_size,
                                                    const unsigned int composition,
                                                    const std::vector<double> &phase_function_values,
                                                    const std::vector<unsigned int> &n_phase_transitions_per_composition) const
      {
        const GrainBoundarySlidingParameters p = compute_slide_parameters(composition,
                                                                          phase_function_values,
                                                                          n_phase_transitions_per_composition);

        // Power law grain-size sensitive creep equation:
        // viscosity = 0.5 * A^(-1/n) * d^m * exp[(E + PV)/(nRT)] * edot_ii^[(1 - n)/n]
        // A: prefactor, edot_ii: square root of second invariant of deviatoric strain rate tensor,
        // E: activation energy, d: grain size, m: grain size exponenet, P: pressure,
        // V; activation volume, n: stress exponent, R: gas constant, T: temperature.
        double viscosity_grainboundarysliding = 0.5 * std::pow(p.prefactor, -1.0 / p.stress_exponent) *
                                                std::pow(grain_size, p.grain_size_exponent) *
                                                std::exp((p.activation_energy + pressure * p.activation_volume) /
                                                         (constants::gas_constant * temperature * p.stress_exponent)) *
                                                std::pow(strain_rate, (1.0 - p.stress_exponent) / p.stress_exponent);

        Assert (viscosity_grainboundarysliding > 0.0,
                ExcMessage ("Negative grain boundary sliding viscosity detected. This is unphysical and should not happen. "
                            "Check for negative parameters. Temperature and pressure are "
                            + Utilities::to_string(temperature) + " K, " + Utilities::to_string(pressure) + " Pa. "));

        // In ice other deformation mechanisms become dominant at lower or higher strain rates
        // and temperatures so it is therefore both reasonable and desirable to require
        // the single-mechanism viscosity to be smaller than std::sqrt(max_double).
        viscosity_grainboundarysliding = std::min(viscosity_grainboundarysliding, std::sqrt(std::numeric_limits<double>::max()));

        return viscosity_grainboundarysliding;
      }

      template <int dim>
      void
      GrainBoundarySliding<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Grain boundary sliding");
        {
          prm.declare_entry ("Prefactors for grain boundary sliding", "3.9e-19.2",
                             Patterns::Anything(),
                             "Here we use the default values for ice as given in Goldsby & Kohlstedt, 2001");
          prm.declare_entry ("Stress exponents for grain boundary sliding", "1.8",
                             Patterns::List(Patterns::Double(0.)),
                             "Here we use the default values for ice as given in Goldsby & Kohlstedt, 2001.");
          prm.declare_entry ("Grain size exponents for grain boundary sliding", "1.4.",
                             Patterns::Anything(),
                             "Here we use the default values for ice as given in Goldsby & Kohlstedt, 2001.");
          prm.declare_entry ("Activation energies for grain boundary sliding", "49",
                             Patterns::Anything(),
                             "Here we use the default values for ice as given in Goldsby & Kohlstedt, 2001"
                             "for T > 262 k."
                             "Units: \\si{\\joule\\per\\mole}.");
          prm.declare_entry ("Activation volumes for grain boundary sliding", "13e-6",
                             Patterns::Anything(),
                             "Here we use the default values for ice as given in Goldsby & Kohlstedt, 2001."
                             "Units: \\si{\\meter\\cubed\\per\\mole}.");
          prm.declare_entry ("Grain size", "26e-6", Patterns::Double (0.),
                             "Here we use the default values for ice as given in Goldsby & Kohlstedt, 2001."
                             "Units: \\si{\\meter}.");
          prm.leave_subsection();
        }
      }



      template <int dim>
      void
      GrainBoundarySliding<dim>::parse_parameters (ParameterHandler &prm,
                                                   const std::unique_ptr<std::vector<unsigned int>> &expected_n_phases_per_composition)
      {
        prm.enter_subsection("Grain boundary sliding");
        {
          // Retrieve the list of composition names
          std::vector<std::string> compositional_field_names = this->introspection().get_composition_names();

          // Retrieve the list of names of fields that represent chemical compositions, and not, e.g.,
          // plastic strain
          std::vector<std::string> chemical_field_names = this->introspection().chemical_composition_field_names();

          // Establish that a background field is required here
          compositional_field_names.insert(compositional_field_names.begin(), "background");
          chemical_field_names.insert(chemical_field_names.begin(), "background");

          // Make options file for parsing maps to double arrays
          Utilities::MapParsing::Options options(chemical_field_names, "Prefactors for grain boundary sliding");
          options.list_of_allowed_keys = compositional_field_names;
          options.allow_multiple_values_per_key = true;
          if (expected_n_phases_per_composition)
            {
              options.n_values_per_key = *expected_n_phases_per_composition;

              // check_values_per_key is required to be true to duplicate single values
              // if they are to be used for all phases associated with a given key.
              options.check_values_per_key = true;
            }

          // Read parameters, each of size of number of composition + number of phases + 1
          prefactors = Utilities::MapParsing::parse_map_to_double_array(prm.get("Prefactors for grain boundary sliding"),
                                                                        options);

          options.property_name = "Stress exponents for grain boundary sliding";
          stress_exponents = Utilities::MapParsing::parse_map_to_double_array(prm.get("Stress exponents for grain boundary sliding"),
                                                                              options);

          options.property_name = "Grain size exponents for grain boundary sliding";
          grain_size_exponents = Utilities::MapParsing::parse_map_to_double_array(prm.get("Grain size exponents for grain boundary sliding"),
                                                                                  options);

          options.property_name = "Activation energies for grain boundary sliding";
          activation_energies = Utilities::MapParsing::parse_map_to_double_array(prm.get("Activation energies for grain boundary sliding"),
                                                                                 options);

          options.property_name = "Activation volumes for grain boundary sliding";
          activation_volumes = Utilities::MapParsing::parse_map_to_double_array(prm.get("Activation volumes for grain boundary sliding"),
                                                                                options);

          fixed_grain_size = prm.get_double("Grain size");

          // Check that there are no entries set to zero,
          // for example because the entry is for a field
          // that is masked anyway, like strain. Despite
          // these compositions being masked, their viscosities
          // are computed anyway and this will lead to division by zero.
          for (const double prefactor : prefactors)
            AssertThrow(prefactor > 0.,
                        ExcMessage("The grain boundary sliding prefactor should be larger than zero."));
          prm.leave_subsection();
        }
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
    template class GrainBoundarySliding<dim>; \
  }

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
