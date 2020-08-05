/*
  Copyright (C) 2011 - 2020 by the authors of the ASPECT code.

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


#include <aspect/material_model/latent_heat.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/gravity_model/interface.h>


namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    void
    LatentHeat<dim>::
    evaluate(const MaterialModelInputs<dim> &in,
             MaterialModelOutputs<dim> &out) const
    {
      for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
        {
          const double temperature = in.temperature[i];
          const double pressure = in.pressure[i];
          const std::vector<double> composition = in.composition[i];
          const Point<dim> position = in.position[i];

          // Assign constant material properties
          {
            out.specific_heat[i] = reference_specific_heat;
            out.thermal_conductivities[i] = k_value;
            out.thermal_expansion_coefficients[i] = thermal_alpha;
            out.compressibilities[i] = reference_compressibility;
          }

          // Calculate Viscosity
          {
            const double reference_temperature = (this->include_adiabatic_heating()
                                                  ?
                                                  this->get_adiabatic_conditions().temperature(in.position[i])
                                                  :
                                                  reference_T);

            const double delta_temp = temperature-reference_temperature;
            const double T_dependence = ( thermal_viscosity_exponent == 0.0
                                          ?
                                          0.0
                                          :
                                          thermal_viscosity_exponent * delta_temp / reference_temperature );

            double visc_temperature_dependence = std::max(std::min(std::exp(-T_dependence),1e2),1e-2);

            if (std::isnan(visc_temperature_dependence))
              visc_temperature_dependence = 1.0;

            double visc_composition_dependence = 1.0;
            if ((composition_viscosity_prefactor != 1.0) && (composition.size() > 0))
              {
                // geometric interpolation
                out.viscosities[i] = (pow(10, ((1-composition[0]) * log10(eta*visc_temperature_dependence)
                                               + composition[0] * log10(eta*composition_viscosity_prefactor*visc_temperature_dependence))));
              }
            else
              out.viscosities[i] = visc_composition_dependence * visc_temperature_dependence * eta;
          }

          // Calculate density
          // and phase dependence of viscosity
          {
            // first, calculate temperature dependence of density
            double density_temperature_dependence = 1.0;
            if (this->include_adiabatic_heating ())
              {
                // temperature dependence is 1 - alpha * (T - T(adiabatic))
                density_temperature_dependence -= (temperature - this->get_adiabatic_conditions().temperature(position))
                                                  * thermal_alpha;
              }
            else
              density_temperature_dependence -= temperature * thermal_alpha;

            // second, calculate composition dependence of density
            // constant density difference between peridotite and eclogite
            const double density_composition_dependence = composition.size()>0
                                                          ?
                                                          compositional_delta_rho * composition[0]
                                                          :
                                                          0.0;

            // third, calculate the density (and viscosity) differences due to phase
            // transitions (temperature- and pressure dependence included).
            // the phase function gives the percentage of material that has
            // already undergone the phase transition to the higher-pressure material
            // (this is done individual for each transitions and summed up
            // in the end)
            // this means, that there are no actual density or viscosity "jumps", but
            // gradual transitions between the materials
            double phase_dependence = 0.0;
            double viscosity_phase_dependence = 1.0;


            // Loop through phase transitions
            for (unsigned int phase=0; phase<phase_function.n_phase_transitions(); ++phase)
              {
                const double depth = this->get_geometry_model().depth(position);
                const double pressure_depth_derivative = (depth > 0.0)
                                                         ?
                                                         pressure / depth
                                                         :
                                                         this->get_gravity_model().gravity_vector(in.position[i]).norm() * reference_rho;

                const MaterialUtilities::PhaseFunctionInputs<dim> phase_in(temperature,
                                                                           pressure,
                                                                           depth,
                                                                           pressure_depth_derivative,
                                                                           phase);

                const double phaseFunction = phase_function.compute_value(phase_in);

                // Note that for the densities we have a list of jumps, so the index used
                // in the loop corresponds to the index of the phase transition. For the
                // viscosities we have a list of prefactors, which has one more entry
                // for the first layer, so we have to use phase+1 as the index.
                if (composition.size()==0)
                  {
                    phase_dependence += phaseFunction * density_jumps[phase];
                    viscosity_phase_dependence *= 1. + phaseFunction * (phase_prefactors[phase+1]-1.);
                  }
                else if (composition.size()>0)
                  {
                    if (transition_phases[phase] == 0)     // 1st compositional field
                      phase_dependence += phaseFunction * density_jumps[phase] * (1.0 - composition[0]);
                    else if (transition_phases[phase] == 1) // 2nd compositional field
                      phase_dependence += phaseFunction * density_jumps[phase] * composition[0];

                    viscosity_phase_dependence *= 1. + phaseFunction * (phase_prefactors[phase]-1.);
                  }
              }

            // fourth, pressure dependence of density
            const double kappa = reference_compressibility;
            const double pressure_dependence = reference_rho * kappa * (pressure - this->get_surface_pressure());

            // in the end, all the influences are added up
            out.densities[i] = (reference_rho + density_composition_dependence + pressure_dependence + phase_dependence)
                               * density_temperature_dependence;
            out.viscosities[i] = std::max(min_viscosity, std::min(max_viscosity, out.viscosities[i] * viscosity_phase_dependence));
          }

          // Calculate entropy derivative
          {
            double entropy_gradient_pressure = 0.0;
            double entropy_gradient_temperature = 0.0;
            const double rho = out.densities[i];

            if (this->get_adiabatic_conditions().is_initialized() && this->include_latent_heat())
              for (unsigned int phase=0; phase<phase_function.n_phase_transitions(); ++phase)
                {
                  const double depth = this->get_geometry_model().depth(in.position[i]);
                  const double pressure_depth_derivative = (pressure > 0)
                                                           ?
                                                           depth / pressure
                                                           :
                                                           this->get_gravity_model().gravity_vector(in.position[i]).norm() * reference_rho;

                  const MaterialUtilities::PhaseFunctionInputs<dim> phase_in(temperature,
                                                                             pressure,
                                                                             depth,
                                                                             pressure_depth_derivative,
                                                                             phase);

                  const double PhaseFunctionDerivative = phase_function.compute_derivative(phase_in);
                  const double clapeyron_slope = phase_function.get_transition_slope(phase);

                  double entropy_change = 0.0;
                  if (composition.size()==0)      // only one compositional field
                    entropy_change = clapeyron_slope * density_jumps[phase] / (rho * rho);
                  else
                    {
                      if (transition_phases[phase] == 0)     // 1st compositional field
                        entropy_change = clapeyron_slope * density_jumps[phase] / (rho * rho) * (1.0 - composition[0]);
                      else if (transition_phases[phase] == 1) // 2nd compositional field
                        entropy_change = clapeyron_slope * density_jumps[phase] / (rho * rho) * composition[0];
                    }
                  // we need DeltaS * DX/Dpressure_deviation for the pressure derivative
                  // and - DeltaS * DX/Dpressure_deviation * gamma for the temperature derivative
                  entropy_gradient_pressure += PhaseFunctionDerivative * entropy_change;
                  entropy_gradient_temperature -= PhaseFunctionDerivative * entropy_change * clapeyron_slope;
                }

            out.entropy_derivative_pressure[i] = entropy_gradient_pressure;
            out.entropy_derivative_temperature[i] = entropy_gradient_temperature;
          }

          // Assign reaction terms
          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            out.reaction_terms[i][c] = 0.0;

        }
    }


    template <int dim>
    double
    LatentHeat<dim>::
    reference_viscosity () const
    {
      return eta;
    }




    template <int dim>
    bool
    LatentHeat<dim>::
    is_compressible () const
    {
      if (reference_compressibility > 0)
        return true;
      else
        return false;
    }



    template <int dim>
    void
    LatentHeat<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Latent heat");
        {
          prm.declare_entry ("Reference density", "3300.",
                             Patterns::Double (0.),
                             "Reference density $\\rho_0$. "
                             "Units: \\si{\\kilogram\\per\\meter\\cubed}.");
          prm.declare_entry ("Reference temperature", "293.",
                             Patterns::Double (0.),
                             "The reference temperature $T_0$. Units: \\si{\\kelvin}.");
          prm.declare_entry ("Viscosity", "5e24",
                             Patterns::Double (0.),
                             "The value of the constant viscosity. "
                             "Units: \\si{\\pascal\\second}.");
          prm.declare_entry ("Composition viscosity prefactor", "1.0",
                             Patterns::Double (0.),
                             "A linear dependency of viscosity on composition. Dimensionless prefactor.");
          prm.declare_entry ("Thermal viscosity exponent", "0.0",
                             Patterns::Double (0.),
                             "The temperature dependence of viscosity. Dimensionless exponent.");
          prm.declare_entry ("Thermal conductivity", "2.38",
                             Patterns::Double (0.),
                             "The value of the thermal conductivity $k$. "
                             "Units: \\si{\\watt\\per\\meter\\per\\kelvin}.");
          prm.declare_entry ("Reference specific heat", "1250.",
                             Patterns::Double (0.),
                             "The value of the specific heat $C_p$. "
                             "Units: \\si{\\joule\\per\\kelvin\\per\\kilogram}.");
          prm.declare_entry ("Thermal expansion coefficient", "4e-5",
                             Patterns::Double (0.),
                             "The value of the thermal expansion coefficient $\\beta$. "
                             "Units: \\si{\\per\\kelvin}.");
          prm.declare_entry ("Compressibility", "5.124e-12",
                             Patterns::Double (0.),
                             "The value of the compressibility $\\kappa$. "
                             "Units: \\si{\\per\\pascal}.");
          prm.declare_entry ("Density differential for compositional field 1", "0.",
                             Patterns::Double(),
                             "If compositional fields are used, then one would frequently want "
                             "to make the density depend on these fields. In this simple material "
                             "model, we make the following assumptions: if no compositional fields "
                             "are used in the current simulation, then the density is simply the usual "
                             "one with its linear dependence on the temperature. If there are compositional "
                             "fields, then the density only depends on the first one in such a way that "
                             "the density has an additional term of the kind $+\\Delta \\rho \\; c_1(\\mathbf x)$. "
                             "This parameter describes the value of $\\Delta \\rho$. "
                             "Units: \\si{\\kilogram\\per\\meter\\cubed}/unit change in composition.");
          prm.declare_entry ("Phase transition density jumps", "",
                             Patterns::List (Patterns::Double (0.)),
                             "A list of density jumps at each phase transition. A positive value means "
                             "that the density increases with depth. The corresponding entry in "
                             "Corresponding phase for density jump determines if the density jump occurs "
                             "in peridotite, eclogite or none of them."
                             "List must have the same number of entries as Phase transition depths. "
                             "Units: \\si{\\kilogram\\per\\meter\\cubed}.");
          prm.declare_entry ("Corresponding phase for density jump", "",
                             Patterns::List (Patterns::Integer(0)),
                             "A list of phases, which correspond to the Phase transition density jumps. "
                             "The density jumps occur only in the phase that is given by this phase value. "
                             "0 stands for the 1st compositional fields, 1 for the second compositional field "
                             "and -1 for none of them. "
                             "List must have the same number of entries as Phase transition depths. "
                             "Units: \\si{\\pascal\\per\\kelvin}.");
          prm.declare_entry ("Viscosity prefactors", "",
                             Patterns::List (Patterns::Double (0.)),
                             "A list of prefactors for the viscosity for each phase. The reference "
                             "viscosity will be multiplied by this factor to get the corresponding "
                             "viscosity for each phase. "
                             "List must have one more entry than Phase transition depths. "
                             "Units: non-dimensional.");
          prm.declare_entry ("Minimum viscosity", "1e19",
                             Patterns::Double (0.),
                             "Limit for the minimum viscosity in the model. "
                             "Units: Pa \\, s.");
          prm.declare_entry ("Maximum viscosity", "1e24",
                             Patterns::Double (0.),
                             "Limit for the maximum viscosity in the model. "
                             "Units: Pa \\, s.");

          MaterialUtilities::PhaseFunction<dim>::declare_parameters(prm);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    LatentHeat<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Latent heat");
        {
          reference_rho              = prm.get_double ("Reference density");
          reference_T                = prm.get_double ("Reference temperature");
          eta                        = prm.get_double ("Viscosity");
          composition_viscosity_prefactor = prm.get_double ("Composition viscosity prefactor");
          thermal_viscosity_exponent = prm.get_double ("Thermal viscosity exponent");
          k_value                    = prm.get_double ("Thermal conductivity");
          reference_specific_heat    = prm.get_double ("Reference specific heat");
          thermal_alpha              = prm.get_double ("Thermal expansion coefficient");
          reference_compressibility  = prm.get_double ("Compressibility");
          compositional_delta_rho    = prm.get_double ("Density differential for compositional field 1");
          min_viscosity              = prm.get_double ("Minimum viscosity");
          max_viscosity              = prm.get_double ("Maximum viscosity");

          phase_function.initialize_simulator (this->get_simulator());
          phase_function.parse_parameters (prm);

          density_jumps = Utilities::string_to_double
                          (Utilities::split_string_list(prm.get ("Phase transition density jumps")));
          transition_phases = Utilities::string_to_int
                              (Utilities::split_string_list(prm.get ("Corresponding phase for density jump")));
          phase_prefactors = Utilities::string_to_double
                             (Utilities::split_string_list(prm.get ("Viscosity prefactors")));

          const unsigned int n_transitions = phase_function.n_phase_transitions();
          if (density_jumps.size() != n_transitions ||
              transition_phases.size() != n_transitions ||
              phase_prefactors.size() != n_transitions+1)
            AssertThrow(false, ExcMessage("Error: At least one list that provides input parameters for phase "
                                          "transitions has the wrong size. The phase function object reports that "
                                          "there are " + std::to_string(n_transitions) + " transitions, "
                                          "therefore the material model expects " + std::to_string(n_transitions) +
                                          " density jumps and corresponding phases, and "
                                          + std::to_string(n_transitions+1) + " viscosity prefactors."));

          // as the phase viscosity prefactors are all applied multiplicatively on top of each other,
          // we have to scale them here so that they are relative factors in comparison to the product
          // of the prefactors of all phase above the current one
          for (unsigned int phase=1; phase<phase_prefactors.size(); ++phase)
            {
              phase_prefactors[phase] /= phase_prefactors[phase-1];
            }

          if (thermal_viscosity_exponent!=0.0 && reference_T == 0.0)
            AssertThrow(false, ExcMessage("Error: Material model latent heat with Thermal viscosity exponent can not have reference_T=0."));

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Declare dependencies on solution variables
      this->model_dependence.viscosity = NonlinearDependence::temperature | NonlinearDependence::compositional_fields;
      this->model_dependence.density = NonlinearDependence::temperature | NonlinearDependence::pressure | NonlinearDependence::compositional_fields;
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
    ASPECT_REGISTER_MATERIAL_MODEL(LatentHeat,
                                   "latent heat",
                                   "A material model that includes phase transitions "
                                   "and the possibility that latent heat is released "
                                   "or absorbed when material crosses one of the "
                                   "phase transitions of up to two different materials "
                                   "(compositional fields). "
                                   "This model implements a standard approximation "
                                   "of the latent heat terms following Christensen \\& Yuen, 1985. "
                                   "The change of entropy is calculated as "
                                   "$Delta S = \\gamma \\frac{\\Delta\\rho}{\\rho^2}$ with the "
                                   "Clapeyron slope $\\gamma$ and the density change $\\Delta\\rho$ "
                                   "of the phase transition being input parameters. "
                                   "The model employs an analytic phase function in the form "
                                   "$X=\\frac{1}{2} \\left( 1 + \\tanh \\left( \\frac{\\Delta p}{\\Delta p_0} \\right) \\right)$ "
                                   "with $\\Delta p = p - p_{\\text{transition}} - \\gamma \\left( T - T_{\\text{transition}} \\right)$ "
                                   "and $\\Delta p_0$ being the pressure difference over the width "
                                   "of the phase transition (specified as input parameter).")
  }
}
