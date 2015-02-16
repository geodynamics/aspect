/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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


#include <aspect/material_model/latent_heat.h>
#include <deal.II/base/parameter_handler.h>

using namespace dealii;

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    double
    LatentHeat<dim>::
    phase_function (const Point<dim> &position,
                    const double temperature,
                    const double pressure,
                    const int phase) const
    {
      // if we already have the adiabatic conditions, we can use them
      if (this->get_adiabatic_conditions().is_initialized())
        {
          // first, get the pressure at which the phase transition occurs normally
          // and get the pressure change in the range of the phase transition
          const Point<dim,double> transition_point = this->get_geometry_model().representative_point(transition_depths[phase]);
          const Point<dim,double> transition_plus_width = this->get_geometry_model().representative_point(transition_depths[phase] + transition_widths[phase]);
          const Point<dim,double> transition_minus_width = this->get_geometry_model().representative_point(transition_depths[phase] - transition_widths[phase]);
          const double transition_pressure = this->get_adiabatic_conditions().pressure(transition_point);
          const double pressure_width = 0.5 * (this->get_adiabatic_conditions().pressure(transition_plus_width)
                                               - this->get_adiabatic_conditions().pressure(transition_minus_width));

          // then calculate the deviation from the transition point (both in temperature
          // and in pressure)
          double pressure_deviation = pressure - transition_pressure
                                      - transition_slopes[phase] * (temperature - transition_temperatures[phase]);

          // last, calculate the percentage of material that has undergone the transition
          // (also in dependence of the phase transition width - this is an input parameter)
          double phase_func;
          // use delta function for width = 0
          if (transition_widths[phase]==0)
            (pressure_deviation > 0) ? phase_func = 1 : phase_func = 0;
          else
            phase_func = 0.5*(1.0 + std::tanh(pressure_deviation / pressure_width));
          return phase_func;
        }
      // if we do not have the adiabatic conditions, we have to use the depth instead
      // this is less precise, because we do not have the exact pressure gradient, instead we use pressure/depth
      // (this is for calculating e.g. the density in the adiabatic profile)
      else
        {
          double depth = this->get_geometry_model().depth(position);
          double depth_deviation = (pressure > 0
                                    ?
                                    depth - transition_depths[phase]
                                    - transition_slopes[phase] * (depth / pressure) * (temperature - transition_temperatures[phase])
                                    :
                                    depth - transition_depths[phase]
                                    - transition_slopes[phase] / (this->get_gravity_model().gravity_vector(position).norm() * reference_rho)
                                    * (temperature - transition_temperatures[phase]));
          double phase_func;
          // use delta function for width = 0
          if (transition_widths[phase]==0)
            (depth_deviation > 0) ? phase_func = 1 : phase_func = 0;
          else
            phase_func = 0.5*(1.0 + std::tanh(depth_deviation / transition_widths[phase]));
          return phase_func;
        }
    }

    template <int dim>
    double
    LatentHeat<dim>::
    phase_function_derivative (const Point<dim> &position,
                               const double temperature,
                               const double pressure,
                               const int phase) const
    {
      // we already should have the adiabatic conditions here
      AssertThrow (this->get_adiabatic_conditions().is_initialized(),
                   ExcMessage("need adiabatic conditions to incorporate phase transitions"));

      // first, get the pressure at which the phase transition occurs normally
      const Point<dim,double> transition_point = this->get_geometry_model().representative_point(transition_depths[phase]);
      const Point<dim,double> transition_plus_width = this->get_geometry_model().representative_point(transition_depths[phase] + transition_widths[phase]);
      const Point<dim,double> transition_minus_width = this->get_geometry_model().representative_point(transition_depths[phase] - transition_widths[phase]);
      const double transition_pressure = this->get_adiabatic_conditions().pressure(transition_point);
      const double pressure_width = 0.5 * (this->get_adiabatic_conditions().pressure(transition_plus_width)
                                           - this->get_adiabatic_conditions().pressure(transition_minus_width));

      // then calculate the deviation from the transition point (both in temperature
      // and in pressure)
      double pressure_deviation = pressure - transition_pressure
                                  - transition_slopes[phase] * (temperature - transition_temperatures[phase]);

      // last, calculate the analytical derivative of the phase function
      if (transition_widths[phase]==0)
        return 0;
      else
        return 0.5 / pressure_width * (1.0 - std::tanh(pressure_deviation / pressure_width)
                                       * std::tanh(pressure_deviation / pressure_width));
    }

    template <int dim>
    double
    LatentHeat<dim>::
    viscosity (const double temperature,
               const double pressure,
               const std::vector<double> &composition,       /*composition*/
               const SymmetricTensor<2,dim> &,
               const Point<dim> &position) const
    {
      const double delta_temp = temperature-reference_T;
      double temperature_dependence = std::max(std::min(std::exp(-thermal_viscosity_exponent*delta_temp/reference_T),1e2),1e-2);

      if (std::isnan(temperature_dependence))
        temperature_dependence = 1.0;

      double composition_dependence = 1.0;
      if ((composition_viscosity_prefactor != 1.0) && (composition.size() > 0))
        {
          //geometric interpolation
          return (pow(10, ((1-composition[0]) * log10(eta*temperature_dependence)
                           + composition[0] * log10(eta*composition_viscosity_prefactor*temperature_dependence))));
        }

      return composition_dependence * temperature_dependence * eta;
    }


    template <int dim>
    double
    LatentHeat<dim>::
    reference_viscosity () const
    {
      return eta;
    }

    template <int dim>
    double
    LatentHeat<dim>::
    reference_density () const
    {
      return reference_rho;
    }

    template <int dim>
    double
    LatentHeat<dim>::
    reference_thermal_expansion_coefficient () const
    {
      return thermal_alpha;
    }

    template <int dim>
    double
    LatentHeat<dim>::
    specific_heat (const double,
                   const double,
                   const std::vector<double> &, /*composition*/
                   const Point<dim> &) const
    {
      return reference_specific_heat;
    }

    template <int dim>
    double
    LatentHeat<dim>::
    reference_cp () const
    {
      return reference_specific_heat;
    }

    template <int dim>
    double
    LatentHeat<dim>::
    thermal_conductivity (const double,
                          const double,
                          const std::vector<double> &, /*composition*/
                          const Point<dim> &position) const
    {
      return k_value;
    }

    template <int dim>
    double
    LatentHeat<dim>::
    reference_thermal_diffusivity () const
    {
      return k_value/(reference_rho*reference_specific_heat);
    }

    template <int dim>
    double
    LatentHeat<dim>::
    density (const double temperature,
             const double pressure,
             const std::vector<double> &compositional_fields, /*composition*/
             const Point<dim> &position) const
    {
      // first, calculate temperature dependence of density
      double temperature_dependence = 1.0;
      if (this->include_adiabatic_heating ())
        {
          // temperature dependence is 1 - alpha * (T - T(adiabatic))
          if (this->get_adiabatic_conditions().is_initialized())
            temperature_dependence -= (temperature - this->get_adiabatic_conditions().temperature(position))
                                      * thermal_expansion_coefficient(temperature, pressure, compositional_fields, position);
        }
      else
        temperature_dependence -= temperature * thermal_expansion_coefficient(temperature, pressure, compositional_fields, position);

      // second, calculate composition dependence of density
      // constant density difference between peridotite and eclogite
      const double composition_dependence = compositional_fields.size()>0
                                            ?
                                            compositional_delta_rho * compositional_fields[0]
                                            :
                                            0.0;

      // third, calculate the density differences due to phase transitions (temperature-
      // and pressure dependence included)
      // the phase function gives the percentage of material that has
      // already undergone the phase transition to the higher-pressure material
      // (this is done individual for each transitions and summed up
      // in the end)
      // this means, that there are no actual density "jumps", but gradual
      // transition between the materials
      double phase_dependence = 0.0;
      if (compositional_fields.size()==0)      //only one field
        for (unsigned int i=0; i<transition_depths.size(); ++i)
          {
            const double phaseFunction = phase_function (position,
                                                         temperature,
                                                         pressure,
                                                         i);

            phase_dependence += phaseFunction * density_jumps[i];
          }
      else if (compositional_fields.size()>0)
        for (unsigned int i=0; i<transition_depths.size(); ++i)
          {
            const double phaseFunction = phase_function (position,
                                                         temperature,
                                                         pressure,
                                                         i);
            if (transition_phases[i] == 0)     // 1st compositional field
              phase_dependence += phaseFunction * density_jumps[i] * (1.0 - compositional_fields[0]);
            else if (transition_phases[i] == 1) // 2nd compositional field
              phase_dependence += phaseFunction * density_jumps[i] * compositional_fields[0];
          }

      // fourth, pressure dependence of density
      const double kappa = compressibility(temperature,pressure,compositional_fields,position);
      const double pressure_dependence = reference_rho * kappa * (pressure - this->get_surface_pressure());

      // in the end, all the influences are added up
      return (reference_rho + composition_dependence + pressure_dependence + phase_dependence) * temperature_dependence;
    }


    template <int dim>
    double
    LatentHeat<dim>::
    thermal_expansion_coefficient (const double temperature,
                                   const double,
                                   const std::vector<double> &, /*composition*/
                                   const Point<dim> &position) const
    {
      return thermal_alpha;
    }


    template <int dim>
    double
    LatentHeat<dim>::
    compressibility (const double,
                     const double,
                     const std::vector<double> &, /*composition*/
                     const Point<dim> &) const
    {
      return reference_compressibility;
    }


    template <int dim>
    double
    LatentHeat<dim>::
    entropy_derivative (const double temperature,
                        const double pressure,
                        const std::vector<double> &compositional_fields,
                        const Point<dim> &position,
                        const NonlinearDependence::Dependence dependence) const
    {
      double entropy_gradient = 0.0;
      const double rho = density (temperature, pressure, compositional_fields, position);
      if (this->get_adiabatic_conditions().is_initialized() && this->include_latent_heat())
        for (unsigned int phase=0; phase<transition_depths.size(); ++phase)
          {
            //calculate derivative of the phase function
            const double PhaseFunctionDerivative = phase_function_derivative(position,
                                                                             temperature,
                                                                             pressure,
                                                                             phase);

            // calculate the change of entropy across the phase transition
            double entropy_change = 0.0;
            if (compositional_fields.size()==0)      // only one compositional field
              entropy_change = transition_slopes[phase] * density_jumps[phase] / (rho * rho);
            else
              {
                if (transition_phases[phase] == 0)     // 1st compositional field
                  entropy_change = transition_slopes[phase] * density_jumps[phase] / (rho * rho) * (1.0 - compositional_fields[0]);
                else if (transition_phases[phase] == 1) // 2nd compositional field
                  entropy_change = transition_slopes[phase] * density_jumps[phase] / (rho * rho) * compositional_fields[0];
              }
            // we need DeltaS * DX/Dpressure_deviation for the pressure derivative
            // and - DeltaS * DX/Dpressure_deviation * gamma for the temperature derivative
            if (dependence == NonlinearDependence::pressure)
              entropy_gradient += PhaseFunctionDerivative * entropy_change;
            else if (dependence == NonlinearDependence::temperature)
              entropy_gradient -= PhaseFunctionDerivative * entropy_change * transition_slopes[phase];
            else
              AssertThrow(false, ExcMessage("not implemented"));
          }

      return entropy_gradient;
    }


    template <int dim>
    bool
    LatentHeat<dim>::
    viscosity_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      if ((dependence & NonlinearDependence::temperature) != NonlinearDependence::none)
        return true;
      else if ((dependence & NonlinearDependence::compositional_fields) != NonlinearDependence::none)
        return true;
      else
        return false;
    }


    template <int dim>
    bool
    LatentHeat<dim>::
    density_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      if ((dependence & NonlinearDependence::temperature) != NonlinearDependence::none)
        return true;
      else if ((dependence & NonlinearDependence::pressure) != NonlinearDependence::none)
        return true;
      else if ((dependence & NonlinearDependence::compositional_fields) != NonlinearDependence::none)
        return true;
      else
        return false;
    }

    template <int dim>
    bool
    LatentHeat<dim>::
    compressibility_depends_on (const NonlinearDependence::Dependence) const
    {
      return false;
    }

    template <int dim>
    bool
    LatentHeat<dim>::
    specific_heat_depends_on (const NonlinearDependence::Dependence) const
    {
      return false;
    }

    template <int dim>
    bool
    LatentHeat<dim>::
    thermal_conductivity_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      return false;
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
          prm.declare_entry ("Reference density", "3300",
                             Patterns::Double (0),
                             "Reference density $\\rho_0$. Units: $kg/m^3$.");
          prm.declare_entry ("Reference temperature", "293",
                             Patterns::Double (0),
                             "The reference temperature $T_0$. Units: $K$.");
          prm.declare_entry ("Viscosity", "5e24",
                             Patterns::Double (0),
                             "The value of the constant viscosity. Units: $kg/m/s$.");
          prm.declare_entry ("Composition viscosity prefactor", "1.0",
                             Patterns::Double (0),
                             "A linear dependency of viscosity on composition. Dimensionless prefactor.");
          prm.declare_entry ("Thermal viscosity exponent", "0.0",
                             Patterns::Double (0),
                             "The temperature dependence of viscosity. Dimensionless exponent.");
          prm.declare_entry ("Thermal conductivity", "2.38",
                             Patterns::Double (0),
                             "The value of the thermal conductivity $k$. "
                             "Units: $W/m/K$.");
          prm.declare_entry ("Reference specific heat", "1250",
                             Patterns::Double (0),
                             "The value of the specific heat $cp$. "
                             "Units: $J/kg/K$.");
          prm.declare_entry ("Thermal expansion coefficient", "4e-5",
                             Patterns::Double (0),
                             "The value of the thermal expansion coefficient $\\beta$. "
                             "Units: $1/K$.");
          prm.declare_entry ("Compressibility", "5.124e-12",
                             Patterns::Double (0),
                             "The value of the compressibility $\\kappa$. "
                             "Units: $1/Pa$.");
          prm.declare_entry ("Density differential for compositional field 1", "0",
                             Patterns::Double(),
                             "If compositional fields are used, then one would frequently want "
                             "to make the density depend on these fields. In this simple material "
                             "model, we make the following assumptions: if no compositional fields "
                             "are used in the current simulation, then the density is simply the usual "
                             "one with its linear dependence on the temperature. If there are compositional "
                             "fields, then the density only depends on the first one in such a way that "
                             "the density has an additional term of the kind $+\\Delta \\rho \\; c_1(\\mathbf x)$. "
                             "This parameter describes the value of $\\Delta \\rho$. Units: $kg/m^3/\\textrm{unit "
                             "change in composition}$.");
          prm.declare_entry ("Phase transition depths", "",
                             Patterns::List (Patterns::Double(0)),
                             "A list of depths where phase transitions occur. Values must "
                             "monotonically increase. "
                             "Units: $m$.");
          prm.declare_entry ("Phase transition widths", "",
                             Patterns::List (Patterns::Double(0)),
                             "A list of widths for each phase transition. The phase functions "
                             "are scaled with these values, leading to a jump betwen phases "
                             "for a value of zero and a gradual transition for larger values. "
                             "List must have the same number of entries as Phase transition depths. "
                             "Units: $m$.");
          prm.declare_entry ("Phase transition temperatures", "",
                             Patterns::List (Patterns::Double(0)),
                             "A list of temperatures where phase transitions occur. Higher or lower "
                             "temperatures lead to phase transition ocurring in smaller or greater "
                             "depths than given in Phase transition depths, depending on the "
                             "Clapeyron slope given in Phase transition Clapeyron slopes. "
                             "List must have the same number of entries as Phase transition depths. "
                             "Units: $K$.");
          prm.declare_entry ("Phase transition Clapeyron slopes", "",
                             Patterns::List (Patterns::Double()),
                             "A list of Clapeyron slopes for each phase transition. A positive "
                             "Clapeyron slope indicates that the phase transition will occur in "
                             "a greater depth, if the temperature is higher than the one given in "
                             "Phase transition temperatures and in a smaller depth, if the "
                             "temperature is smaller than the one given in Phase transition temperatures. "
                             "For negative slopes the other way round. "
                             "List must have the same number of entries as Phase transition depths. "
                             "Units: $Pa/K$.");
          prm.declare_entry ("Phase transition density jumps", "",
                             Patterns::List (Patterns::Double(0)),
                             "A list of density jumps at each phase transition. A positive value means "
                             "that the density increases with depth. The corresponding entry in "
                             "Corresponding phase for density jump determines if the density jump occurs "
                             "in peridotite, eclogite or none of them."
                             "List must have the same number of entries as Phase transition depths. "
                             "Units: $kg/m^3$.");
          prm.declare_entry ("Corresponding phase for density jump", "",
                             Patterns::List (Patterns::Integer(0)),
                             "A list of phases, which correspond to the Phase transition density jumps. "
                             "The density jumps occur only in the phase that is given by this phase value. "
                             "0 stands for the 1st compositional fields, 1 for the second compositional field "
                             "and -1 for none of them. "
                             "List must have the same number of entries as Phase transition depths. "
                             "Units: $Pa/K$.");
          prm.declare_entry ("Viscosity prefactors", "",
                             Patterns::List (Patterns::Double(0)),
                             "A list of prefactors for the viscosity for each phase. The reference "
                             "viscosity will be multiplied by this factor to get the corresponding "
                             "viscosity for each phase. "
                             "List must have one more entry than Phase transition depths. "
                             "Units: non-dimensional.");
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

          transition_depths = Utilities::string_to_double
                              (Utilities::split_string_list(prm.get ("Phase transition depths")));
          transition_widths= Utilities::string_to_double
                             (Utilities::split_string_list(prm.get ("Phase transition widths")));
          transition_temperatures = Utilities::string_to_double
                                    (Utilities::split_string_list(prm.get ("Phase transition temperatures")));
          transition_slopes = Utilities::string_to_double
                              (Utilities::split_string_list(prm.get ("Phase transition Clapeyron slopes")));
          density_jumps = Utilities::string_to_double
                          (Utilities::split_string_list(prm.get ("Phase transition density jumps")));
          transition_phases = Utilities::string_to_int
                              (Utilities::split_string_list(prm.get ("Corresponding phase for density jump")));
          phase_prefactors = Utilities::string_to_double
                             (Utilities::split_string_list(prm.get ("Viscosity prefactors")));

          if (transition_widths.size() != transition_depths.size() ||
              transition_temperatures.size() != transition_depths.size() ||
              transition_slopes.size() != transition_depths.size() ||
              density_jumps.size() != transition_depths.size() ||
              transition_phases.size() != transition_depths.size() ||
              phase_prefactors.size() != transition_depths.size()+1)
            AssertThrow(false, ExcMessage("Error: At least one list that gives input parameters for the phase transitions has the wrong size."));

          if (thermal_viscosity_exponent!=0.0 && reference_T == 0.0)
            AssertThrow(false, ExcMessage("Error: Material model latent heat with Thermal viscosity exponent can not have reference_T=0."));
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
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
                                   "of the latent heat terms following Christensen \\& Yuen, 1986. "
                                   "The change of entropy is calculated as "
                                   "$Delta S = \\gamma \\frac{\\Delta\\rho}{\\rho^2}$ with the "
                                   "Clapeyron slope $\\gamma$ and the density change $\\Delta\\rho$ "
                                   "of the phase transition being input parameters. "
                                   "The model employs an analytic phase function in the form "
                                   "$X=0.5 \\left( 1 + \\tanh \\left( \\frac{\\Delta p}{\\Delta p_0} \\right) \\right)$ "
                                   "with $\\Delta p = p - p_{transition} - \\gamma \\left( T - T_{transition} \\right)$ "
                                   "and $\\Delta p_0$ being the pressure difference over the width "
                                   "of the phase transition (specified as input parameter).")
  }
}
