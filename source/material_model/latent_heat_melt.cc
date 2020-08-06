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


#include <aspect/material_model/latent_heat_melt.h>
#include <aspect/adiabatic_conditions/interface.h>


namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    void
    LatentHeatMelt<dim>::

    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
        {

          out.entropy_derivative_pressure[i] = 0;
          out.entropy_derivative_temperature[i] = 0;

          for (unsigned int c=0; c < in.composition[i].size(); ++c )
            {
              out.reaction_terms[i][c]            = 0;
            }

          std::vector<double> composition = in.composition[i];

          const double reference_temperature = (this->include_adiabatic_heating()
                                                ?
                                                this->get_adiabatic_conditions().temperature(in.position[i])
                                                :
                                                reference_T);

          const double delta_temp = in.temperature[i] - reference_temperature;
          const double T_dependence = ( thermal_viscosity_exponent == 0.0
                                        ?
                                        0.0
                                        :
                                        thermal_viscosity_exponent * delta_temp / reference_temperature );


          double temperature_dependence = std::max (std::min ( std::exp( - T_dependence ), 1e2 ), 1e-2 );

          if (std::isnan(temperature_dependence))
            temperature_dependence = 1.0;

          double composition_dependence = 1.0;
          if ((composition_viscosity_prefactor != 1.0) && (composition.size() > 0))
            {
              //geometric interpolation
              out.viscosities[i] = (pow (10, ( (1 - composition[0] ) * log10 ( eta * temperature_dependence )
                                               + composition[0] * log10 ( eta * composition_viscosity_prefactor * temperature_dependence))));
            }
          else
            {
              out.viscosities[i] = composition_dependence * temperature_dependence * eta;
            }

          out.specific_heat[i] = reference_specific_heat;
          out.thermal_conductivities[i] = k_value;
          out.thermal_expansion_coefficients[i] = thermal_alpha;
          out.compressibilities[i] = reference_compressibility;

          temperature_dependence = 1.0;

          if (this->include_adiabatic_heating ())
            {
              // temperature dependence is 1 - alpha * (T - T(adiabatic))
              temperature_dependence -= (in.temperature[i] - this->get_adiabatic_conditions().temperature(in.position[i])) * thermal_alpha;
            }
          else
            {
              temperature_dependence -= in.temperature[i] * thermal_alpha;
            }
          // second, calculate composition dependence of density
          // constant density difference between peridotite and eclogite
          composition_dependence = composition.size()>0
                                   ?
                                   compositional_delta_rho * composition[0]
                                   :
                                   0.0;

          // third, pressure dependence of density
          const double kappa = reference_compressibility;
          const double pressure_dependence = reference_rho * kappa * (in.pressure[i] - this->get_surface_pressure());

          // fourth, melt fraction dependence
          double melt_dependence = (1.0 - relative_melt_density)
                                   * melt_fraction(in.temperature[i], in.pressure[i], in.composition[i], in.position[i]);

          // in the end, all the influences are added up
          out.densities[i] = (reference_rho + composition_dependence + pressure_dependence) * temperature_dependence
                             * (1.0 - melt_dependence);


          out.entropy_derivative_pressure[i] = entropy_derivative (in.temperature[i], in.pressure[i], composition,
                                                                   in.position[i], NonlinearDependence::pressure) ; // for pressure dependence

          out.entropy_derivative_temperature[i] = entropy_derivative (in.temperature[i], in.pressure[i], composition,
                                                                      in.position[i], NonlinearDependence::temperature) ; // for temperature dependence

        }
    }


    template <int dim>
    double
    LatentHeatMelt<dim>::
    reference_viscosity () const
    {
      return eta;
    }


    template <int dim>
    double
    LatentHeatMelt<dim>::
    entropy_derivative (const double temperature,
                        const double pressure,
                        const std::vector<double> &compositional_fields,
                        const Point<dim> &position,
                        const NonlinearDependence::Dependence dependence) const
    {
      double entropy_gradient = 0.0;

      // calculate latent heat of melting
      // we need the change of melt fraction in dependence of pressure and temperature

      // for peridotite after Katz, 2003
      const double T_solidus        = A1 + 273.15
                                      + A2 * pressure
                                      + A3 * pressure * pressure;
      const double T_lherz_liquidus = B1 + 273.15
                                      + B2 * pressure
                                      + B3 * pressure * pressure;
      const double T_liquidus       = C1 + 273.15
                                      + C2 * pressure
                                      + C3 * pressure * pressure;

      const double dT_solidus_dp        = A2 + 2 * A3 * pressure;
      const double dT_lherz_liquidus_dp = B2 + 2 * B3 * pressure;
      const double dT_liquidus_dp       = C2 + 2 * C3 * pressure;

      const double peridotite_fraction = (this->n_compositional_fields()>0
                                          ?
                                          1.0 - compositional_fields[0]
                                          :
                                          1.0);

      if (temperature > T_solidus && temperature < T_liquidus && pressure < 1.3e10)
        {
          // melt fraction when clinopyroxene is still present
          double melt_fraction_derivative_temperature
            = beta * pow((temperature - T_solidus)/(T_lherz_liquidus - T_solidus),beta-1)
              / (T_lherz_liquidus - T_solidus);

          double melt_fraction_derivative_pressure
            = beta * pow((temperature - T_solidus)/(T_lherz_liquidus - T_solidus),beta-1)
              * (dT_solidus_dp * (temperature - T_lherz_liquidus)
                 + dT_lherz_liquidus_dp * (T_solidus - temperature))
              / pow(T_lherz_liquidus - T_solidus,2);

          // melt fraction after melting of all clinopyroxene
          const double R_cpx = r1 + r2 * pressure;
          const double F_max = M_cpx / R_cpx;

          if (peridotite_melt_fraction(temperature, pressure, compositional_fields, position) > F_max)
            {
              const double T_max = std::pow(F_max,1.0/beta) * (T_lherz_liquidus - T_solidus) + T_solidus;
              const double dF_max_dp = - M_cpx * std::pow(r1 + r2 * pressure,-2) * r2;
              const double dT_max_dp = dT_solidus_dp
                                       + 1.0/beta * std::pow(F_max,1.0/beta - 1.0) * dF_max_dp * (T_lherz_liquidus - T_solidus)
                                       + std::pow(F_max,1.0/beta) * (dT_lherz_liquidus_dp - dT_solidus_dp);

              melt_fraction_derivative_temperature
                = (1.0 - F_max) * beta * std::pow((temperature - T_max)/(T_liquidus - T_max),beta-1)
                  / (T_liquidus - T_max);

              melt_fraction_derivative_pressure
                = dF_max_dp
                  - dF_max_dp * std::pow((temperature - T_max)/(T_liquidus - T_max),beta)
                  + (1.0 - F_max) * beta * std::pow((temperature - T_max)/(T_liquidus - T_max),beta-1)
                  * (dT_max_dp * (T_max - T_liquidus) - (dT_liquidus_dp - dT_max_dp) * (temperature - T_max)) / std::pow(T_liquidus - T_max, 2);
            }

          double melt_fraction_derivative = 0;
          if (dependence == NonlinearDependence::temperature)
            melt_fraction_derivative = melt_fraction_derivative_temperature;
          else if (dependence == NonlinearDependence::pressure)
            melt_fraction_derivative = melt_fraction_derivative_pressure;
          else
            AssertThrow(false, ExcMessage("Error in calculating melt fraction derivative: not implemented"));

          entropy_gradient += melt_fraction_derivative * peridotite_melting_entropy_change * peridotite_fraction;
        }


      // for melting of pyroxenite after Sobolev et al., 2011
      if (this->n_compositional_fields()>0)
        {
          // calculate change of entropy for melting all material
          const double X = pyroxenite_melt_fraction(temperature, pressure, compositional_fields, position);

          // calculate change of melt fraction in dependence of pressure and temperature
          const double T_melting = D1 + 273.15
                                   + D2 * pressure
                                   + D3 * pressure * pressure;
          const double dT_melting_dp = 2*D3*pressure + D2;
          const double discriminant = E1*E1/(E2*E2*4) + (temperature-T_melting)/E2;

          double melt_fraction_derivative = 0.0;
          if (temperature > T_melting && X < F_px_max && pressure < 1.3e10)
            {
              if (dependence == NonlinearDependence::temperature)
                melt_fraction_derivative = -1.0/(2*E2 * sqrt(discriminant));
              else if (dependence == NonlinearDependence::pressure)
                melt_fraction_derivative = (dT_melting_dp)/(2*E2 * sqrt(discriminant));
              else
                AssertThrow(false, ExcMessage("Error in calculating melt fraction derivative: not implemented"));
            }

          entropy_gradient += melt_fraction_derivative * pyroxenite_melting_entropy_change * compositional_fields[0];
        }

      return entropy_gradient;
    }


    template <int dim>
    void
    LatentHeatMelt<dim>::
    melt_fractions (const MaterialModel::MaterialModelInputs<dim> &in,
                    std::vector<double> &melt_fractions) const
    {
      for (unsigned int q=0; q<in.n_evaluation_points(); ++q)
        melt_fractions[q] = melt_fraction(in.temperature[q],
                                          std::max(0.0, in.pressure[q]),
                                          in.composition[q],
                                          in.position[q]);
    }

    template <int dim>
    double
    LatentHeatMelt<dim>::
    peridotite_melt_fraction (const double temperature,
                              const double pressure,
                              const std::vector<double> &,
                              const Point<dim> &) const
    {
      // anhydrous melting of peridotite after Katz, 2003
      const double T_solidus  = A1 + 273.15
                                + A2 * pressure
                                + A3 * pressure * pressure;
      const double T_lherz_liquidus = B1 + 273.15
                                      + B2 * pressure
                                      + B3 * pressure * pressure;
      const double T_liquidus = C1 + 273.15
                                + C2 * pressure
                                + C3 * pressure * pressure;

      // melt fraction for peridotite with clinopyroxene
      double peridotite_melt_fraction;
      if (temperature < T_solidus || pressure > 1.3e10)
        peridotite_melt_fraction = 0.0;
      else if (temperature > T_lherz_liquidus)
        peridotite_melt_fraction = 1.0;
      else
        peridotite_melt_fraction = std::pow((temperature - T_solidus) / (T_lherz_liquidus - T_solidus),beta);

      // melt fraction after melting of all clinopyroxene
      const double R_cpx = r1 + r2 * pressure;
      const double F_max = M_cpx / R_cpx;

      if (peridotite_melt_fraction > F_max && temperature < T_liquidus)
        {
          const double T_max = std::pow(F_max,1/beta) * (T_lherz_liquidus - T_solidus) + T_solidus;
          peridotite_melt_fraction = F_max + (1 - F_max) * pow((temperature - T_max) / (T_liquidus - T_max),beta);
        }
      return peridotite_melt_fraction;

    }

    template <int dim>
    double
    LatentHeatMelt<dim>::
    pyroxenite_melt_fraction (const double temperature,
                              const double pressure,
                              const std::vector<double> &,
                              const Point<dim> &) const
    {
      // melting of pyroxenite after Sobolev et al., 2011
      const double T_melting = D1 + 273.15
                               + D2 * pressure
                               + D3 * pressure * pressure;

      const double discriminant = E1*E1/(E2*E2*4) + (temperature-T_melting)/E2;

      double pyroxenite_melt_fraction;
      if (temperature < T_melting || pressure > 1.3e10)
        pyroxenite_melt_fraction = 0.0;
      else if (discriminant < 0)
        pyroxenite_melt_fraction = F_px_max;
      else
        pyroxenite_melt_fraction = -E1/(2*E2) - std::sqrt(discriminant);

      return pyroxenite_melt_fraction;
    }

    template <int dim>
    double
    LatentHeatMelt<dim>::
    melt_fraction (const double temperature,
                   const double pressure,
                   const std::vector<double> &composition, /*composition*/
                   const Point<dim> &position) const
    {
      return (this->n_compositional_fields()>0
              ?
              pyroxenite_melt_fraction(temperature, pressure, composition, position)
              * composition[0]
              +
              peridotite_melt_fraction(temperature, pressure, composition, position)
              * (1.0 - composition[0])
              :
              peridotite_melt_fraction(temperature, pressure, composition, position));

    }



    template <int dim>
    bool
    LatentHeatMelt<dim>::
    is_compressible () const
    {
      if (reference_compressibility > 0)
        return true;
      else
        return false;
    }



    template <int dim>
    void
    LatentHeatMelt<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Latent heat melt");
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
                             "The value of the thermal expansion coefficient $\\alpha_s$. "
                             "Units: \\si{\\per\\kelvin}.");
          prm.declare_entry ("Thermal expansion coefficient of melt", "6.8e-5",
                             Patterns::Double (0.),
                             "The value of the thermal expansion coefficient $\\alpha_f$. "
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
          prm.declare_entry ("A1", "1085.7",
                             Patterns::Double (),
                             "Constant parameter in the quadratic "
                             "function that approximates the solidus "
                             "of peridotite. "
                             "Units: \\si{\\degreeCelsius}.");
          prm.declare_entry ("A2", "1.329e-7",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the quadratic function that approximates "
                             "the solidus of peridotite. "
                             "Units: \\si{\\degreeCelsius\\per\\pascal}.");
          prm.declare_entry ("A3", "-5.1e-18",
                             Patterns::Double (),
                             "Prefactor of the quadratic pressure term "
                             "in the quadratic function that approximates "
                             "the solidus of peridotite. "
                             "Units: \\si{\\degreeCelsius\\per\\pascal\\squared}.");
          prm.declare_entry ("B1", "1475.0",
                             Patterns::Double (),
                             "Constant parameter in the quadratic "
                             "function that approximates the lherzolite "
                             "liquidus used for calculating the fraction "
                             "of peridotite-derived melt. "
                             "Units: \\si{\\degreeCelsius}.");
          prm.declare_entry ("B2", "8.0e-8",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the quadratic function that approximates "
                             "the  lherzolite liquidus used for "
                             "calculating the fraction of peridotite-"
                             "derived melt. "
                             "Units: \\si{\\degreeCelsius\\per\\pascal}.");
          prm.declare_entry ("B3", "-3.2e-18",
                             Patterns::Double (),
                             "Prefactor of the quadratic pressure term "
                             "in the quadratic function that approximates "
                             "the  lherzolite liquidus used for "
                             "calculating the fraction of peridotite-"
                             "derived melt. "
                             "Units: \\si{\\degreeCelsius\\per\\pascal\\squared}.");
          prm.declare_entry ("C1", "1780.0",
                             Patterns::Double (),
                             "Constant parameter in the quadratic "
                             "function that approximates the liquidus "
                             "of peridotite. "
                             "Units: \\si{\\degreeCelsius}.");
          prm.declare_entry ("C2", "4.50e-8",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the quadratic function that approximates "
                             "the liquidus of peridotite. "
                             "Units: \\si{\\degreeCelsius\\per\\pascal}.");
          prm.declare_entry ("C3", "-2.0e-18",
                             Patterns::Double (),
                             "Prefactor of the quadratic pressure term "
                             "in the quadratic function that approximates "
                             "the liquidus of peridotite. "
                             "Units: \\si{\\degreeCelsius\\per\\pascal\\squared}.");
          prm.declare_entry ("r1", "0.5",
                             Patterns::Double (),
                             "Constant in the linear function that "
                             "approximates the clinopyroxene reaction "
                             "coefficient. "
                             "Units: non-dimensional.");
          prm.declare_entry ("r2", "8e-11",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the linear function that approximates "
                             "the clinopyroxene reaction coefficient. "
                             "Units: \\si{\\per\\pascal}.");
          prm.declare_entry ("beta", "1.5",
                             Patterns::Double (),
                             "Exponent of the melting temperature in "
                             "the melt fraction calculation. "
                             "Units: non-dimensional.");
          prm.declare_entry ("Peridotite melting entropy change", "-300.",
                             Patterns::Double (),
                             "The entropy change for the phase transition "
                             "from solid to melt of peridotite. "
                             "Units: \\si{\\joule\\per\\kelvin\\per\\kilogram}.");
          prm.declare_entry ("Mass fraction cpx", "0.15",
                             Patterns::Double (),
                             "Mass fraction of clinopyroxene in the "
                             "peridotite to be molten. "
                             "Units: non-dimensional.");
          prm.declare_entry ("D1", "976.0",
                             Patterns::Double (),
                             "Constant parameter in the quadratic "
                             "function that approximates the solidus "
                             "of pyroxenite. "
                             "Units: \\si{\\degreeCelsius}.");
          prm.declare_entry ("D2", "1.329e-7",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the quadratic function that approximates "
                             "the solidus of pyroxenite. "
                             "Note that this factor is different from the "
                             "value given in Sobolev, 2011, because they use "
                             "the potential temperature whereas we use the "
                             "absolute temperature. "
                             "Units: \\si{\\degreeCelsius\\per\\pascal}.");
          prm.declare_entry ("D3", "-5.1e-18",
                             Patterns::Double (),
                             "Prefactor of the quadratic pressure term "
                             "in the quadratic function that approximates "
                             "the solidus of pyroxenite. "
                             "Units: \\si{\\degreeCelsius\\per\\pascal\\squared}.");
          prm.declare_entry ("E1", "663.8",
                             Patterns::Double (),
                             "Prefactor of the linear depletion term "
                             "in the quadratic function that approximates "
                             "the melt fraction of pyroxenite. "
                             "Units: \\si{\\degreeCelsius\\per\\pascal}.");
          prm.declare_entry ("E2", "-611.4",
                             Patterns::Double (),
                             "Prefactor of the quadratic depletion term "
                             "in the quadratic function that approximates "
                             "the melt fraction of pyroxenite. "
                             "Units: \\si{\\degreeCelsius\\per\\pascal\\squared}.");
          prm.declare_entry ("Pyroxenite melting entropy change", "-400.",
                             Patterns::Double (),
                             "The entropy change for the phase transition "
                             "from solid to melt of pyroxenite. "
                             "Units: \\si{\\joule\\per\\kelvin\\per\\kilogram}.");
          prm.declare_entry ("Maximum pyroxenite melt fraction", "0.5429",
                             Patterns::Double (),
                             "Maximum melt fraction of pyroxenite "
                             "in this parameterization. At higher "
                             "temperatures peridotite begins to melt.");
          prm.declare_entry ("Relative density of melt", "0.9",
                             Patterns::Double (),
                             "The relative density of melt compared to the "
                             "solid material. This means, the density change "
                             "upon melting is this parameter times the density "
                             "of solid material."
                             "Units: non-dimensional.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    LatentHeatMelt<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Latent heat melt");
        {
          reference_rho              = prm.get_double ("Reference density");
          reference_T                = prm.get_double ("Reference temperature");
          eta                        = prm.get_double ("Viscosity");
          composition_viscosity_prefactor = prm.get_double ("Composition viscosity prefactor");
          thermal_viscosity_exponent = prm.get_double ("Thermal viscosity exponent");
          k_value                    = prm.get_double ("Thermal conductivity");
          reference_specific_heat    = prm.get_double ("Reference specific heat");
          thermal_alpha              = prm.get_double ("Thermal expansion coefficient");
          melt_thermal_alpha         = prm.get_double ("Thermal expansion coefficient of melt");
          reference_compressibility  = prm.get_double ("Compressibility");
          compositional_delta_rho    = prm.get_double ("Density differential for compositional field 1");

          if (thermal_viscosity_exponent!=0.0 && reference_T == 0.0)
            AssertThrow(false, ExcMessage("Error: Material model latent heat with Thermal viscosity exponent can not have reference_T=0."));

          A1              = prm.get_double ("A1");
          A2              = prm.get_double ("A2");
          A3              = prm.get_double ("A3");
          B1              = prm.get_double ("B1");
          B2              = prm.get_double ("B2");
          B3              = prm.get_double ("B3");
          C1              = prm.get_double ("C1");
          C2              = prm.get_double ("C2");
          C3              = prm.get_double ("C3");
          r1              = prm.get_double ("r1");
          r2              = prm.get_double ("r2");
          beta            = prm.get_double ("beta");
          peridotite_melting_entropy_change
            = prm.get_double ("Peridotite melting entropy change");

          M_cpx           = prm.get_double ("Mass fraction cpx");
          D1              = prm.get_double ("D1");
          D2              = prm.get_double ("D2");
          D3              = prm.get_double ("D3");
          E1              = prm.get_double ("E1");
          E2              = prm.get_double ("E2");
          pyroxenite_melting_entropy_change
            = prm.get_double ("Pyroxenite melting entropy change");

          F_px_max        = prm.get_double ("Maximum pyroxenite melt fraction");
          relative_melt_density = prm.get_double ("Relative density of melt");
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
    ASPECT_REGISTER_MATERIAL_MODEL(LatentHeatMelt,
                                   "latent heat melt",
                                   "A material model that includes the latent heat of melting "
                                   "for two materials: peridotite and pyroxenite. The melting model "
                                   "for peridotite is taken from Katz et al., 2003 (A new "
                                   "parameterization of hydrous mantle melting) and the one for "
                                   "pyroxenite from Sobolev et al., 2011 (Linking mantle plumes, "
                                   "large igneous provinces and environmental catastrophes). "
                                   "The model assumes a constant entropy change for melting 100\\% "
                                   "of the material, which can be specified in the input file. "
                                   "The partial derivatives of entropy with respect to temperature "
                                   "and pressure required for calculating the latent heat consumption "
                                   "are then calculated as product of this constant entropy change, "
                                   "and the respective derivative of the function the describes the "
                                   "melt fraction. This is linearly averaged with respect to the "
                                   "fractions of the two materials present. "
                                   "If no compositional fields are specified in the input file, the "
                                   "model assumes that the material is peridotite. If compositional "
                                   "fields are specified, the model assumes that the first compositional "
                                   "field is the fraction of pyroxenite and the rest of the material "
                                   "is peridotite. "
                                   "\n\n"
                                   "Otherwise, this material model has a temperature- and pressure-"
                                   "dependent density and viscosity and the density and thermal "
                                   "expansivity depend on the melt fraction present. "
                                   "It is possible to extent this model to include a melt fraction "
                                   "dependence of all the material parameters by calling the "
                                   "function melt_fraction in the calculation of the respective "
                                   "parameter. "
                                   "However, melt and solid move with the same velocity and "
                                   "melt extraction is not taken into account (batch melting). ")
  }
}
