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


#include <aspect/material_model/reaction_model/katz2003_mantle_melting.h>
#include <aspect/utilities.h>
#include <deal.II/base/parameter_handler.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace ReactionModel
    {
      template <int dim>
      Katz2003MantleMelting<dim>::Katz2003MantleMelting()
        = default;

      template <int dim>
      double
      Katz2003MantleMelting<dim>::
      melt_fraction (const double temperature,
                     const double pressure) const
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
        const double R_cpx = r1 + r2 * std::max(0.0, pressure);
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
      Katz2003MantleMelting<dim>::
      entropy_change (const double temperature,
                      const double pressure,
                      const double maximum_melt_fraction,
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
            const double R_cpx = r1 + r2 * std::max(0.0, pressure);
            const double F_max = M_cpx / R_cpx;

            if (melt_fraction(temperature, pressure) > F_max)
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
              AssertThrow(false, ExcMessage("not implemented"));

            if (melt_fraction(temperature, pressure) >= maximum_melt_fraction)
              entropy_gradient = melt_fraction_derivative * peridotite_melting_entropy_change;
          }
        return entropy_gradient;
      }

      template <int dim>
      void
      Katz2003MantleMelting<dim>::declare_parameters (ParameterHandler &prm)
      {
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
        prm.declare_entry ("Mass fraction cpx", "0.15",
                           Patterns::Double (),
                           "Mass fraction of clinopyroxene in the "
                           "peridotite to be molten. "
                           "Units: non-dimensional.");
        prm.declare_entry ("Peridotite melting entropy change", "-300.",
                           Patterns::Double (),
                           "The entropy change for the phase transition "
                           "from solid to melt of peridotite. "
                           "Units: \\si{\\joule\\per\\kelvin\\per\\kilogram}.");
      }

      template <int dim>
      void
      Katz2003MantleMelting<dim>::parse_parameters (ParameterHandler &prm)
      {
        A1    = prm.get_double ("A1");
        A2    = prm.get_double ("A2");
        A3    = prm.get_double ("A3");
        B1    = prm.get_double ("B1");
        B2    = prm.get_double ("B2");
        B3    = prm.get_double ("B3");
        C1    = prm.get_double ("C1");
        C2    = prm.get_double ("C2");
        C3    = prm.get_double ("C3");
        r1    = prm.get_double ("r1");
        r2    = prm.get_double ("r2");
        beta  = prm.get_double ("beta");
        M_cpx = prm.get_double ("Mass fraction cpx");
        peridotite_melting_entropy_change
          = prm.get_double ("Peridotite melting entropy change");
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
  namespace ReactionModel \
  { \
    template class Katz2003MantleMelting<dim>; \
  }

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
