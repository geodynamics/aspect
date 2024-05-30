/*
  Copyright (C) 2015 - 2022 by the authors of the ASPECT code.

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


#include <aspect/material_model/melt_model/katz2003_mantle_melting.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/utilities.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/numerics/fe_field_function.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace MeltModel
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
  namespace MeltModel \
  { \
    template class Katz2003MantleMelting<dim>; \
  }

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
