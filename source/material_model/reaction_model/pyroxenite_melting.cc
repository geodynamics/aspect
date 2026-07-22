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



#include <aspect/material_model/reaction_model/pyroxenite_melting.h>
#include <aspect/utilities.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <deal.II/base/parameter_handler.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace ReactionModel
    {
      // melting of pyroxenite after Sobolev et al., 2011. 
      // Originally implemented in melt fraction postprocessor.
      template <int dim>
      double
      PyroxeniteMelting<dim>::
      melt_fraction (const double temperature,
                     const double pressure) const
      {
        const double T_melting = D1 + 273.15
                                  + D2 * pressure
                                  + D3 * pressure * pressure;

        const double discriminant = E1*E1/(E2*E2*4) + (temperature-T_melting)/E2;

        double pyroxenite_melt_fraction;
        if (temperature < T_melting || pressure > 1.3e10)
          pyroxenite_melt_fraction = 0.0;
        else if (discriminant < 0)
          pyroxenite_melt_fraction = 0.5429;
        else
          pyroxenite_melt_fraction = -E1/(2*E2) - std::sqrt(discriminant);

        return pyroxenite_melt_fraction;
      }


      template <int dim>
      void
      PyroxeniteMelting<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("D1", "976.0",
                            Patterns::Double (),
                            "Constant parameter in the quadratic "
                            "function that approximates the solidus "
                            "of pyroxenite. "
                            "Units: $^\\circ\\text{C}$.");
        prm.declare_entry ("D2", "1.329e-7",
                            Patterns::Double (),
                            "Prefactor of the linear pressure term "
                            "in the quadratic function that approximates "
                            "the solidus of pyroxenite. "
                            "Note that this factor is different from the "
                            "value given in Sobolev, 2011, because they use "
                            "the potential temperature whereas we use the "
                            "absolute temperature. "
                            "$\\frac{^\\circ\\text{C}}{\\text{Pa}}$.");
        prm.declare_entry ("D3", "-5.1e-18",
                            Patterns::Double (),
                            "Prefactor of the quadratic pressure term "
                            "in the quadratic function that approximates "
                            "the solidus of pyroxenite. "
                            "$\\frac{^\\circ\\text{C}}{\\text{Pa}^2}$.");
        prm.declare_entry ("E1", "663.8",
                            Patterns::Double (),
                            "Prefactor of the linear depletion term "
                            "in the quadratic function that approximates "
                            "the melt fraction of pyroxenite. "
                            "$\\frac{^\\circ\\text{C}}{\\text{Pa}}$.");
        prm.declare_entry ("E2", "-611.4",
                            Patterns::Double (),
                            "Prefactor of the quadratic depletion term "
                            "in the quadratic function that approximates "
                            "the melt fraction of pyroxenite. "
                            "$\\frac{^\\circ\\text{C}}{\\text{Pa}^2}$.");


      }


      template <int dim>
      void
      PyroxeniteMelting<dim>::parse_parameters (ParameterHandler &prm)
      {
        D1              = prm.get_double ("D1");
        D2              = prm.get_double ("D2");
        D3              = prm.get_double ("D3");
        E1              = prm.get_double ("E1");
        E2              = prm.get_double ("E2");
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    namespace ReactionModel
    {
#define INSTANTIATE(dim) \
  template class PyroxeniteMelting<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)
#undef INSTANTIATE
    }
  }
}
