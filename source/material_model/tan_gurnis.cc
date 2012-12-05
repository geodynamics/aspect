/*
  Copyright (C) 2011, 2012 by the authors of the ASPECT code.

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
/*  $Id$  */


#include <aspect/material_model/tan_gurnis.h>
#include <deal.II/base/parameter_handler.h>

using namespace dealii;

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    TanGurnis<dim>::TanGurnis()
    {
      a=0; // 0 or 2

      //BA:
      //Di=0;gamma=10000; //=inf
      //EBA:
      //Di=0.5;gamma=inf;
      //TALA:
      Di=0.5;
      gamma=1.0;

      wavenumber=1;
    }


    template <int dim>
    double
    TanGurnis<dim>::
    viscosity (const double,
               const double,
               const std::vector<double> &,       /*composition*/
               const SymmetricTensor<2,dim> &,
               const Point<dim> &pos) const
    {
      const double depth = 1.0-pos(dim-1);
      return exp(a*depth);
    }


    template <int dim>
    double
    TanGurnis<dim>::
    reference_viscosity () const
    {
      return 1.0;
    }



    template <int dim>
    double
    TanGurnis<dim>::
    reference_density () const
    {
      return 1.0;
    }



    template <int dim>
    double
    TanGurnis<dim>::
    reference_thermal_expansion_coefficient () const
    {
      return 1.0;
    }



    template <int dim>
    double
    TanGurnis<dim>::
    specific_heat (const double,
                   const double,
                   const std::vector<double> &, /*composition*/
                   const Point<dim> &) const
    {
      return 1250;
    }



    template <int dim>
    double
    TanGurnis<dim>::
    reference_cp () const
    {
      return 1250;
    }



    template <int dim>
    double
    TanGurnis<dim>::
    thermal_conductivity (const double,
                          const double,
                          const std::vector<double> &, /*composition*/
                          const Point<dim> &) const
    {
      return 2e-5;
    }



    template <int dim>
    double
    TanGurnis<dim>::
    reference_thermal_diffusivity () const
    {
      return k_value/(reference_rho*reference_specific_heat);
    }



    template <int dim>
    double
    TanGurnis<dim>::
    density (const double,
             const double,
             const std::vector<double> &, /*composition*/
             const Point<dim> &pos) const
    {
      const double depth = 1.0-pos(dim-1);
      const double temperature = sin(numbers::PI*pos(dim-1))*cos(numbers::PI*wavenumber*pos(0));
      return -1.0*temperature*exp(Di/gamma*(depth));
    }



    template <int dim>
    double
    TanGurnis<dim>::
    thermal_expansion_coefficient (const double,
                                   const double,
                                   const std::vector<double> &, /*composition*/
                                   const Point<dim> &) const
    {
      return thermal_alpha;
    }



    template <int dim>
    double
    TanGurnis<dim>::
    compressibility (const double temperature,
                     const double pressure,
                     const std::vector<double> &compositional_fields,
                     const Point<dim> &pos) const
    {
      return Di/gamma / density(temperature, pressure, compositional_fields, pos);
    }



    template <int dim>
    bool
    TanGurnis<dim>::
    viscosity_depends_on (const NonlinearDependence::Dependence) const
    {
      return false;
    }



    template <int dim>
    bool
    TanGurnis<dim>::
    density_depends_on (const NonlinearDependence::Dependence) const
    {
      return false;
    }



    template <int dim>
    bool
    TanGurnis<dim>::
    compressibility_depends_on (const NonlinearDependence::Dependence) const
    {
      return false;
    }



    template <int dim>
    bool
    TanGurnis<dim>::
    specific_heat_depends_on (const NonlinearDependence::Dependence) const
    {
      return false;
    }



    template <int dim>
    bool
    TanGurnis<dim>::
    thermal_conductivity_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      return false;
    }



    template <int dim>
    bool
    TanGurnis<dim>::
    is_compressible () const
    {
      return true;
    }



    template <int dim>
    double
    TanGurnis<dim>::
    parameter_a() const
    {
      return a;
    }



    template <int dim>
    double
    TanGurnis<dim>::
    parameter_wavenumber() const
    {
      return wavenumber;
    }



    template <int dim>
    double
    TanGurnis<dim>::
    parameter_Di() const
    {
      return Di;
    }



    template <int dim>
    double
    TanGurnis<dim>::
    parameter_gamma() const
    {
      return gamma;
    }



    template <int dim>
    void
    TanGurnis<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Tan Gurnis model");
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
          prm.declare_entry ("Thermal conductivity", "4.7",
                             Patterns::Double (0),
                             "The value of the thermal conductivity $k$. "
                             "Units: $W/m/K$.");
          prm.declare_entry ("Reference specific heat", "1250",
                             Patterns::Double (0),
                             "The value of the specific heat $cp$. "
                             "Units: $J/kg/K$.");
          prm.declare_entry ("Thermal expansion coefficient", "2e-5",
                             Patterns::Double (0),
                             "The value of the thermal expansion coefficient $\\beta$. "
                             "Units: $1/K$.");
          prm.declare_entry ("a", "0",
                             Patterns::Double (0),
                             "");
          prm.declare_entry ("Di", "0.5",
                             Patterns::Double (0),
                             "");
          prm.declare_entry ("gamma", "1",
                             Patterns::Double (0),
                             "");
          prm.declare_entry ("wavenumber", "1",
                             Patterns::Double (0),
                             "");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    TanGurnis<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Tan Gurnis model");
        {
          reference_rho     = prm.get_double ("Reference density");
          reference_T = prm.get_double ("Reference temperature");
          eta                   = prm.get_double ("Viscosity");
          k_value               = prm.get_double ("Thermal conductivity");
          reference_specific_heat = prm.get_double ("Reference specific heat");
          thermal_alpha = prm.get_double ("Thermal expansion coefficient");
          a = prm.get_double("a");
          Di = prm.get_double("Di");
          gamma = prm.get_double("gamma");
          wavenumber = prm.get_double("wavenumber");
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
    ASPECT_REGISTER_MATERIAL_MODEL(TanGurnis,
                                   "Tan Gurnis",
                                   "A simple compressible material model based on a benchmark"
                                   " from the paper of Tan/Gurnis (2007). This does not use the"
                                   " temperature equation, but has a hardcoded temperature.")
  }
}
