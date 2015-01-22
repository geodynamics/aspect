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



#include <aspect/material_model/simple_compressible.h>
#include <deal.II/base/parameter_handler.h>

using namespace dealii;

namespace aspect
{
  namespace MaterialModel
  {

    template <int dim>
    double
    SimpleCompressible<dim>::
    viscosity (const double,
               const double,
               const std::vector<double> &,
               const SymmetricTensor<2,dim> &,
               const Point<dim> &) const
    {
      return eta;
    }


    template <int dim>
    double
    SimpleCompressible<dim>::
    reference_viscosity () const
    {
      return eta;
    }

    template <int dim>
    double
    SimpleCompressible<dim>::
    reference_density () const
    {
      return reference_rho;
    }

    template <int dim>
    double
    SimpleCompressible<dim>::
    reference_thermal_expansion_coefficient () const
    {
      return thermal_alpha;
    }

    template <int dim>
    double
    SimpleCompressible<dim>::
    specific_heat (const double,
                   const double,
                   const std::vector<double> &, /*composition*/
                   const Point<dim> &) const
    {
      return reference_specific_heat;
    }

    template <int dim>
    double
    SimpleCompressible<dim>::
    reference_cp () const
    {
      return reference_specific_heat;
    }

    template <int dim>
    double
    SimpleCompressible<dim>::
    thermal_conductivity (const double,
                          const double,
                          const std::vector<double> &, /*composition*/
                          const Point<dim> &) const
    {
      return k_value;
    }

    template <int dim>
    double
    SimpleCompressible<dim>::
    reference_thermal_diffusivity () const
    {
      return k_value/(reference_rho*reference_specific_heat);
    }

    template <int dim>
    double
    SimpleCompressible<dim>::
    density (const double temperature,
             const double pressure,
             const std::vector<double> &,
             const Point<dim> &position) const
    {
      double rho = reference_rho * std::exp(reference_compressibility * (pressure - this->get_surface_pressure()));

      if (this->get_adiabatic_conditions().is_initialized())
        rho *= (1 - thermal_alpha * (temperature - this->get_adiabatic_conditions().temperature(position)));

      return rho;
    }


    template <int dim>
    double
    SimpleCompressible<dim>::
    thermal_expansion_coefficient (const double,
                                   const double,
                                   const std::vector<double> &,
                                   const Point<dim> &) const
    {
      return thermal_alpha;
    }


    template <int dim>
    double
    SimpleCompressible<dim>::
    compressibility (const double,
                     const double,
                     const std::vector<double> &,
                     const Point<dim> &) const
    {
      // compressibility = 1/rho drho/dp
      return reference_compressibility;
    }

    template <int dim>
    bool
    SimpleCompressible<dim>::
    viscosity_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      return false;
    }


    template <int dim>
    bool
    SimpleCompressible<dim>::
    density_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      // compare this with the implementation of the density() function
      // to see the dependencies
      if ((dependence & NonlinearDependence::temperature) != NonlinearDependence::none)
        return (thermal_alpha != 0);
      else if ((dependence & NonlinearDependence::pressure) != NonlinearDependence::none)
        return (reference_compressibility != 0);
      else
        return false;
    }

    template <int dim>
    bool
    SimpleCompressible<dim>::
    compressibility_depends_on (const NonlinearDependence::Dependence) const
    {
      return false;
    }

    template <int dim>
    bool
    SimpleCompressible<dim>::
    specific_heat_depends_on (const NonlinearDependence::Dependence) const
    {
      return false;
    }

    template <int dim>
    bool
    SimpleCompressible<dim>::
    thermal_conductivity_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      return false;
    }


    template <int dim>
    bool
    SimpleCompressible<dim>::
    is_compressible () const
    {
      return (reference_compressibility != 0);
    }



    template <int dim>
    void
    SimpleCompressible<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Simple compressible model");
        {
          prm.declare_entry ("Reference density", "3300",
                             Patterns::Double (0),
                             "Reference density $\\rho_0$. Units: $kg/m^3$.");
          prm.declare_entry ("Viscosity", "1e21",
                             Patterns::Double (0),
                             "The value of the constant viscosity $\\eta_0$. Units: $kg/m/s$.");
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
          prm.declare_entry ("Reference compressibility", "4e-12",
                             Patterns::Double (0),
                             "The value of the reference compressibility. "
                             "Units: $1/Pa$.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    SimpleCompressible<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Simple compressible model");
        {
          reference_rho              = prm.get_double ("Reference density");
          eta                        = prm.get_double ("Viscosity");
          k_value                    = prm.get_double ("Thermal conductivity");
          reference_specific_heat    = prm.get_double ("Reference specific heat");
          thermal_alpha              = prm.get_double ("Thermal expansion coefficient");
          reference_compressibility  = prm.get_double ("Reference compressibility");
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
    ASPECT_REGISTER_MATERIAL_MODEL(SimpleCompressible,
                                   "simple compressible",
                                   "A material model that has constant values "
                                   "for all coefficients but the density. The defaults for all "
                                   "coefficients are chosen to be similar to what is believed to be correct "
                                   "for Earth's mantle. All of the values that define this model are read "
                                   "from a section ``Material model/Simple compressible model'' in the input file, see "
                                   "Section~\\ref{parameters:Material_20model/Simple_20compressible_20model}."
                                   "\n\n"
                                   "This model uses the following equations for the density: "
                                   "\\begin{align}"
                                   "  \\rho(p,T) = \\rho_0"
                                   "              \\left(1-\\alpha (T-T_a)\\right) "
                                   "              \\exp{\\beta (P-P_0))}"
                                   "\\end{align}")
  }
}
