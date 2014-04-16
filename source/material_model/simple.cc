/*
  Copyright (C) 2011 - 2014 by the authors of the ASPECT code.

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


#include <aspect/material_model/simple.h>
#include <deal.II/base/parameter_handler.h>

using namespace dealii;

namespace aspect
{
  namespace MaterialModel
  {

    template <int dim>
    double
    Simple<dim>::
    viscosity (const double temperature,
               const double,
               const std::vector<double> &composition,       /*composition*/
               const SymmetricTensor<2,dim> &,
               const Point<dim> &) const
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
    Simple<dim>::
    reference_viscosity () const
    {
      return eta;
    }

    template <int dim>
    double
    Simple<dim>::
    reference_density () const
    {
      return reference_rho;
    }

    template <int dim>
    double
    Simple<dim>::
    reference_thermal_expansion_coefficient () const
    {
      return thermal_alpha;
    }

    template <int dim>
    double
    Simple<dim>::
    specific_heat (const double,
                   const double,
                   const std::vector<double> &, /*composition*/
                   const Point<dim> &) const
    {
      return reference_specific_heat;
    }

    template <int dim>
    double
    Simple<dim>::
    reference_cp () const
    {
      return reference_specific_heat;
    }

    template <int dim>
    double
    Simple<dim>::
    thermal_conductivity (const double,
                          const double,
                          const std::vector<double> &, /*composition*/
                          const Point<dim> &) const
    {
      return k_value;
    }

    template <int dim>
    double
    Simple<dim>::
    reference_thermal_diffusivity () const
    {
      return k_value/(reference_rho*reference_specific_heat);
    }

    template <int dim>
    double
    Simple<dim>::
    density (const double temperature,
             const double,
             const std::vector<double> &compositional_fields, /*composition*/
             const Point<dim> &) const
    {
      const double c = compositional_fields.size()>0?
                       std::max(0.0, compositional_fields[0])
                       :
                       0.0;
      return reference_rho * (1 - thermal_alpha * (temperature - reference_T))
             + compositional_delta_rho * c;
    }


    template <int dim>
    double
    Simple<dim>::
    thermal_expansion_coefficient (const double temperature,
                                   const double,
                                   const std::vector<double> &, /*composition*/
                                   const Point<dim> &) const
    {
      return thermal_alpha;
    }


    template <int dim>
    double
    Simple<dim>::
    compressibility (const double,
                     const double,
                     const std::vector<double> &, /*composition*/
                     const Point<dim> &) const
    {
      return 0.0;
    }

    template <int dim>
    bool
    Simple<dim>::
    viscosity_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      // compare this with the implementation of the viscosity() function
      // to see the dependencies
      if (((dependence & NonlinearDependence::temperature) != NonlinearDependence::none)
          &&
          (thermal_viscosity_exponent != 0))
        return true;
      else if (((dependence & NonlinearDependence::compositional_fields) != NonlinearDependence::none)
               &&
               (composition_viscosity_prefactor != 1.0))
        return true;
      else
        return false;
    }


    template <int dim>
    bool
    Simple<dim>::
    density_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      // compare this with the implementation of the density() function
      // to see the dependencies
      if (((dependence & NonlinearDependence::temperature) != NonlinearDependence::none)
          &&
          (thermal_alpha != 0))
        return true;
      else if (((dependence & NonlinearDependence::compositional_fields) != NonlinearDependence::none)
               &&
               (compositional_delta_rho != 0))
        return true;
      else
        return false;
    }

    template <int dim>
    bool
    Simple<dim>::
    compressibility_depends_on (const NonlinearDependence::Dependence) const
    {
      return false;
    }

    template <int dim>
    bool
    Simple<dim>::
    specific_heat_depends_on (const NonlinearDependence::Dependence) const
    {
      return false;
    }

    template <int dim>
    bool
    Simple<dim>::
    thermal_conductivity_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      return false;
    }


    template <int dim>
    bool
    Simple<dim>::
    is_compressible () const
    {
      return false;
    }



    template <int dim>
    void
    Simple<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Simple model");
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
                             "A linear dependency of viscosity on the first compositional field. "
                             "Dimensionless prefactor. With a value of 1.0 (the default) the "
                             "viscosity does not depend on the composition.");
          prm.declare_entry ("Thermal viscosity exponent", "0.0",
                             Patterns::Double (0),
                             "The temperature dependence of viscosity. Dimensionless exponent.");
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
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    Simple<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Simple model");
        {
          reference_rho              = prm.get_double ("Reference density");
          reference_T                = prm.get_double ("Reference temperature");
          eta                        = prm.get_double ("Viscosity");
          composition_viscosity_prefactor = prm.get_double ("Composition viscosity prefactor");
          thermal_viscosity_exponent = prm.get_double ("Thermal viscosity exponent");
          k_value                    = prm.get_double ("Thermal conductivity");
          reference_specific_heat    = prm.get_double ("Reference specific heat");
          thermal_alpha              = prm.get_double ("Thermal expansion coefficient");
          compositional_delta_rho    = prm.get_double ("Density differential for compositional field 1");

          if (thermal_viscosity_exponent!=0.0 && reference_T == 0.0)
            AssertThrow(false, ExcMessage("Error: Material model simple with Thermal viscosity exponent can not have reference_T=0."));
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
    ASPECT_REGISTER_MATERIAL_MODEL(Simple,
                                   "simple",
                                   "A simple material model that has constant values "
                                   "for all coefficients but the density and viscosity. "
                                   "This model uses the formulation that assumes an incompressible"
                                   " medium despite the fact that the density follows the law "
                                   "$\\rho(T)=\\rho_0(1-\\beta(T-T_{\\text{ref}})$. "
                                   "The temperature dependency of viscosity is "
                                   " switched off by default and follows the formula"
                                   "$\\eta(T)=\\eta_0*e^{\\eta_T*\\Delta T / T_{\\text{ref}})}$."
                                   "The value for the components of this formula and additional "
                                   "parameters are read from the parameter file in subsection "
                                   "'Simple model'.")
  }
}
