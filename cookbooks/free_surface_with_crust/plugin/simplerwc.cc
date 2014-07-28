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


#include "./simplerwc.h"
#include <deal.II/base/parameter_handler.h>
#include <iostream>
#include <cmath>

using namespace dealii;

namespace aspect
{
  namespace MaterialModel
  {



    
    
    template <int dim>
    bool
    Simplerwc<dim>::
    viscosity_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      return false;
    }


    template <int dim>
    bool
    Simplerwc<dim>::
    density_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      return false;
    }

    template <int dim>
    bool
    Simplerwc<dim>::
    compressibility_depends_on (const NonlinearDependence::Dependence) const
    {
      return false;
    }

    template <int dim>
    bool
    Simplerwc<dim>::
    specific_heat_depends_on (const NonlinearDependence::Dependence) const
    {
      return false;
    }

    template <int dim>
    bool
    Simplerwc<dim>::
    thermal_conductivity_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      return false;
    }


    template <int dim>
    bool
    Simplerwc<dim>::
    is_compressible () const
    {
      return false;
    }

    template <int dim>
    double
    Simplerwc<dim>::
    reference_viscosity () const
    {
      return eta_L;
    }
 
        

    template <int dim>
    double
    Simplerwc<dim>::
    reference_density () const
    {
      return reference_rho;
    }

    template <int dim>
    void
    Simplerwc<dim>::
    evaluate(const typename Interface<dim>::MaterialModelInputs &in, typename Interface<dim>::MaterialModelOutputs &out ) const
    {
      
      double z;
      for (unsigned int i=0; i<in.position.size(); ++i)
        { 
          z = in.position[i][1];
          out.viscosities[i] = (z>jump_height) ?  eta_U: eta_L ;
                     
          out.densities[i] = reference_rho * (1.0 - thermal_alpha * (in.temperature[i] - reference_T));
          out.thermal_expansion_coefficients[i] = thermal_alpha;
          out.specific_heat[i] = reference_specific_heat;
          out.thermal_conductivities[i] = k_value;
          out.compressibilities[i] = 0.0;
        }
    }


    template <int dim>
    void
   Simplerwc<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Simplerwc model");
        {
          prm.declare_entry ("Reference density", "3300",
                             Patterns::Double (0),
                             "Reference density $\\rho_0$. Units: $kg/m^3$.");
          prm.declare_entry ("Reference temperature", "293",
                             Patterns::Double (0),
                             "The reference temperature $T_0$. The reference temperature is used "
                             "in the density formula. Units: $K$.");
          prm.declare_entry ("Lower viscosity", "5e24",
                             Patterns::Double (0),
                             "The value of the viscosity $\\eta$L. Units: $kg/m/s$.");
          prm.declare_entry ("Upper viscosity", "5e24",
                             Patterns::Double (0),
	                     "The value of the viscosity in the top section $\\eta$U. Units: $kg/m/s$.");
          prm.declare_entry ("Jump height", "100000",
                              Patterns::Double (0),
                              "The height at which the viscosity changes. Units: m.");
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

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    Simplerwc<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Simplerwc model");
        {
          reference_rho              = prm.get_double ("Reference density");
          reference_T                = prm.get_double ("Reference temperature");
          eta_L                      = prm.get_double ("Lower viscosity");
          eta_U                      = prm.get_double ("Upper viscosity");
          jump_height                 = prm.get_double ("Jump height"); 
          k_value                    = prm.get_double ("Thermal conductivity");
          reference_specific_heat    = prm.get_double ("Reference specific heat");
          thermal_alpha              = prm.get_double ("Thermal expansion coefficient");
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
    ASPECT_REGISTER_MATERIAL_MODEL(Simplerwc,
                                   "simplerwc",
                                   "A material model that has constant values "
                                   "except for density, which depends linearly on temperature: "
                                   "\\begin{align}"
                                   "  \\rho(p,T) &= \\left(1-\\alpha (T-T_0)\\right)\\rho_0."
                                   "\\end{align}"
                                   "\n\n"
                                   "\\note{This material model fills the role the ``simple'' material "
                                   "model was originally intended to fill, before the latter acquired "
                                   "all sorts of complicated temperature and compositional dependencies.}")
  }
}
