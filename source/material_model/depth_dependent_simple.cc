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


#include <aspect/material_model/depth_dependent_simple.h>
#include <deal.II/base/parameter_handler.h>

using namespace dealii;

namespace aspect
{
  namespace MaterialModel
  {

    template <int dim>
    double
    DepthDependentSimple<dim>::
    viscosity (const double temperature,
               const double,
               const std::vector<double> &composition,
               const SymmetricTensor<2,dim> &,
               const Point<dim> &position) const
    {
      const double delta_temp = temperature-reference_T;
      const double temperature_dependence = (reference_T > 0
                                             ?
                                             std::max(std::min(std::exp(-thermal_viscosity_exponent*delta_temp/reference_T),
                                                               1e2),
                                                      1e-2)
                                             :
                                             1.0);

      double depth_viscosity_prefactor = 1.0;
      const double depth = this->get_geometry_model().depth(position);
      Point<1> dpoint(depth);
      depth_viscosity_prefactor = depth_dependence.value(dpoint);
      Assert( depth_viscosity_prefactor > 0.0, ExcMessage("Depth viscosity prefactor must be greater than zero") );

      double composition_dependence = 1.0;
      if ((composition_viscosity_prefactor != 1.0) && (composition.size() > 0))
        {
          //geometric interpolation
          return (pow(10, ((1-composition[0]) * log10(eta*depth_viscosity_prefactor*temperature_dependence)
                           + composition[0] * log10(eta*depth_viscosity_prefactor*composition_viscosity_prefactor*temperature_dependence))));
        }

      return depth_viscosity_prefactor * composition_dependence * temperature_dependence * eta;
    }


    template <int dim>
    double
    DepthDependentSimple<dim>::
    reference_viscosity () const
    {
      return eta;
    }

    template <int dim>
    double
    DepthDependentSimple<dim>::
    reference_density () const
    {
      return reference_rho;
    }

    template <int dim>
    double
    DepthDependentSimple<dim>::
    reference_thermal_expansion_coefficient () const
    {
      return thermal_alpha;
    }

    template <int dim>
    double
    DepthDependentSimple<dim>::
    specific_heat (const double,
                   const double,
                   const std::vector<double> &, /*composition*/
                   const Point<dim> &) const
    {
      return reference_specific_heat;
    }

    template <int dim>
    double
    DepthDependentSimple<dim>::
    reference_cp () const
    {
      return reference_specific_heat;
    }

    template <int dim>
    double
    DepthDependentSimple<dim>::
    thermal_conductivity (const double,
                          const double,
                          const std::vector<double> &, /*composition*/
                          const Point<dim> &) const
    {
      return k_value;
    }

    template <int dim>
    double
    DepthDependentSimple<dim>::
    reference_thermal_diffusivity () const
    {
      return k_value/(reference_rho*reference_specific_heat);
    }

    template <int dim>
    double
    DepthDependentSimple<dim>::
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
    DepthDependentSimple<dim>::
    thermal_expansion_coefficient (const double,
                                   const double,
                                   const std::vector<double> &, /*composition*/
                                   const Point<dim> &) const
    {
      return thermal_alpha;
    }


    template <int dim>
    double
    DepthDependentSimple<dim>::
    compressibility (const double,
                     const double,
                     const std::vector<double> &, /*composition*/
                     const Point<dim> &) const
    {
      return 0.0;
    }

    template <int dim>
    bool
    DepthDependentSimple<dim>::
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
    DepthDependentSimple<dim>::
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
    DepthDependentSimple<dim>::
    compressibility_depends_on (const NonlinearDependence::Dependence) const
    {
      return false;
    }

    template <int dim>
    bool
    DepthDependentSimple<dim>::
    specific_heat_depends_on (const NonlinearDependence::Dependence) const
    {
      return false;
    }

    template <int dim>
    bool
    DepthDependentSimple<dim>::
    thermal_conductivity_depends_on (const NonlinearDependence::Dependence) const
    {
      return false;
    }


    template <int dim>
    bool
    DepthDependentSimple<dim>::
    is_compressible () const
    {
      return false;
    }



    template <int dim>
    void
    DepthDependentSimple<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Depth dependent simple model");
        {
          prm.declare_entry ("Reference density", "3300",
                             Patterns::Double (0),
                             "Reference density $\\rho_0$. Units: $kg/m^3$.");
          prm.declare_entry ("Reference temperature", "293",
                             Patterns::Double (0),
                             "The reference temperature $T_0$. The reference temperature is used "
                             "in both the density and viscosity formulas. Units: $K$.");
          prm.declare_entry ("Viscosity", "5e24",
                             Patterns::Double (0),
                             "The value of the constant viscosity $\\eta_0$. This viscosity may be "
                             "modified by both temperature and compositional dependencies. Units: $kg/m/s$.");
          prm.declare_entry ("Composition viscosity prefactor", "1.0",
                             Patterns::Double (0),
                             "A linear dependency of viscosity on the first compositional field. "
                             "Dimensionless prefactor. With a value of 1.0 (the default) the "
                             "viscosity does not depend on the composition. See the general documentation "
                             "of this model for a formula that states the dependence of the "
                             "viscosity on this factor, which is called $\\xi$ there.");
          prm.declare_entry ("Thermal viscosity exponent", "0.0",
                             Patterns::Double (0),
                             "The temperature dependence of viscosity. Dimensionless exponent. "
                             "See the general documentation "
                             "of this model for a formula that states the dependence of the "
                             "viscosity on this factor, which is called $\\beta$ there.");
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
                             "to make the density depend on these fields. In this depth dependent simple material "
                             "model, we make the following assumptions: if no compositional fields "
                             "are used in the current simulation, then the density is simply the usual "
                             "one with its linear dependence on the temperature. If there are compositional "
                             "fields, then the density only depends on the first one in such a way that "
                             "the density has an additional term of the kind $+\\Delta \\rho \\; c_1(\\mathbf x)$. "
                             "This parameter describes the value of $\\Delta \\rho$. Units: $kg/m^3/\\textrm{unit "
                             "change in composition}$.");
          prm.enter_subsection("Viscosity depth prefactor");
          {
            Functions::ParsedFunction<1>::declare_parameters(prm,1);
            prm.declare_entry("Function expression","1.0");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    DepthDependentSimple<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Depth dependent simple model");
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
            AssertThrow(false, ExcMessage("Error: Material model Depth dependent simple with Thermal viscosity exponent can not have reference_T=0."));

          prm.enter_subsection("Viscosity depth prefactor");
          {
            try
              {
                depth_dependence.parse_parameters(prm);
              }
            catch (...)
              {
                std::cerr << "FunctionParser failed to parse\n"
                          << "\t'Viscosity depth prefactor.Function'\n"
                          << "with expression\n"
                          <<"\t'" << prm.get("Function expression") << "'";
                throw;
              }
          }
          prm.leave_subsection();
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
    ASPECT_REGISTER_MATERIAL_MODEL(DepthDependentSimple,
                                   "depth dependent simple",
                                   "A material model that includes temperature and depth dependence of viscosity. "
                                   "Note that this material model extends the Simple material model. "
                                   "The Depth dependent simple model allows for temperature-, depth-, "
                                   "and composition-dependent viscosity."
                                   " All of the values that define this model are read "
                                   "from a section ``Material model/Depth dependent simple model'' in the input file, see "
                                   "Section~\\ref{parameters:Material_20model/Depth_20dependent_20simple_20model}."
                                   "\n\n"
                                   "This model uses the following set of equations for the two coefficients that "
                                   "are non-constant: "
                                   "\\begin{align}"
                                   "  \\eta(z,T,\\mathfrak c) &= \\eta_z(a) \\tau(T) \\zeta(\\mathfrak c) \\eta_0, \\\\"
                                   "  \\rho(z,T,\\mathfrak c) &= \\left(1-\\alpha (T-T_0)\\right)\\rho_0 + \\Delta\\rho \\; c_0,"
                                   "\\end{align}"
                                   "where $z$ is depth, $\\eta_z$ is a depth-dependent viscosity prefactor defined "
                                   "by the user (assumed to be uniformly equal to 1.0 if not defined)."
                                   "$c_0$ is the first component of the compositional vector "
                                   "$\\mathfrak c$ if the model uses compositional fields, or zero otherwise. "
                                   "\n\n"
                                   "The temperature pre-factor for the viscosity formula above is "
                                   "defined as "
                                   "\\begin{align}"
                                   "  \\tau(T) &= H\\left(e^{\\beta (T-T_0)/T_0}\\right),"
                                   "  \\qquad\\qquad H(x) = \\begin{cases}"
                                   "                            10^{-2} & \\text{if}\\; x<10^{-2}, \\\\"
                                   "                            x & \\text{if}\\; 10^{-2}\\le x \\le 10^2, \\\\"
                                   "                            10^{2} & \\text{if}\\; x>10^{2}, \\\\"
                                   "                         \\end{cases}"
                                   "\\end{align} "
                                   "where $\\beta$ corresponds to the input parameter ``Thermal viscosity exponent'' "
                                   "and $T_0$ to the parameter ``Reference temperature''. If you set $T_0=0$ "
                                   "in the input file, the thermal pre-factor $\\tau(T)=1$."
                                   "\n\n"
                                   "The compositional pre-factor for the viscosity is defined as "
                                   "\\begin{align}"
                                   "  \\zeta(\\mathfrak c) &= \\xi^{c_0}"
                                   "\\end{align} "
                                   "if the model has compositional fields and equals one otherwise. $\\xi$ "
                                   "corresponds to the parameter ``Composition viscosity prefactor'' in the "
                                   "input file."
                                   "if the model has depth-dependent viscosity, the depth-dependence is specified"
                                   "in the input file using a ParsedFunction in the subsection to the "
                                   "Depth dependent simple model: ``Viscosity depth prefactor''"
                                   "\n\n"
                                   "Finally, in the formula for the density, $\\Delta\\rho$ "
                                   "corresponds to the parameter ``Density differential for compositional field 1''."
                                   "\n\n"
                                   "Note that this model uses the formulation that assumes an incompressible "
                                   "medium despite the fact that the density follows the law "
                                   "$\\rho(T)=\\rho_0(1-\\beta(T-T_{\\text{ref}}))$. "
                                   "\n\n")
  }
}
