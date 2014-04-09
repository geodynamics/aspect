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


#include <aspect/material_model/multicomponent.h>
#include <deal.II/base/parameter_handler.h>

using namespace dealii;

namespace aspect
{
  namespace MaterialModel
  {

    template <int dim>
    double
    Multicomponent<dim>::
    viscosity (const double temperature,
               const double,
               const std::vector<double> &composition,       /*composition*/
               const SymmetricTensor<2,dim> &,
               const Point<dim> &p) const
    {
      double visc = eta;
      if ((composition_viscosity_prefactor.empty() == false) && (composition.size() > 0))
        {
          double sum_comp = 0;
          for(unsigned int i=0; i< composition.size(); ++i)
            sum_comp += std::max(0.0, composition[i]);


         
          switch (ViscosityAveraging)
          {
            case Arithmetic:
              {
                double arithmetic_mean = (1.0 - sum_comp);
                for(unsigned int i=0; i< composition.size(); ++i)
                  arithmetic_mean += composition[i]*composition_viscosity_prefactor[i];
                visc = arithmetic_mean * eta;
              }
              break;

            case Harmonic:
              {
                double harmonic_mean = (1.0-sum_comp)/eta;
                for(unsigned int i=0; i< composition.size(); ++i)
                  harmonic_mean += composition[i]/(composition_viscosity_prefactor[i]*eta);
                visc = 1.0/harmonic_mean;
              }
              break;
 
            case Geometric:
              //TODO Make this make sense
              //double geometric_mean = 0.0;
              //for(unsigned int i=0; i < composition.size(); ++i)
              //  geometric_mean += composition[i]*std::log(composition_viscosity_prefactor[i]);
              //visc = std::exp(geometric_mean);
              break;
 
            case MaximumComposition:
              {
                double max_composition_visc = 1.0;
                double max_comp = std::max(0.0,1.0-sum_comp);
                for(unsigned int i=0; i< composition.size(); ++i)
                  if (composition[i] > max_comp) {max_composition_visc = composition_viscosity_prefactor[i]; max_comp=composition[i];}
                visc = eta*max_composition_visc;
              }
          }
      }
      AssertThrow(visc>0.0, ExcInternalError());

      return visc;
 
    }


    template <int dim>
    double
    Multicomponent<dim>::
    reference_viscosity () const
    {
      return eta;
    }

    template <int dim>
    double
    Multicomponent<dim>::
    reference_density () const
    {
      return reference_rho;
    }

    template <int dim>
    double
    Multicomponent<dim>::
    reference_thermal_expansion_coefficient () const
    {
      return thermal_alpha;
    }

    template <int dim>
    double
    Multicomponent<dim>::
    specific_heat (const double,
                   const double,
                   const std::vector<double> &, /*composition*/
                   const Point<dim> &) const
    {
      return reference_specific_heat;
    }

    template <int dim>
    double
    Multicomponent<dim>::
    reference_cp () const
    {
      return reference_specific_heat;
    }

    template <int dim>
    double
    Multicomponent<dim>::
    thermal_conductivity (const double,
                          const double,
                          const std::vector<double> &, /*composition*/
                          const Point<dim> &) const
    {
      return k_value;
    }

    template <int dim>
    double
    Multicomponent<dim>::
    reference_thermal_diffusivity () const
    {
      return k_value/(reference_rho*reference_specific_heat);
    }

    template <int dim>
    double
    Multicomponent<dim>::
    density (const double temperature,
             const double,
             const std::vector<double> &compositional_fields, /*composition*/
             const Point<dim> &) const
    {
      double temperature_factor= (1.0 - thermal_alpha * (temperature - reference_T));

      double composition_difference = 0;
      for(unsigned int i=0; i< compositional_fields.size(); ++i)
        composition_difference += compositional_fields[i]*composition_delta_rho[i];

      AssertThrow((reference_rho + composition_difference)*temperature_factor>0.0, ExcInternalError());

      return (reference_rho + composition_difference)*temperature_factor;
    }


    template <int dim>
    double
    Multicomponent<dim>::
    thermal_expansion_coefficient (const double temperature,
                                   const double,
                                   const std::vector<double> &, /*composition*/
                                   const Point<dim> &) const
    {
      return thermal_alpha;
    }


    template <int dim>
    double
    Multicomponent<dim>::
    compressibility (const double,
                     const double,
                     const std::vector<double> &, /*composition*/
                     const Point<dim> &) const
    {
      return 0.0;
    }

    template <int dim>
    bool
    Multicomponent<dim>::
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
               (composition_viscosity_prefactor.empty() == false))
        return true;
      else
        return false;
    }


    template <int dim>
    bool
    Multicomponent<dim>::
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
               (composition_delta_rho.empty() == false))
        return true;
      else
        return false;
    }

    template <int dim>
    bool
    Multicomponent<dim>::
    compressibility_depends_on (const NonlinearDependence::Dependence) const
    {
      return false;
    }

    template <int dim>
    bool
    Multicomponent<dim>::
    specific_heat_depends_on (const NonlinearDependence::Dependence) const
    {
      return false;
    }

    template <int dim>
    bool
    Multicomponent<dim>::
    thermal_conductivity_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      return false;
    }


    template <int dim>
    bool
    Multicomponent<dim>::
    is_compressible () const
    {
      return false;
    }



    template <int dim>
    void
    Multicomponent<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Multicomponent");
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
          prm.declare_entry ("Density differential for compositional fields", "",
                             Patterns::List(Patterns::Double()),
                             "If compositional fields are used, then one would frequently want "
                             "change in compositional densities.  Expects a comma separated list "
                             " of density differences [kg/m^3] with the order corresponding to  "
                             "the ordering of the compositional fields");
          prm.declare_entry ("Compositional viscosity prefactors", "",
                             Patterns::List(Patterns::Double(0)),
                             "Dimensionless prefactor for the viscosities of the compositional fields"
                             " Expects a comma separated list with prefactors that multiply the reference"
                             " viscosity, and with an order corresponding to  "
                             "the ordering of the compositional fields");
          prm.declare_entry("Viscosity averaging scheme", "Harmonic",
                             Patterns::Selection("Arithmetic|Harmonic|Geometric|Maximum composition"),
                             "When more than one compositional field is present at a point "
                             "with different viscosities, we need to come up with an average "
                             "viscosity at that point.  Select a weighted harmonic, arithmetic, "
                             "or geometric average");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    Multicomponent<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Multicomponent");
        {
          reference_rho              = prm.get_double ("Reference density");
          reference_T                = prm.get_double ("Reference temperature");
          eta                        = prm.get_double ("Viscosity");
          thermal_viscosity_exponent = prm.get_double ("Thermal viscosity exponent");
          k_value                    = prm.get_double ("Thermal conductivity");
          reference_specific_heat    = prm.get_double ("Reference specific heat");
          thermal_alpha              = prm.get_double ("Thermal expansion coefficient");
 
          if (prm.get ("Viscosity averaging scheme") == "Harmonic")
            ViscosityAveraging = Harmonic; 
          else if (prm.get ("Viscosity averaging scheme") == "Arithmetic")
            ViscosityAveraging = Arithmetic; 
          else if (prm.get ("Viscosity averaging scheme") == "Geometric")
            ViscosityAveraging = Geometric; 
          else if (prm.get ("Viscosity averaging scheme") == "Maximum composition")
            ViscosityAveraging = MaximumComposition; 

          if (thermal_viscosity_exponent!=0.0 && reference_T == 0.0)
            AssertThrow(false, ExcMessage("Error: Material model simple with Thermal viscosity exponent can not have reference_T=0."));
              
          composition_delta_rho = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Density differential for compositional fields")));
          composition_viscosity_prefactor = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Compositional viscosity prefactors")));

          if (composition_delta_rho.size() != composition_viscosity_prefactor.size())
            AssertThrow(false, ExcMessage("Error: Must have same number of depths and viscosity prefactors"));
          //TODO: Make sure that the list of compositional density differences and viscosity prefactors
          //are the same as the number of compositional fields
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
    ASPECT_REGISTER_MATERIAL_MODEL(Multicomponent,
                                   "multicomponent",
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
                                   "'Multicomponent model'.")
  }
}
