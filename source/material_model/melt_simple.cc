/*
  Copyright (C) 2015 - 2024 by the authors of the ASPECT code.

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


#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/material_model/melt_simple.h>
#include <aspect/material_model/reaction_model/katz2003_mantle_melting.h>
#include <aspect/utilities.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/numerics/fe_field_function.h>


namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    double
    MeltSimple<dim>::
    reference_darcy_coefficient () const
    {
      // 0.01 = 1% melt
      return katz2003_model.reference_darcy_coefficient();
    }

    template <int dim>
    bool
    MeltSimple<dim>::
    is_compressible () const
    {
      return model_is_compressible;
    }

    template <int dim>
    void
    MeltSimple<dim>::
    melt_fractions (const MaterialModel::MaterialModelInputs<dim> &in,
                    std::vector<double> &melt_fractions) const
    {
      for (unsigned int q=0; q<in.n_evaluation_points(); ++q)
        melt_fractions[q] = katz2003_model.melt_fraction(in.temperature[q],
                                                         this->get_adiabatic_conditions().pressure(in.position[q]));
    }


    template <int dim>
    void
    MeltSimple<dim>::initialize ()
    {
      if (this->include_melt_transport())
        {
          AssertThrow(this->get_parameters().use_operator_splitting,
                      ExcMessage("The material model ``Melt simple'' can only be used with operator splitting!"));
          AssertThrow(this->introspection().compositional_name_exists("peridotite"),
                      ExcMessage("Material model Melt simple only works if there is a "
                                 "compositional field called peridotite."));
          AssertThrow(this->introspection().compositional_name_exists("porosity"),
                      ExcMessage("Material model Melt simple with melt transport only "
                                 "works if there is a compositional field called porosity."));
        }
    }


    template <int dim>
    void
    MeltSimple<dim>::
    evaluate(const typename Interface<dim>::MaterialModelInputs &in, typename Interface<dim>::MaterialModelOutputs &out) const
    {
      for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
        {
          // calculate density first, we need it for the reaction term
          // first, calculate temperature dependence of density
          double temperature_dependence = 1.0;
          if (this->include_adiabatic_heating ())
            {
              // temperature dependence is 1 - alpha * (T - T(adiabatic))
              temperature_dependence -= (in.temperature[i] - this->get_adiabatic_conditions().temperature(in.position[i]))
                                        * thermal_expansivity;
            }
          else
            temperature_dependence -= (in.temperature[i] - reference_T) * thermal_expansivity;

          // calculate composition dependence of density
          const double delta_rho = this->introspection().compositional_name_exists("peridotite")
                                   ?
                                   depletion_density_change * in.composition[i][this->introspection().compositional_index_for_name("peridotite")]
                                   :
                                   0.0;
          out.densities[i] = (reference_rho_solid + delta_rho)
                             * temperature_dependence * std::exp(compressibility * (in.pressure[i] - this->get_surface_pressure()));

          out.viscosities[i] = eta_0;
          out.thermal_expansion_coefficients[i] = thermal_expansivity;
          out.specific_heat[i] = reference_specific_heat;
          out.thermal_conductivities[i] = thermal_conductivity;
          out.compressibilities[i] = compressibility;

          double visc_temperature_dependence = 1.0;
          if (this->include_adiabatic_heating ())
            {
              const double delta_temp = in.temperature[i]-this->get_adiabatic_conditions().temperature(in.position[i]);
              visc_temperature_dependence = std::max(std::min(std::exp(-thermal_viscosity_exponent*delta_temp/this->get_adiabatic_conditions().temperature(in.position[i])),1e4),1e-4);
            }
          else
            {
              const double delta_temp = in.temperature[i]-reference_T;
              const double T_dependence = (thermal_viscosity_exponent == 0.0
                                           ?
                                           0.0
                                           :
                                           thermal_viscosity_exponent*delta_temp/reference_T);
              visc_temperature_dependence = std::max(std::min(std::exp(-T_dependence),1e4),1e-4);
            }
          out.viscosities[i] *= visc_temperature_dependence;

        }

      katz2003_model.calculate_reaction_rate_outputs(in, out);
      katz2003_model.calculate_fluid_outputs(in, out, reference_T);
    }


    template <int dim>
    void
    MeltSimple<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Melt simple");
        {
          // Melt Fraction Parameters
          ReactionModel::Katz2003MantleMelting<dim>::declare_parameters(prm);


          prm.declare_entry ("Use full compressibility", "false",
                             Patterns::Bool (),
                             "If the compressibility should be used everywhere in the code "
                             "(if true), changing the volume of material when the density changes, "
                             "or only in the momentum conservation and advection equations "
                             "(if false).");
          prm.declare_entry ("Thermal expansion coefficient", "2e-5",
                             Patterns::Double (0.),
                             "The value of the thermal expansion coefficient $\\beta$. "
                             "Units: \\si{\\per\\kelvin}.");
          prm.declare_entry ("Reference shear viscosity", "5e20",
                             Patterns::Double (0.),
                             "The value of the constant viscosity $\\eta_0$ of the solid matrix. "
                             "This viscosity may be modified by both temperature and porosity "
                             "dependencies. Units: \\si{\\pascal\\second}.");
          prm.declare_entry ("Reference specific heat", "1250.",
                             Patterns::Double (0.),
                             "The value of the specific heat $C_p$. "
                             "Units: \\si{\\joule\\per\\kelvin\\per\\kilogram}.");
          prm.declare_entry ("Thermal conductivity", "4.7",
                             Patterns::Double (0.),
                             "The value of the thermal conductivity $k$. "
                             "Units: \\si{\\watt\\per\\meter\\per\\kelvin}.");
          prm.declare_entry ("Solid compressibility", "0.0",
                             Patterns::Double (0.),
                             "The value of the compressibility of the solid matrix. "
                             "Units: \\si{\\per\\pascal}.");
          prm.declare_entry ("Thermal viscosity exponent", "0.0",
                             Patterns::Double (0.),
                             "The temperature dependence of the shear viscosity. Dimensionless exponent. "
                             "See the general documentation "
                             "of this model for a formula that states the dependence of the "
                             "viscosity on this factor, which is called $\\beta$ there.");
          prm.declare_entry ("Reference temperature", "293.",
                             Patterns::Double (0.),
                             "The reference temperature $T_0$. The reference temperature is used "
                             "in both the density and viscosity formulas. Units: \\si{\\kelvin}.");
          prm.declare_entry ("Depletion density change", "0.0",
                             Patterns::Double (),
                             "The density contrast between material with a depletion of 1 and a "
                             "depletion of zero. Negative values indicate lower densities of "
                             "depleted material. Depletion is indicated by the compositional "
                             "field with the name peridotite. Not used if this field does not "
                             "exist in the model. "
                             "Units: \\si{\\kilogram\\per\\meter\\cubed}.");
          prm.declare_entry ("Reference solid density", "3000.",
                             Patterns::Double (0.),
                             "Reference density of the solid $\\rho_{s,0}$. "
                             "Units: \\si{\\kilogram\\per\\meter\\cubed}.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    MeltSimple<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Melt simple");
        {
          model_is_compressible      = prm.get_bool ("Use full compressibility");
          reference_specific_heat    = prm.get_double ("Reference specific heat");
          eta_0                      = prm.get_double ("Reference shear viscosity");
          thermal_expansivity        = prm.get_double ("Thermal expansion coefficient");
          thermal_conductivity       = prm.get_double ("Thermal conductivity");
          compressibility            = prm.get_double ("Solid compressibility");
          thermal_viscosity_exponent = prm.get_double ("Thermal viscosity exponent");
          reference_T                = prm.get_double ("Reference temperature");
          depletion_density_change   = prm.get_double ("Depletion density change");
          reference_rho_solid        = prm.get_double ("Reference solid density");



          if (thermal_viscosity_exponent!=0.0 && reference_T == 0.0)
            AssertThrow(false, ExcMessage("Error: Material model Melt simple with Thermal viscosity exponent can not have reference_T=0."));

          // Melt model
          katz2003_model.initialize_simulator (this->get_simulator());
          katz2003_model.parse_parameters(prm);

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    MeltSimple<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      if (this->get_parameters().use_operator_splitting && out.template get_additional_output<ReactionRateOutputs<dim>>() == nullptr)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std::make_unique<MaterialModel::ReactionRateOutputs<dim>> (n_points, this->n_compositional_fields()));
        }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(MeltSimple,
                                   "melt simple",
                                   "A material model that implements a simple formulation of the "
                                   "material parameters required for the modeling of melt transport, "
                                   "including a source term for the porosity according to the melting "
                                   "model for dry peridotite of \\cite{katz:etal:2003}. This also includes a "
                                   "computation of the latent heat of melting (if the `latent heat' "
                                   "heating model is active)."
                                   "\n\n"
                                   "Most of the material properties are constant, except for the shear, "
                                   "viscosity $\\eta$, the compaction viscosity $\\xi$, and the "
                                   "permeability $k$, which depend on the porosity; and the solid and melt "
                                   "densities, which depend on temperature and pressure:\n "
                                   "$\\eta(\\phi,T) = \\eta_0 e^{\\alpha(\\phi-\\phi_0)} e^{-\\beta(T-T_0)/T_0}$, "
                                   "$\\xi(\\phi,T) = \\xi_0 \\frac{\\phi_0}{\\phi} e^{-\\beta(T-T_0)/T_0}$, "
                                   "$k=k_0 \\phi^n (1-\\phi)^m$, "
                                   "$\\rho=\\rho_0 (1 - \\alpha (T - T_{\\text{adi}})) e^{\\kappa p}$."
                                   "\n\n"
                                   "The model is compressible only if this is specified in the input file, "
                                   "and contains compressibility for both solid and melt.")
  }
}
