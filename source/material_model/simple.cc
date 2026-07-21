/*
  Copyright (C) 2011 - 2023 by the authors of the ASPECT code.

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


#include <aspect/material_model/simple.h>
#include <aspect/material_model/equation_of_state/interface.h>


namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    void
    Simple<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      // The Simple model has up to one compositional field (plus one background field)
      // that can influence the density
      const unsigned int n_compositions_for_eos = std::min(this->n_compositional_fields()+1, 2u);
      EquationOfStateOutputs<dim> eos_outputs (n_compositions_for_eos);

      for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
        {
          const double delta_temp = in.temperature[i]-reference_T;
          const double temperature_dependence
            = (reference_T > 0
               ?
               std::max(std::min(std::exp(-thermal_viscosity_exponent *
                                          delta_temp/reference_T),
                                 maximum_thermal_prefactor),
                        minimum_thermal_prefactor)
               :
               1.0);

          out.viscosities[i] = ((composition_viscosity_prefactor != 1.0) && (in.composition[i].size()>0))
                               ?
                               // Geometric interpolation
                               std::pow(10.0, ((1-in.composition[i][0]) * std::log10(eta *
                                                                                     temperature_dependence)
                                               + in.composition[i][0] * std::log10(eta *
                                                                                   composition_viscosity_prefactor *
                                                                                   temperature_dependence)))
                               :
                               temperature_dependence * eta;

          equation_of_state.evaluate(in, i, eos_outputs);

          // except for the density, all material properties are constant across compositions
          out.thermal_expansion_coefficients[i] = eos_outputs.thermal_expansion_coefficients[0];
          out.specific_heat[i] = eos_outputs.specific_heat_capacities[0];
          out.thermal_conductivities[i] = k_value;
          out.compressibilities[i] = eos_outputs.compressibilities[0];
          out.entropy_derivative_pressure[i] = eos_outputs.entropy_derivative_pressure[0];
          out.entropy_derivative_temperature[i] = eos_outputs.entropy_derivative_temperature[0];

          // Change in composition due to chemical reactions at the
          // given positions. The term reaction_terms[i][c] is the
          // change in compositional field c at point i.
          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            out.reaction_terms[i][c] = 0.0;

          std::vector<double> volume_fractions (n_compositions_for_eos, 1.0);
          if (in.composition[i].size()>0)
            {
              volume_fractions[1] = std::max(0.0, in.composition[i][0]);
              volume_fractions[0] = 1.0 - volume_fractions[1];
            }

          out.densities[i] = MaterialUtilities::average_value(volume_fractions, eos_outputs.densities, MaterialUtilities::arithmetic);
        }
    }



    template <int dim>
    bool
    Simple<dim>::
    is_compressible () const
    {
      return equation_of_state.is_compressible ();
    }



    template <int dim>
    void
    Simple<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Simple model");
        {
          EquationOfState::LinearizedIncompressible<dim>::declare_parameters (prm, 1);

          prm.declare_entry ("Reference temperature", "293.",
                             Patterns::Double (0.),
                             "The reference temperature $T_0$. The reference temperature is used "
                             "in both the density and viscosity formulas. Units: \\si{\\kelvin}.");
          prm.declare_entry ("Viscosity", "5e24",
                             Patterns::Double (0.),
                             "The value of the constant viscosity $\\eta_0$. This viscosity may be "
                             "modified by both temperature and compositional dependencies. "
                             "Units: \\si{\\pascal\\second}.");
          prm.declare_entry ("Composition viscosity prefactor", "1.0",
                             Patterns::Double (0.),
                             "A linear dependency of viscosity on the first compositional field. "
                             "Dimensionless prefactor. With a value of 1.0 (the default) the "
                             "viscosity does not depend on the composition. See the general documentation "
                             "of this model for a formula that states the dependence of the "
                             "viscosity on this factor, which is called $\\xi$ there.");
          prm.declare_entry ("Thermal viscosity exponent", "0.0",
                             Patterns::Double (0.),
                             "The temperature dependence of viscosity. Dimensionless exponent. "
                             "See the general documentation "
                             "of this model for a formula that states the dependence of the "
                             "viscosity on this factor, which is called $\\beta$ there.");
          prm.declare_entry("Maximum thermal prefactor","1.0e2",
                            Patterns::Double (0.),
                            "The maximum value of the viscosity prefactor associated with temperature "
                            "dependence.");
          prm.declare_entry("Minimum thermal prefactor","1.0e-2",
                            Patterns::Double (0.),
                            "The minimum value of the viscosity prefactor associated with temperature "
                            "dependence.");
          prm.declare_entry ("Thermal conductivity", "4.7",
                             Patterns::Double (0.),
                             "The value of the thermal conductivity $k$. "
                             "Units: \\si{\\watt\\per\\meter\\per\\kelvin}.");
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
          equation_of_state.parse_parameters (prm, 1);

          reference_T                = prm.get_double ("Reference temperature");
          eta                        = prm.get_double ("Viscosity");
          composition_viscosity_prefactor = prm.get_double ("Composition viscosity prefactor");
          thermal_viscosity_exponent = prm.get_double ("Thermal viscosity exponent");
          maximum_thermal_prefactor       = prm.get_double ("Maximum thermal prefactor");
          minimum_thermal_prefactor       = prm.get_double ("Minimum thermal prefactor");
          if ( maximum_thermal_prefactor == 0.0 ) maximum_thermal_prefactor = std::numeric_limits<double>::max();
          if ( minimum_thermal_prefactor == 0.0 ) minimum_thermal_prefactor = std::numeric_limits<double>::min();

          k_value                    = prm.get_double ("Thermal conductivity");

          if (thermal_viscosity_exponent!=0.0 && reference_T == 0.0)
            AssertThrow(false, ExcMessage("Error: Material model simple with Thermal viscosity exponent can not have reference_T=0."));
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Declare dependencies on solution variables
      this->model_dependence.compressibility = NonlinearDependence::none;
      this->model_dependence.specific_heat = NonlinearDependence::none;
      this->model_dependence.thermal_conductivity = NonlinearDependence::none;
      this->model_dependence.viscosity = NonlinearDependence::none;
      this->model_dependence.density = NonlinearDependence::temperature|NonlinearDependence::compositional_fields;

      if (thermal_viscosity_exponent != 0)
        this->model_dependence.viscosity |= NonlinearDependence::temperature;
      if (composition_viscosity_prefactor != 1.0)
        this->model_dependence.viscosity |= NonlinearDependence::compositional_fields;
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
                                   "A material model that has constant values "
                                   "for all coefficients but the density and viscosity. The defaults for all "
                                   "coefficients are chosen to be similar to what is believed to be correct "
                                   "for Earth's mantle. All of the values that define this model are read "
                                   "from a section ``Material model/Simple model'' in the input file, see "
                                   "Section~\\ref{parameters:Material_20model/Simple_20model}."
                                   "\n\n"
                                   "This model uses the following set of equations for the two coefficients that "
                                   "are non-constant: "
                                   "\\begin{align}"
                                   "  \\eta(p,T,\\mathfrak c) &= \\tau(T) \\zeta(\\mathfrak c) \\eta_0, \\\\"
                                   "  \\rho(p,T,\\mathfrak c) &= \\left(1-\\alpha (T-T_0)\\right)\\rho_0 + \\Delta\\rho \\; c_0,"
                                   "\\end{align}"
                                   "where $c_0$ is the first component of the compositional vector "
                                   "$\\mathfrak c$ if the model uses compositional fields, or zero otherwise. "
                                   "\n\n"
                                   "The temperature pre-factor for the viscosity formula above is "
                                   "defined as "
                                   "\\begin{align}"
                                   "  \\tau(T) &= H\\left(e^{-\\beta (T-T_0)/T_0}\\right),"
                                   "\\intertext{with} "
                                   "  \\qquad\\qquad H(x) &= \\begin{cases}"
                                   "                            \\tau_{\\text{min}} & \\text{if}\\; x<\\tau_{\\text{min}}, \\\\"
                                   "                            x & \\text{if}\\; 10^{-2}\\le x \\le 10^2, \\\\"
                                   "                            \\tau_{\\text{max}} & \\text{if}\\; x>\\tau_{\\text{max}}, \\\\"
                                   "                         \\end{cases}"
                                   "\\end{align} "
                                   "where $x=e^{-\\beta (T-T_0)/T_0}$, "
                                   "$\\beta$ corresponds to the input parameter ``Thermal viscosity exponent'', "
                                   "and $T_0$ to the parameter ``Reference temperature''. If you set $T_0=0$ "
                                   "in the input file, the thermal pre-factor $\\tau(T)=1$. The parameters $\\tau_{\\text{min}}$ "
                                   "and $\\tau_{\\text{max}}$ set the minimum and maximum values of the temperature pre-factor "
                                   "and are set using ``Maximum thermal prefactor'' and ``Minimum thermal prefactor''. "
                                   "Specifying a value of 0.0 for the minimum or maximum values will disable pre-factor limiting."
                                   "\n\n"
                                   "The compositional pre-factor for the viscosity is defined as "
                                   "$ \\zeta(\\mathfrak c) = \\xi^{c_0}$ "
                                   "if the model has compositional fields and equals one otherwise. $\\xi$ "
                                   "corresponds to the parameter ``Composition viscosity prefactor'' in the "
                                   "input file."
                                   "\n\n"
                                   "Finally, in the formula for the density, $\\alpha$ corresponds to the "
                                   "``Thermal expansion coefficient'' and "
                                   "$\\Delta\\rho$ "
                                   "corresponds to the parameter ``Density differential for compositional field 1''."
                                   "\n\n"
                                   "Note that this model uses the formulation that assumes an incompressible "
                                   "medium despite the fact that the density follows the law "
                                   "$\\rho(T)=\\rho_0(1-\\alpha(T-T_{\\text{ref}}))$. "
                                   "\n\n"
                                   ":::{note}\n"
                                   "Despite its name, this material model is not exactly ``simple'', "
                                   "as indicated by the formulas above. While it was originally intended "
                                   "to be simple, it has over time acquired all sorts of temperature "
                                   "and compositional dependencies that weren't initially intended. "
                                   "Consequently, there is now a ``simpler'' material model that now fills "
                                   "the role the current model was originally intended to fill.\n"
                                   ":::")
  }
}
