/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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


#include <aspect/material_model/composition_reaction.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/material_model/equation_of_state/interface.h>


namespace aspect
{
  namespace MaterialModel
  {

    template <int dim>
    void
    CompositionReaction<dim>::
    evaluate(const MaterialModelInputs<dim> &in,
             MaterialModelOutputs<dim> &out) const
    {
      ReactionRateOutputs<dim> *reaction_rate_out = out.template get_additional_output<ReactionRateOutputs<dim>>();

      // The Composition reaction model has up to two compositional fields (plus one background field)
      // that can influence the density
      const unsigned int n_compositions_for_eos = std::min(this->n_compositional_fields()+1, 3u);
      EquationOfStateOutputs<dim> eos_outputs (n_compositions_for_eos);

      for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
        {
          const double temperature = in.temperature[i];
          const std::vector<double> &composition = in.composition[i];
          const double delta_temp = temperature-reference_T;
          double temperature_dependence = std::max(std::min(std::exp(-thermal_viscosity_exponent*delta_temp/reference_T),1e2),1e-2);

          if (std::isnan(temperature_dependence))
            temperature_dependence = 1.0;

          switch (composition.size())
            {
              case 0:
                out.viscosities[i] = temperature_dependence * eta;
                break;
              case 1:
                // geometric interpolation
                out.viscosities[i] = (std::pow(10, ((1-composition[0]) * std::log10(eta*temperature_dependence)
                                                    + composition[0] * std::log10(eta*composition_viscosity_prefactor_1*temperature_dependence))));
                break;
              default:
                out.viscosities[i] = (std::pow(10, ((1 - 0.5*composition[0] - 0.5*composition[1]) * std::log10(eta*temperature_dependence)
                                                    + 0.5 * composition[0] * std::log10(eta*composition_viscosity_prefactor_1*temperature_dependence)
                                                    + 0.5 * composition[1] * std::log10(eta*composition_viscosity_prefactor_2*temperature_dependence))));
                break;
            }

          equation_of_state.evaluate(in, i, eos_outputs);

          std::vector<double> volume_fractions (n_compositions_for_eos, 1.0);
          for (unsigned int c=0; c<n_compositions_for_eos-1; ++c)
            {
              volume_fractions[c+1] = std::max(0.0, in.composition[i][c]);
              volume_fractions[0] -= volume_fractions[c+1];
            }

          out.densities[i] = MaterialUtilities::average_value(volume_fractions, eos_outputs.densities, MaterialUtilities::arithmetic);


          const double depth = this->get_geometry_model().depth(in.position[i]);
          for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
            {
              double delta_C = 0.0;
              switch (c)
                {
                  case 0:
                    if (depth < reaction_depth) delta_C = -composition[0];
                    break;
                  case 1:
                    if (depth < reaction_depth) delta_C = composition[0];
                    break;
                  default:
                    delta_C = 0.0;
                    break;
                }
              out.reaction_terms[i][c] = delta_C;

              // Fill reaction rate outputs instead of the reaction terms if we use operator splitting
              // (and then set the latter to zero).
              if (this->get_parameters().use_operator_splitting)
                {
                  if (reaction_rate_out != nullptr)
                    reaction_rate_out->reaction_rates[i][c] = (this->get_timestep_number() > 0
                                                               ?
                                                               out.reaction_terms[i][c] / this->get_timestep()
                                                               :
                                                               0.0);
                  out.reaction_terms[i][c] = 0.0;
                }
            }

          out.thermal_expansion_coefficients[i] = eos_outputs.thermal_expansion_coefficients[0];
          out.specific_heat[i] = eos_outputs.specific_heat_capacities[0];
          out.thermal_conductivities[i] = k_value;
          out.compressibilities[i] = eos_outputs.compressibilities[0];
        }
    }



    template <int dim>
    bool
    CompositionReaction<dim>::
    is_compressible () const
    {
      return false;
    }



    template <int dim>
    void
    CompositionReaction<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Composition reaction model");
        {
          EquationOfState::LinearizedIncompressible<dim>::declare_parameters (prm, 2);

          prm.declare_entry ("Reference temperature", "293.",
                             Patterns::Double (0.),
                             "The reference temperature $T_0$. Units: \\si{\\kelvin}.");
          prm.declare_entry ("Viscosity", "5e24",
                             Patterns::Double (0.),
                             "The value of the constant viscosity. Units: \\si{\\kilogram\\per\\meter\\per\\second}.");
          prm.declare_entry ("Composition viscosity prefactor 1", "1.0",
                             Patterns::Double (0.),
                             "A linear dependency of viscosity on the first compositional field. "
                             "Dimensionless prefactor. With a value of 1.0 (the default) the "
                             "viscosity does not depend on the composition.");
          prm.declare_entry ("Composition viscosity prefactor 2", "1.0",
                             Patterns::Double (0.),
                             "A linear dependency of viscosity on the second compositional field. "
                             "Dimensionless prefactor. With a value of 1.0 (the default) the "
                             "viscosity does not depend on the composition.");
          prm.declare_entry ("Thermal viscosity exponent", "0.0",
                             Patterns::Double (0.),
                             "The temperature dependence of viscosity. Dimensionless exponent.");
          prm.declare_entry ("Thermal conductivity", "4.7",
                             Patterns::Double (0.),
                             "The value of the thermal conductivity $k$. "
                             "Units: \\si{\\watt\\per\\meter\\per\\kelvin}.");
          prm.declare_entry ("Reaction depth", "0.",
                             Patterns::Double (0.),
                             "Above this depth the compositional fields react: "
                             "The first field gets converted to the second field. "
                             "Units: \\si{\\meter}.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    CompositionReaction<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Composition reaction model");
        {
          equation_of_state.parse_parameters (prm, 2);

          reference_T                = prm.get_double ("Reference temperature");
          eta                        = prm.get_double ("Viscosity");
          composition_viscosity_prefactor_1 = prm.get_double ("Composition viscosity prefactor 1");
          composition_viscosity_prefactor_2 = prm.get_double ("Composition viscosity prefactor 2");
          thermal_viscosity_exponent = prm.get_double ("Thermal viscosity exponent");
          k_value                    = prm.get_double ("Thermal conductivity");
          reaction_depth             = prm.get_double ("Reaction depth");

          if (thermal_viscosity_exponent!=0.0 && reference_T == 0.0)
            AssertThrow(false, ExcMessage("Error: Material model composition reaction with Thermal viscosity exponent can not have reference_T=0."));
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Declare dependencies on solution variables
      this->model_dependence.viscosity = NonlinearDependence::none;
      this->model_dependence.density = NonlinearDependence::temperature | NonlinearDependence::compositional_fields;
      this->model_dependence.compressibility = NonlinearDependence::none;
      this->model_dependence.specific_heat = NonlinearDependence::none;
      this->model_dependence.thermal_conductivity = NonlinearDependence::none;

      if (thermal_viscosity_exponent != 0)
        this->model_dependence.viscosity |= NonlinearDependence::temperature;
      if ((composition_viscosity_prefactor_1 != 1.0) ||
          (composition_viscosity_prefactor_2 != 1.0))
        this->model_dependence.viscosity |= NonlinearDependence::compositional_fields;
    }


    template <int dim>
    void
    CompositionReaction<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      if (this->get_parameters().use_operator_splitting
          && out.template get_additional_output<ReactionRateOutputs<dim>>() == nullptr)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std::make_unique<MaterialModel::ReactionRateOutputs<dim>> (n_points,
                                                                        this->n_compositional_fields()));
        }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(CompositionReaction,
                                   "composition reaction",
                                   "A material model that behaves in the same way as "
                                   "the simple material model, but includes two compositional "
                                   "fields and a reaction between them. Above a depth given "
                                   "in the input file, the first fields gets converted to the "
                                   "second field. ")
  }
}
