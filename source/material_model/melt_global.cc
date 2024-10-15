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


#include <aspect/material_model/melt_global.h>
#include <aspect/adiabatic_conditions/interface.h>

#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/numerics/fe_field_function.h>



namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    double
    MeltGlobal<dim>::
    reference_darcy_coefficient () const
    {
      // 0.01 = 1% melt
      return reference_permeability * Utilities::fixed_power<3>(0.01) / eta_f;
    }

    template <int dim>
    bool
    MeltGlobal<dim>::
    is_compressible () const
    {
      return false;
    }

    template <int dim>
    double
    MeltGlobal<dim>::
    melt_fraction (const double temperature,
                   const double pressure,
                   const double depletion) const
    {
      const double T_solidus  = surface_solidus
                                + pressure_solidus_change * pressure
                                + std::max(depletion_solidus_change * depletion, -200.0);
      const double T_liquidus = T_solidus + 500.0;

      double melt_fraction;
      if (temperature < T_solidus)
        melt_fraction = 0.0;
      else if (temperature > T_liquidus)
        melt_fraction = 1.0;
      else
        melt_fraction = (temperature - T_solidus) / (T_liquidus - T_solidus);

      return melt_fraction;
    }


    template <int dim>
    void
    MeltGlobal<dim>::
    melt_fractions (const MaterialModel::MaterialModelInputs<dim> &in,
                    std::vector<double> &melt_fractions) const
    {
      double depletion = 0.0;

      for (unsigned int q=0; q<in.n_evaluation_points(); ++q)
        {
          if (this->include_melt_transport())
            {
              const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
              const unsigned int peridotite_idx = this->introspection().compositional_index_for_name("peridotite");
              depletion = in.composition[q][peridotite_idx] - in.composition[q][porosity_idx];
            }
          melt_fractions[q] = this->melt_fraction(in.temperature[q],
                                                  std::max(0.0, in.pressure[q]),
                                                  depletion);
        }
    }


    template <int dim>
    void
    MeltGlobal<dim>::
    evaluate(const typename Interface<dim>::MaterialModelInputs &in, typename Interface<dim>::MaterialModelOutputs &out) const
    {
      std::vector<double> old_porosity(in.n_evaluation_points());

      ReactionRateOutputs<dim> *reaction_rate_out = out.template get_additional_output<ReactionRateOutputs<dim>>();

      // we want to get the porosity field from the old solution here,
      // because we need a field that is not updated in the nonlinear iterations
      if (this->include_melt_transport() && in.current_cell.state() == IteratorState::valid
          && this->get_timestep_number() > 0 && !this->get_parameters().use_operator_splitting)
        {
          // Prepare the field function

          Functions::FEFieldFunction<dim, LinearAlgebra::BlockVector>

          fe_value(this->get_dof_handler(), this->get_old_solution(), this->get_mapping());

          const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");

          fe_value.set_active_cell(in.current_cell);
          fe_value.value_list(in.position,
                              old_porosity,
                              this->introspection().component_indices.compositional_fields[porosity_idx]);
        }
      else if (this->get_parameters().use_operator_splitting)
        for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
          {
            const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
            old_porosity[i] = in.composition[i][porosity_idx];
          }

      for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
        {
          // calculate density first, we need it for the reaction term
          // temperature dependence of density is 1 - alpha * (T - T(adiabatic))
          double temperature_dependence = 1.0;
          if (this->include_adiabatic_heating ())
            temperature_dependence -= (in.temperature[i] - this->get_adiabatic_conditions().temperature(in.position[i]))
                                      * thermal_expansivity;
          else
            temperature_dependence -= (in.temperature[i] - reference_T) * thermal_expansivity;

          // calculate composition dependence of density
          const double delta_rho = this->introspection().compositional_name_exists("peridotite")
                                   ?
                                   depletion_density_change * in.composition[i][this->introspection().compositional_index_for_name("peridotite")]
                                   :
                                   0.0;
          out.densities[i] = (reference_rho_s + delta_rho) * temperature_dependence
                             * std::exp(compressibility * (in.pressure[i] - this->get_surface_pressure()));

          out.viscosities[i] = eta_0;
          // By default, no melting or freezing --> set all reactions to zero
          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            {
              out.reaction_terms[i][c] = 0.0;

              if (this->get_parameters().use_operator_splitting && reaction_rate_out != nullptr)
                reaction_rate_out->reaction_rates[i][c] = 0.0;
            }

          if (this->include_melt_transport())
            {
              const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
              const double porosity = std::min(1.0, std::max(in.composition[i][porosity_idx],0.0));

              // calculate viscosity based on local melt
              out.viscosities[i] *= std::exp(- alpha_phi * porosity);

              if (include_melting_and_freezing && (in.requests_property(MaterialProperties::reaction_terms) ||
                                                   in.requests_property(MaterialProperties::reaction_rates) ||
                                                   in.requests_property(MaterialProperties::viscosity)))
                {
                  Assert(this->get_timestep_number()<=1 || std::isfinite(in.strain_rate[i].norm()),
                         ExcMessage("Invalid strain_rate in the MaterialModelInputs. This is likely because it was "
                                    "not filled by the caller."));
                  const double trace_strain_rate =
                    (this->get_timestep_number() > 1) ?
                    (trace(in.strain_rate[i]))
                    :
                    numbers::signaling_nan<double>();

                  const unsigned int peridotite_idx = this->introspection().compositional_index_for_name("peridotite");

                  // Calculate the melting rate as difference between the equilibrium melt fraction
                  // and the solution of the previous time step (or the current solution, in case
                  // operator splitting is used).
                  // The solidus is lowered by previous melting events (fractional melting).
                  const double eq_melt_fraction = melt_fraction(in.temperature[i],
                                                                this->get_adiabatic_conditions().pressure(in.position[i]),
                                                                in.composition[i][peridotite_idx] - in.composition[i][porosity_idx]);
                  double porosity_change = eq_melt_fraction - old_porosity[i];
                  // do not allow negative porosity
                  if (old_porosity[i] + porosity_change < 0)
                    porosity_change = -old_porosity[i];

                  for (unsigned int c = 0; c < in.composition[i].size(); ++c)
                    {
                      if (c == peridotite_idx && this->get_timestep_number() > 1)
                        out.reaction_terms[i][c] = porosity_change - in.composition[i][peridotite_idx] * trace_strain_rate * this->get_timestep();
                      else if (c == porosity_idx && this->get_timestep_number() > 1)
                        out.reaction_terms[i][c] = porosity_change * out.densities[i] / this->get_timestep();
                      else
                        out.reaction_terms[i][c] = 0.0;

                      // fill reaction rate outputs if the model uses operator splitting
                      if (this->get_parameters().use_operator_splitting)
                        {
                          if (reaction_rate_out != nullptr)
                            {
                              if (c == peridotite_idx && this->get_timestep_number() > 0)
                                reaction_rate_out->reaction_rates[i][c] = porosity_change / melting_time_scale - in.composition[i][peridotite_idx] * trace(in.strain_rate[i]);
                              else if (c == porosity_idx && this->get_timestep_number() > 0)
                                reaction_rate_out->reaction_rates[i][c] = porosity_change / melting_time_scale;
                              else
                                reaction_rate_out->reaction_rates[i][c] = 0.0;
                            }
                          out.reaction_terms[i][c] = 0.0;
                        }
                    }

                  // find depletion = peridotite, which might affect shear viscosity:
                  const double depletion_visc = std::min(1.0, std::max(in.composition[i][peridotite_idx], 0.0));

                  // calculate strengthening due to depletion:
                  const double depletion_strengthening = std::min(std::exp(alpha_depletion * depletion_visc), delta_eta_depletion_max);

                  // calculate viscosity change due to local melt and depletion:
                  out.viscosities[i] *= depletion_strengthening;
                }
            }

          out.entropy_derivative_pressure[i]    = 0.0;
          out.entropy_derivative_temperature[i] = 0.0;
          out.thermal_expansion_coefficients[i] = thermal_expansivity;
          out.specific_heat[i] = reference_specific_heat;
          out.thermal_conductivities[i] = thermal_conductivity;
          out.compressibilities[i] = 0.0;

          double visc_temperature_dependence = 1.0;
          if (this->include_adiabatic_heating ())
            {
              const double delta_temp = in.temperature[i]-this->get_adiabatic_conditions().temperature(in.position[i]);
              visc_temperature_dependence = std::max(std::min(std::exp(-thermal_viscosity_exponent*delta_temp/this->get_adiabatic_conditions().temperature(in.position[i])),1e4),1e-4);
            }
          else if (thermal_viscosity_exponent != 0.0)
            {
              const double delta_temp = in.temperature[i]-reference_T;
              visc_temperature_dependence = std::max(std::min(std::exp(-thermal_viscosity_exponent*delta_temp/reference_T),1e4),1e-4);
            }
          out.viscosities[i] *= visc_temperature_dependence;
        }

      // fill melt outputs if they exist
      MeltOutputs<dim> *melt_out = out.template get_additional_output<MeltOutputs<dim>>();

      if (melt_out != nullptr)
        {
          const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");

          for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
            {
              double porosity = std::max(in.composition[i][porosity_idx],0.0);

              melt_out->fluid_viscosities[i] = eta_f;
              melt_out->permeabilities[i] = reference_permeability * Utilities::fixed_power<3>(porosity) * Utilities::fixed_power<2>(1.0-porosity);
              melt_out->fluid_density_gradients[i] = Tensor<1,dim>();

              // temperature dependence of density is 1 - alpha * (T - T(adiabatic))
              double temperature_dependence = 1.0;
              if (this->include_adiabatic_heating ())
                temperature_dependence -= (in.temperature[i] - this->get_adiabatic_conditions().temperature(in.position[i]))
                                          * thermal_expansivity;
              else
                temperature_dependence -= (in.temperature[i] - reference_T) * thermal_expansivity;
              melt_out->fluid_densities[i] = reference_rho_f * temperature_dependence
                                             * std::exp(melt_compressibility * (in.pressure[i] - this->get_surface_pressure()));

              melt_out->compaction_viscosities[i] = xi_0 * std::exp(- alpha_phi * porosity);

              double visc_temperature_dependence = 1.0;
              if (this->include_adiabatic_heating ())
                {
                  const double delta_temp = in.temperature[i]-this->get_adiabatic_conditions().temperature(in.position[i]);
                  visc_temperature_dependence = std::max(std::min(std::exp(-thermal_bulk_viscosity_exponent*delta_temp/this->get_adiabatic_conditions().temperature(in.position[i])),1e4),1e-4);
                }
              else if (thermal_viscosity_exponent != 0.0)
                {
                  const double delta_temp = in.temperature[i]-reference_T;
                  visc_temperature_dependence = std::max(std::min(std::exp(-thermal_bulk_viscosity_exponent*delta_temp/reference_T),1e4),1e-4);
                }
              melt_out->compaction_viscosities[i] *= visc_temperature_dependence;
            }
        }
    }



    template <int dim>
    void
    MeltGlobal<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Melt global");
        {
          prm.declare_entry ("Reference solid density", "3000.",
                             Patterns::Double (0.),
                             "Reference density of the solid $\\rho_{s,0}$. "
                             "Units: \\si{\\kilogram\\per\\meter\\cubed}.");
          prm.declare_entry ("Reference melt density", "2500.",
                             Patterns::Double (0.),
                             "Reference density of the melt/fluid$\\rho_{f,0}$. "
                             "Units: \\si{\\kilogram\\per\\meter\\cubed}.");
          prm.declare_entry ("Reference temperature", "293.",
                             Patterns::Double (0.),
                             "The reference temperature $T_0$. The reference temperature is used "
                             "in both the density and viscosity formulas. Units: \\si{\\kelvin}.");
          prm.declare_entry ("Reference shear viscosity", "5e20",
                             Patterns::Double (0.),
                             "The value of the constant viscosity $\\eta_0$ of the solid matrix. "
                             "This viscosity may be modified by both temperature and porosity "
                             "dependencies. Units: \\si{\\pascal\\second}.");
          prm.declare_entry ("Reference bulk viscosity", "1e22",
                             Patterns::Double (0.),
                             "The value of the constant bulk viscosity $\\xi_0$ of the solid matrix. "
                             "This viscosity may be modified by both temperature and porosity "
                             "dependencies. Units: \\si{\\pascal\\second}.");
          prm.declare_entry ("Reference melt viscosity", "10.",
                             Patterns::Double (0.),
                             "The value of the constant melt viscosity $\\eta_f$. Units: \\si{\\pascal\\second}.");
          prm.declare_entry ("Exponential melt weakening factor", "27.",
                             Patterns::Double (0.),
                             "The porosity dependence of the viscosity. Units: dimensionless.");
          prm.declare_entry ("Thermal viscosity exponent", "0.0",
                             Patterns::Double (0.),
                             "The temperature dependence of the shear viscosity. Dimensionless exponent. "
                             "See the general documentation "
                             "of this model for a formula that states the dependence of the "
                             "viscosity on this factor, which is called $\\beta$ there.");
          prm.declare_entry ("Thermal bulk viscosity exponent", "0.0",
                             Patterns::Double (0.),
                             "The temperature dependence of the bulk viscosity. Dimensionless exponent. "
                             "See the general documentation "
                             "of this model for a formula that states the dependence of the "
                             "viscosity on this factor, which is called $\\beta$ there.");
          prm.declare_entry ("Thermal conductivity", "4.7",
                             Patterns::Double (0.),
                             "The value of the thermal conductivity $k$. "
                             "Units: \\si{\\watt\\per\\meter\\per\\kelvin}.");
          prm.declare_entry ("Reference specific heat", "1250.",
                             Patterns::Double (0.),
                             "The value of the specific heat $C_p$. "
                             "Units: \\si{\\joule\\per\\kelvin\\per\\kilogram}.");
          prm.declare_entry ("Thermal expansion coefficient", "2e-5",
                             Patterns::Double (0.),
                             "The value of the thermal expansion coefficient $\\beta$. "
                             "Units: \\si{\\per\\kelvin}.");
          prm.declare_entry ("Reference permeability", "1e-8",
                             Patterns::Double(),
                             "Reference permeability of the solid host rock."
                             "Units: \\si{\\meter\\squared}.");
          prm.declare_entry ("Depletion density change", "0.0",
                             Patterns::Double (),
                             "The density contrast between material with a depletion of 1 and a "
                             "depletion of zero. Negative values indicate lower densities of "
                             "depleted material. Depletion is indicated by the compositional "
                             "field with the name peridotite. Not used if this field does not "
                             "exist in the model. "
                             "Units: \\si{\\kilogram\\per\\meter\\cubed}.");
          prm.declare_entry ("Surface solidus", "1300.",
                             Patterns::Double (0.),
                             "Solidus for a pressure of zero. "
                             "Units: \\si{\\kelvin}.");
          prm.declare_entry ("Depletion solidus change", "200.0",
                             Patterns::Double (),
                             "The solidus temperature change for a depletion of 100\\%. For positive "
                             "values, the solidus gets increased for a positive peridotite field "
                             "(depletion) and lowered for a negative peridotite field (enrichment). "
                             "Scaling with depletion is linear. Only active when fractional melting "
                             "is used. "
                             "Units: \\si{\\kelvin}.");
          prm.declare_entry ("Pressure solidus change", "6e-8",
                             Patterns::Double (),
                             "The linear solidus temperature change with pressure. For positive "
                             "values, the solidus gets increased for positive pressures. "
                             "Units: \\si{\\per\\pascal}.");
          prm.declare_entry ("Solid compressibility", "0.0",
                             Patterns::Double (0.),
                             "The value of the compressibility of the solid matrix. "
                             "Units: \\si{\\per\\pascal}.");
          prm.declare_entry ("Melt compressibility", "0.0",
                             Patterns::Double (0.),
                             "The value of the compressibility of the melt. "
                             "Units: \\si{\\per\\pascal}.");
          prm.declare_entry ("Melt bulk modulus derivative", "0.0",
                             Patterns::Double (0.),
                             "The value of the pressure derivative of the melt bulk "
                             "modulus. "
                             "Units: None.");
          prm.declare_entry ("Include melting and freezing", "true",
                             Patterns::Bool (),
                             "Whether to include melting and freezing (according to a simplified "
                             "linear melting approximation in the model (if true), or not (if "
                             "false).");
          prm.declare_entry ("Melting time scale for operator splitting", "1e3",
                             Patterns::Double (0.),
                             "In case the operator splitting scheme is used, the porosity field can not "
                             "be set to a new equilibrium melt fraction instantly, but the model has to "
                             "provide a melting time scale instead. This time scale defines how fast melting "
                             "happens, or more specifically, the parameter defines the time after which "
                             "the deviation of the porosity from the equilibrium melt fraction will be "
                             "reduced to a fraction of $1/e$. So if the melting time scale is small compared "
                             "to the time step size, the reaction will be so fast that the porosity is very "
                             "close to the equilibrium melt fraction after reactions are computed. Conversely, "
                             "if the melting time scale is large compared to the time step size, almost no "
                             "melting and freezing will occur."
                             "\n\n"
                             "Also note that the melting time scale has to be larger than or equal to the reaction "
                             "time step used in the operator splitting scheme, otherwise reactions can not be "
                             "computed. If the model does not use operator splitting, this parameter is not used. "
                             "Units: yr or s, depending on the ``Use years "
                             "in output instead of seconds'' parameter.");
          prm.declare_entry ("Exponential depletion strengthening factor", "0.0",
                             Patterns::Double (0.),
                             "$\\alpha_F$: exponential dependency of viscosity on the depletion "
                             "field $F$ (called peridotite). "
                             "Dimensionless factor. With a value of 0.0 (the default) the "
                             "viscosity does not depend on the depletion. The effective viscosity increase"
                             "due to depletion is defined as $std::exp( \\alpha_F * F)$. "
                             "Rationale: melting dehydrates the source rock by removing most of the volatiles,"
                             "and makes it stronger. Hirth and Kohlstedt (1996) report typical values around a "
                             "factor 100 to 1000 viscosity contrast between wet and dry rocks, although some "
                             "experimental studies report a smaller (factor 10) contrast (e.g. Fei et al., 2013).");
          prm.declare_entry ("Maximum Depletion viscosity change", "1.0e3",
                             Patterns::Double (0.),
                             "$\\Delta \\eta_{F,max}$: maximum depletion strengthening of viscosity. "
                             "Rationale: melting dehydrates the source rock by removing most of the volatiles,"
                             "and makes it stronger. Hirth and Kohlstedt (1996) report typical values around a "
                             "factor 100 to 1000 viscosity contrast between wet and dry rocks, although some "
                             "experimental studies report a smaller (factor 10) contrast (e.g. Fei et al., 2013).");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    MeltGlobal<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Melt global");
        {
          reference_rho_s                   = prm.get_double ("Reference solid density");
          reference_rho_f                   = prm.get_double ("Reference melt density");
          reference_T                       = prm.get_double ("Reference temperature");
          eta_0                             = prm.get_double ("Reference shear viscosity");
          xi_0                              = prm.get_double ("Reference bulk viscosity");
          eta_f                             = prm.get_double ("Reference melt viscosity");
          reference_permeability            = prm.get_double ("Reference permeability");
          thermal_viscosity_exponent        = prm.get_double ("Thermal viscosity exponent");
          thermal_bulk_viscosity_exponent   = prm.get_double ("Thermal bulk viscosity exponent");
          thermal_conductivity              = prm.get_double ("Thermal conductivity");
          reference_specific_heat           = prm.get_double ("Reference specific heat");
          thermal_expansivity               = prm.get_double ("Thermal expansion coefficient");
          alpha_phi                         = prm.get_double ("Exponential melt weakening factor");
          depletion_density_change          = prm.get_double ("Depletion density change");
          surface_solidus                   = prm.get_double ("Surface solidus");
          depletion_solidus_change          = prm.get_double ("Depletion solidus change");
          pressure_solidus_change           = prm.get_double ("Pressure solidus change");
          compressibility                   = prm.get_double ("Solid compressibility");
          melt_compressibility              = prm.get_double ("Melt compressibility");
          include_melting_and_freezing      = prm.get_bool ("Include melting and freezing");
          melting_time_scale                = prm.get_double ("Melting time scale for operator splitting");
          alpha_depletion                   = prm.get_double ("Exponential depletion strengthening factor");
          delta_eta_depletion_max           = prm.get_double ("Maximum Depletion viscosity change");

          if (thermal_viscosity_exponent!=0.0 && reference_T == 0.0)
            AssertThrow(false, ExcMessage("Error: Material model Melt global with Thermal viscosity exponent can not have reference_T=0."));

          if (this->convert_output_to_years() == true)
            melting_time_scale *= year_in_seconds;

          if (this->get_parameters().use_operator_splitting)
            {
              if (this->get_parameters().reaction_solver_type == Parameters<dim>::ReactionSolverType::fixed_step)
                AssertThrow(melting_time_scale >= this->get_parameters().reaction_time_step,
                            ExcMessage("The reaction time step " + Utilities::to_string(this->get_parameters().reaction_time_step)
                                       + " in the operator splitting scheme is too large to compute melting rates! "
                                       "You have to choose it in such a way that it is smaller than the 'Melting time scale for "
                                       "operator splitting' chosen in the material model, which is currently "
                                       + Utilities::to_string(melting_time_scale) + "."));
              AssertThrow(melting_time_scale > 0,
                          ExcMessage("The Melting time scale for operator splitting must be larger than 0!"));
              AssertThrow(this->introspection().compositional_name_exists("porosity"),
                          ExcMessage("Material model Melt global with melt transport only "
                                     "works if there is a compositional field called porosity."));
            }

          if (this->include_melt_transport())
            {
              AssertThrow(this->introspection().compositional_name_exists("porosity"),
                          ExcMessage("Material model Melt global with melt transport only "
                                     "works if there is a compositional field called porosity."));
              if (include_melting_and_freezing)
                {
                  AssertThrow(this->introspection().compositional_name_exists("peridotite"),
                              ExcMessage("Material model Melt global only works if there is a "
                                         "compositional field called peridotite."));
                }
            }

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    MeltGlobal<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      if (this->get_parameters().use_operator_splitting
          && out.template get_additional_output<ReactionRateOutputs<dim>>() == nullptr)
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
    ASPECT_REGISTER_MATERIAL_MODEL(MeltGlobal,
                                   "melt global",
                                   "A material model that implements a simple formulation of the "
                                   "material parameters required for the modeling of melt transport, "
                                   "including a source term for the porosity according to a simplified "
                                   "linear melting model similar to \\cite{schmeling2006}:\n"
                                   "$\\phi_{\\text{equilibrium}} = \\frac{T-T_{\\text{sol}}}{T_{\\text{liq}}-T_{\\text{sol}}}$\n"
                                   "with "
                                   "$T_{\\text{sol}} = T_{\\text{sol,0}} + \\Delta T_p \\, p + \\Delta T_c \\, C$ \n"
                                   "$T_{\\text{liq}} = T_{\\text{sol}}  + \\Delta T_{\\text{sol-liq}}$.")
  }
}
