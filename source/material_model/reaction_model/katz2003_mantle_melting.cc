/*
  Copyright (C) 2024 by the authors of the ASPECT code.

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


#include <aspect/material_model/reaction_model/katz2003_mantle_melting.h>
#include <aspect/utilities.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <deal.II/base/parameter_handler.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace ReactionModel
    {


      template <int dim>
      Katz2003MantleMelting<dim>::Katz2003MantleMelting()
        = default;



      template <int dim>
      double
      Katz2003MantleMelting<dim>::
      reference_darcy_coefficient () const
      {
        // 0.01 = 1% melt
        return reference_permeability * Utilities::fixed_power<3>(0.01) / viscosity_fluid;
      }


      template <int dim>
      double
      Katz2003MantleMelting<dim>::
      melt_fraction (const double temperature,
                     const double pressure) const
      {
        // anhydrous melting of peridotite after Katz, 2003
        const double T_solidus  = A1 + 273.15
                                  + A2 * pressure
                                  + A3 * pressure * pressure;
        const double T_lherz_liquidus = B1 + 273.15
                                        + B2 * pressure
                                        + B3 * pressure * pressure;
        const double T_liquidus = C1 + 273.15
                                  + C2 * pressure
                                  + C3 * pressure * pressure;

        // melt fraction for peridotite with clinopyroxene
        double peridotite_melt_fraction;
        if (temperature < T_solidus || pressure > 1.3e10)
          peridotite_melt_fraction = 0.0;
        else if (temperature > T_lherz_liquidus)
          peridotite_melt_fraction = 1.0;
        else
          peridotite_melt_fraction = std::pow((temperature - T_solidus) / (T_lherz_liquidus - T_solidus),beta);

        // melt fraction after melting of all clinopyroxene
        const double R_cpx = r1 + r2 * std::max(0.0, pressure);
        const double F_max = M_cpx / R_cpx;

        if (peridotite_melt_fraction > F_max && temperature < T_liquidus)
          {
            const double T_max = std::pow(F_max,1/beta) * (T_lherz_liquidus - T_solidus) + T_solidus;
            peridotite_melt_fraction = F_max + (1 - F_max) * std::pow((temperature - T_max) / (T_liquidus - T_max),beta);
          }
        return peridotite_melt_fraction;
      }



      template <int dim>
      double
      Katz2003MantleMelting<dim>::
      entropy_change (const double temperature,
                      const double pressure,
                      const double maximum_melt_fraction,
                      const NonlinearDependence::Dependence dependence) const
      {
        double entropy_gradient = 0.0;

        // calculate latent heat of melting
        // we need the change of melt fraction in dependence of pressure and temperature

        // for peridotite after Katz, 2003
        const double T_solidus        = A1 + 273.15
                                        + A2 * pressure
                                        + A3 * pressure * pressure;
        const double T_lherz_liquidus = B1 + 273.15
                                        + B2 * pressure
                                        + B3 * pressure * pressure;
        const double T_liquidus       = C1 + 273.15
                                        + C2 * pressure
                                        + C3 * pressure * pressure;

        const double dT_solidus_dp        = A2 + 2 * A3 * pressure;
        const double dT_lherz_liquidus_dp = B2 + 2 * B3 * pressure;
        const double dT_liquidus_dp       = C2 + 2 * C3 * pressure;

        if (temperature > T_solidus && temperature < T_liquidus && pressure < 1.3e10)
          {
            // melt fraction when clinopyroxene is still present
            double melt_fraction_derivative_temperature
              = beta * std::pow((temperature - T_solidus)/(T_lherz_liquidus - T_solidus),beta-1)
                / (T_lherz_liquidus - T_solidus);

            double melt_fraction_derivative_pressure
              = beta * std::pow((temperature - T_solidus)/(T_lherz_liquidus - T_solidus),beta-1)
                * (dT_solidus_dp * (temperature - T_lherz_liquidus)
                   + dT_lherz_liquidus_dp * (T_solidus - temperature))
                / Utilities::fixed_power<2>(T_lherz_liquidus - T_solidus);

            // melt fraction after melting of all clinopyroxene
            const double R_cpx = r1 + r2 * std::max(0.0, pressure);
            const double F_max = M_cpx / R_cpx;

            if (melt_fraction(temperature, pressure) > F_max)
              {
                const double T_max = std::pow(F_max,1.0/beta) * (T_lherz_liquidus - T_solidus) + T_solidus;
                const double dF_max_dp = - M_cpx * Utilities::fixed_power<-2>(r1 + r2 * pressure) * r2;
                const double dT_max_dp = dT_solidus_dp
                                         + 1.0/beta * std::pow(F_max,1.0/beta - 1.0) * dF_max_dp * (T_lherz_liquidus - T_solidus)
                                         + std::pow(F_max,1.0/beta) * (dT_lherz_liquidus_dp - dT_solidus_dp);

                melt_fraction_derivative_temperature
                  = (1.0 - F_max) * beta * std::pow((temperature - T_max)/(T_liquidus - T_max),beta-1)
                    / (T_liquidus - T_max);

                melt_fraction_derivative_pressure
                  = dF_max_dp
                    - dF_max_dp * std::pow((temperature - T_max)/(T_liquidus - T_max),beta)
                    + (1.0 - F_max) * beta * std::pow((temperature - T_max)/(T_liquidus - T_max),beta-1)
                    * (dT_max_dp * (T_max - T_liquidus) - (dT_liquidus_dp - dT_max_dp) * (temperature - T_max)) / Utilities::fixed_power<2>(T_liquidus - T_max);
              }

            double melt_fraction_derivative = 0;
            if (dependence == NonlinearDependence::temperature)
              melt_fraction_derivative = melt_fraction_derivative_temperature;
            else if (dependence == NonlinearDependence::pressure)
              melt_fraction_derivative = melt_fraction_derivative_pressure;
            else
              AssertThrow(false, ExcMessage("not implemented"));

            if (melt_fraction(temperature, pressure) >= maximum_melt_fraction)
              entropy_gradient = melt_fraction_derivative * peridotite_melting_entropy_change;
          }
        return entropy_gradient;
      }



      template <int dim>
      void
      Katz2003MantleMelting<dim>::
      calculate_reaction_rate_outputs(const typename Interface<dim>::MaterialModelInputs &in,
                                      typename Interface<dim>::MaterialModelOutputs &out) const
      {
        ReactionRateOutputs<dim> *reaction_rate_out = out.template get_additional_output<ReactionRateOutputs<dim>>();

        for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
          {

            if (this->include_melt_transport())
              {
                const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
                const unsigned int peridotite_idx = this->introspection().compositional_index_for_name("peridotite");
                const double old_porosity = in.composition[i][porosity_idx];
                const double maximum_melt_fraction = in.composition[i][peridotite_idx];

                // calculate the melting rate as difference between the equilibrium melt fraction
                // and the solution of the previous time step
                double porosity_change = 0.0;
                if (fractional_melting)
                  {
                    // solidus is lowered by previous melting events (fractional melting)
                    const double solidus_change = (maximum_melt_fraction - old_porosity) * depletion_solidus_change;
                    const double eq_melt_fraction = melt_fraction(in.temperature[i] - solidus_change, this->get_adiabatic_conditions().pressure(in.position[i]));
                    porosity_change = eq_melt_fraction - old_porosity;
                  }
                else
                  {
                    // batch melting
                    porosity_change = melt_fraction(in.temperature[i], this->get_adiabatic_conditions().pressure(in.position[i]))
                                      - std::max(maximum_melt_fraction, 0.0);
                    porosity_change = std::max(porosity_change, 0.0);

                    // freezing of melt below the solidus

                    // If the porosity is larger than the equilibrium melt fraction, melt should freeze again.
                    // Because we do not track the melt composition, we have to use a workaround here for freezing of melt:
                    // We reduce the porosity until either it reaches the equilibrium melt fraction, or the depletion
                    // (peridotite field), which decreases as melt freezes, reaches the same value as the equilibrium
                    // melt fraction, whatever happens earlier. An exception is when the melt fraction is zero; in this case
                    // all melt should freeze.
                    const double eq_melt_fraction = melt_fraction(in.temperature[i], this->get_adiabatic_conditions().pressure(in.position[i]));

                    // If the porosity change is not negative, there is no freezing, and the change in porosity
                    // is covered by the melting relation above.

                    // porosity reaches the equilibrium melt fraction:
                    const double porosity_change_wrt_melt_fraction = std::min(eq_melt_fraction - old_porosity - porosity_change,0.0);

                    // depletion reaches the equilibrium melt fraction:
                    const double porosity_change_wrt_depletion = std::min((eq_melt_fraction - std::max(maximum_melt_fraction, 0.0))
                                                                          * (1.0 - old_porosity) / (1.0 - maximum_melt_fraction),0.0);
                    double freezing_amount = std::max(porosity_change_wrt_melt_fraction, porosity_change_wrt_depletion);

                    if (eq_melt_fraction == 0.0)
                      freezing_amount = - old_porosity;

                    porosity_change += freezing_amount;

                    // Adapt time scale of freezing with respect to melting.
                    // We have to multiply with the melting time scale here to obtain the porosity change
                    // that happens in the time defined by the melting time scale (as opposed to a rate).
                    // This is important because we want to perform some checks on this quantity (for example,
                    // we want to make sure that this change does not lead to a negative porosity, see below).
                    // Later on, the overall porosity change is then divided again by the melting time scale
                    // to obtain the rate of melting or freezing, which is used in the operator splitting scheme.
                    if (porosity_change < 0 )
                      porosity_change *= freezing_rate * melting_time_scale;

                  }

                // remove melt that gets close to the surface
                if (this->get_geometry_model().depth(in.position[i]) < extraction_depth)
                  porosity_change = -old_porosity * (in.position[i](1) - (this->get_geometry_model().maximal_depth() - extraction_depth))/extraction_depth;

                // do not allow negative porosity
                porosity_change = std::max(porosity_change, -old_porosity);

                // because depletion is a volume-based, and not a mass-based property that is advected,
                // additional scaling factors on the right hand side apply
                for (unsigned int c=0; c<in.composition[i].size(); ++c)
                  {
                    // fill reaction rate outputs
                    if (reaction_rate_out != nullptr && in.requests_property(MaterialProperties::reaction_rates))
                      {
                        if (c == peridotite_idx && this->get_timestep_number() > 0)
                          reaction_rate_out->reaction_rates[i][c] = porosity_change / melting_time_scale * (1 - maximum_melt_fraction) / (1 - old_porosity);
                        else if (c == porosity_idx && this->get_timestep_number() > 0)
                          reaction_rate_out->reaction_rates[i][c] = porosity_change / melting_time_scale;
                        else
                          reaction_rate_out->reaction_rates[i][c] = 0.0;
                      }
                    out.reaction_terms[i][c] = 0.0;
                  }

                out.entropy_derivative_pressure[i]    = entropy_change (in.temperature[i], this->get_adiabatic_conditions().pressure(in.position[i]), maximum_melt_fraction, NonlinearDependence::pressure);
                out.entropy_derivative_temperature[i] = entropy_change (in.temperature[i], this->get_adiabatic_conditions().pressure(in.position[i]), maximum_melt_fraction, NonlinearDependence::temperature);
              }
            else
              {
                out.entropy_derivative_pressure[i]    = entropy_change (in.temperature[i], this->get_adiabatic_conditions().pressure(in.position[i]), 0, NonlinearDependence::pressure);
                out.entropy_derivative_temperature[i] = entropy_change (in.temperature[i], this->get_adiabatic_conditions().pressure(in.position[i]), 0, NonlinearDependence::temperature);

                // no melting/freezing is used in the model --> set all reactions to zero
                for (unsigned int c=0; c<in.composition[i].size(); ++c)
                  {
                    out.reaction_terms[i][c] = 0.0;

                    if (reaction_rate_out != nullptr)
                      reaction_rate_out->reaction_rates[i][c] = 0.0;
                  }
              }
          }

      }


      template <int dim>
      void
      Katz2003MantleMelting<dim>::
      calculate_fluid_outputs(const typename Interface<dim>::MaterialModelInputs &in,
                              typename Interface<dim>::MaterialModelOutputs &out,
                              const double reference_T) const
      {
        MeltOutputs<dim> *melt_out = out.template get_additional_output<MeltOutputs<dim>>();
        const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");

        if (melt_out != nullptr)
          {

            for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
              {
                double porosity = std::max(in.composition[i][porosity_idx],0.0);

                melt_out->fluid_viscosities[i] = viscosity_fluid;
                melt_out->permeabilities[i] = reference_permeability * Utilities::fixed_power<3>(porosity) * Utilities::fixed_power<2>(1.0-porosity);

                // first, calculate temperature dependence of density
                double temperature_dependence = 1.0;
                if (this->include_adiabatic_heating ())
                  {
                    // temperature dependence is 1 - alpha * (T - T(adiabatic))
                    temperature_dependence -= (in.temperature[i] - this->get_adiabatic_conditions().temperature(in.position[i]))
                                              * out.thermal_expansion_coefficients[i];
                  }
                else
                  temperature_dependence -= (in.temperature[i] - reference_T) * out.thermal_expansion_coefficients[i];

                // the fluid compressibility includes two parts, a constant compressibility, and a pressure-dependent one
                // this is a simplified formulation, experimental data are often fit to the Birch-Murnaghan equation of state
                const double fluid_compressibility = melt_compressibility / (1.0 + in.pressure[i] * melt_bulk_modulus_derivative * melt_compressibility);

                melt_out->fluid_densities[i] = reference_rho_fluid * std::exp(fluid_compressibility * (in.pressure[i] - this->get_surface_pressure()))
                                               * temperature_dependence;

                melt_out->fluid_density_gradients[i] = melt_out->fluid_densities[i] * melt_out->fluid_densities[i]
                                                       * fluid_compressibility
                                                       * this->get_gravity_model().gravity_vector(in.position[i]);

                const double phi_0 = 0.05;
                porosity = std::max(std::min(porosity,0.995),1e-4);
                melt_out->compaction_viscosities[i] = xi_0 * phi_0 / porosity;

                double visc_temperature_dependence = 1.0;
                if (this->include_adiabatic_heating ())
                  {
                    const double delta_temp = in.temperature[i]-this->get_adiabatic_conditions().temperature(in.position[i]);
                    visc_temperature_dependence = std::max(std::min(std::exp(-thermal_bulk_viscosity_exponent*delta_temp/this->get_adiabatic_conditions().temperature(in.position[i])),1e4),1e-4);
                  }
                else
                  {
                    const double delta_temp = in.temperature[i]-reference_T;
                    const double T_dependence = (thermal_bulk_viscosity_exponent == 0.0
                                                 ?
                                                 0.0
                                                 :
                                                 thermal_bulk_viscosity_exponent*delta_temp/reference_T);
                    visc_temperature_dependence = std::max(std::min(std::exp(-T_dependence),1e4),1e-4);
                  }
                melt_out->compaction_viscosities[i] *= visc_temperature_dependence;


              }
          }

        if (this->include_melt_transport() && in.requests_property(MaterialProperties::viscosity))
          {
            for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
              {
                const double porosity = std::min(1.0, std::max(in.composition[i][porosity_idx],0.0));
                out.viscosities[i] *= std::exp(- alpha_phi * porosity);
              }
          }
      }



      template <int dim>
      void
      Katz2003MantleMelting<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("A1", "1085.7",
                           Patterns::Double (),
                           "Constant parameter in the quadratic "
                           "function that approximates the solidus "
                           "of peridotite. "
                           "Units: \\si{\\degreeCelsius}.");
        prm.declare_entry ("A2", "1.329e-7",
                           Patterns::Double (),
                           "Prefactor of the linear pressure term "
                           "in the quadratic function that approximates "
                           "the solidus of peridotite. "
                           "Units: \\si{\\degreeCelsius\\per\\pascal}.");
        prm.declare_entry ("A3", "-5.1e-18",
                           Patterns::Double (),
                           "Prefactor of the quadratic pressure term "
                           "in the quadratic function that approximates "
                           "the solidus of peridotite. "
                           "Units: \\si{\\degreeCelsius\\per\\pascal\\squared}.");
        prm.declare_entry ("B1", "1475.0",
                           Patterns::Double (),
                           "Constant parameter in the quadratic "
                           "function that approximates the lherzolite "
                           "liquidus used for calculating the fraction "
                           "of peridotite-derived melt. "
                           "Units: \\si{\\degreeCelsius}.");
        prm.declare_entry ("B2", "8.0e-8",
                           Patterns::Double (),
                           "Prefactor of the linear pressure term "
                           "in the quadratic function that approximates "
                           "the  lherzolite liquidus used for "
                           "calculating the fraction of peridotite-"
                           "derived melt. "
                           "Units: \\si{\\degreeCelsius\\per\\pascal}.");
        prm.declare_entry ("B3", "-3.2e-18",
                           Patterns::Double (),
                           "Prefactor of the quadratic pressure term "
                           "in the quadratic function that approximates "
                           "the  lherzolite liquidus used for "
                           "calculating the fraction of peridotite-"
                           "derived melt. "
                           "Units: \\si{\\degreeCelsius\\per\\pascal\\squared}.");
        prm.declare_entry ("C1", "1780.0",
                           Patterns::Double (),
                           "Constant parameter in the quadratic "
                           "function that approximates the liquidus "
                           "of peridotite. "
                           "Units: \\si{\\degreeCelsius}.");
        prm.declare_entry ("C2", "4.50e-8",
                           Patterns::Double (),
                           "Prefactor of the linear pressure term "
                           "in the quadratic function that approximates "
                           "the liquidus of peridotite. "
                           "Units: \\si{\\degreeCelsius\\per\\pascal}.");
        prm.declare_entry ("C3", "-2.0e-18",
                           Patterns::Double (),
                           "Prefactor of the quadratic pressure term "
                           "in the quadratic function that approximates "
                           "the liquidus of peridotite. "
                           "Units: \\si{\\degreeCelsius\\per\\pascal\\squared}.");
        prm.declare_entry ("r1", "0.5",
                           Patterns::Double (),
                           "Constant in the linear function that "
                           "approximates the clinopyroxene reaction "
                           "coefficient. "
                           "Units: non-dimensional.");
        prm.declare_entry ("r2", "8e-11",
                           Patterns::Double (),
                           "Prefactor of the linear pressure term "
                           "in the linear function that approximates "
                           "the clinopyroxene reaction coefficient. "
                           "Units: \\si{\\per\\pascal}.");
        prm.declare_entry ("beta", "1.5",
                           Patterns::Double (),
                           "Exponent of the melting temperature in "
                           "the melt fraction calculation. "
                           "Units: non-dimensional.");
        prm.declare_entry ("Mass fraction cpx", "0.15",
                           Patterns::Double (),
                           "Mass fraction of clinopyroxene in the "
                           "peridotite to be molten. "
                           "Units: non-dimensional.");
        prm.declare_entry ("Peridotite melting entropy change", "-300.",
                           Patterns::Double (),
                           "The entropy change for the phase transition "
                           "from solid to melt of peridotite. "
                           "Units: \\si{\\joule\\per\\kelvin\\per\\kilogram}.");
        prm.declare_entry ("Reference melt density", "2500.",
                           Patterns::Double (0.),
                           "Reference density of the melt/fluid$\\rho_{f,0}$. "
                           "Units: \\si{\\kilogram\\per\\meter\\cubed}.");
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
        prm.declare_entry ("Thermal bulk viscosity exponent", "0.0",
                           Patterns::Double (0.),
                           "The temperature dependence of the bulk viscosity. Dimensionless exponent. "
                           "See the general documentation "
                           "of this model for a formula that states the dependence of the "
                           "viscosity on this factor, which is called $\\beta$ there.");
        prm.declare_entry ("Melt extraction depth", "1000.0",
                           Patterns::Double (0.),
                           "Depth above that melt will be extracted from the model, "
                           "which is done by a negative reaction term proportional to the "
                           "porosity field. "
                           "Units: \\si{\\meter}.");
        prm.declare_entry ("Melt compressibility", "0.0",
                           Patterns::Double (0.),
                           "The value of the compressibility of the melt. "
                           "Units: \\si{\\per\\pascal}.");
        prm.declare_entry ("Melt bulk modulus derivative", "0.0",
                           Patterns::Double (0.),
                           "The value of the pressure derivative of the melt bulk "
                           "modulus. "
                           "Units: None.");
        prm.declare_entry ("Use fractional melting", "false",
                           Patterns::Bool (),
                           "If fractional melting should be used (if true), including a solidus "
                           "change based on depletion (in this case, the amount of melt that has "
                           "migrated away from its origin), and freezing of melt when it has moved "
                           "to a region with temperatures lower than the solidus; or if batch "
                           "melting should be used (if false), assuming that the melt fraction only "
                           "depends on temperature and pressure, and how much melt has already been "
                           "generated at a given point, but not considering movement of melt in "
                           "the melting parameterization."
                           "\n\n"
                           "Note that melt does not freeze unless the 'Freezing rate' parameter is set "
                           "to a value larger than 0.");
        prm.declare_entry ("Freezing rate", "0.0",
                           Patterns::Double (0.),
                           "Freezing rate of melt when in subsolidus regions. "
                           "If this parameter is set to a number larger than 0.0, it specifies the "
                           "fraction of melt that will freeze per year (or per second, depending on the "
                           "``Use years in output instead of seconds'' parameter), as soon as the porosity "
                           "exceeds the equilibrium melt fraction, and the equilibrium melt fraction "
                           "falls below the depletion. In this case, melt will freeze according to the "
                           "given rate until one of those conditions is not fulfilled anymore. The "
                           "reasoning behind this is that there should not be more melt present than "
                           "the equilibrium melt fraction, as melt production decreases with increasing "
                           "depletion, but the freezing process of melt also reduces the depletion by "
                           "the same amount, and as soon as the depletion falls below the equilibrium "
                           "melt fraction, we expect that material should melt again (no matter how "
                           "much melt is present). This is quite a simplification and not a realistic "
                           "freezing parameterization, but without tracking the melt composition, there "
                           "is no way to compute freezing rates accurately. "
                           "If this parameter is set to zero, no freezing will occur. "
                           "Note that freezing can never be faster than determined by the "
                           "``Melting time scale for operator splitting''. The product of the "
                           "``Freezing rate'' and the ``Melting time scale for operator splitting'' "
                           "defines how fast freezing occurs with respect to melting (if the "
                           "product is 0.5, melting will occur twice as fast as freezing). "
                           "Units: 1/yr or 1/s, depending on the ``Use years "
                           "in output instead of seconds'' parameter.");
        prm.declare_entry ("Melting time scale for operator splitting", "1e3",
                           Patterns::Double (0.),
                           "Because the operator splitting scheme is used, the porosity field can not "
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
                           "computed. "
                           "Units: yr or s, depending on the ``Use years in output instead of seconds'' parameter.");
        prm.declare_entry ("Depletion solidus change", "200.0",
                           Patterns::Double (0.),
                           "The solidus temperature change for a depletion of 100\\%. For positive "
                           "values, the solidus gets increased for a positive peridotite field "
                           "(depletion) and lowered for a negative peridotite field (enrichment). "
                           "Scaling with depletion is linear. Only active when fractional melting "
                           "is used. "
                           "Units: \\si{\\kelvin}.");
        prm.declare_entry ("Reference permeability", "1e-8",
                           Patterns::Double(),
                           "Reference permeability of the solid host rock."
                           "Units: \\si{\\meter\\squared}.");


      }


      template <int dim>
      void
      Katz2003MantleMelting<dim>::parse_parameters (ParameterHandler &prm)
      {
        A1    = prm.get_double ("A1");
        A2    = prm.get_double ("A2");
        A3    = prm.get_double ("A3");
        B1    = prm.get_double ("B1");
        B2    = prm.get_double ("B2");
        B3    = prm.get_double ("B3");
        C1    = prm.get_double ("C1");
        C2    = prm.get_double ("C2");
        C3    = prm.get_double ("C3");
        r1    = prm.get_double ("r1");
        r2    = prm.get_double ("r2");
        beta  = prm.get_double ("beta");
        M_cpx = prm.get_double ("Mass fraction cpx");
        peridotite_melting_entropy_change
          = prm.get_double ("Peridotite melting entropy change");

        reference_rho_fluid        = prm.get_double ("Reference melt density");
        xi_0                       = prm.get_double ("Reference bulk viscosity");
        viscosity_fluid            = prm.get_double ("Reference melt viscosity");
        thermal_bulk_viscosity_exponent = prm.get_double ("Thermal bulk viscosity exponent");
        alpha_phi                  = prm.get_double ("Exponential melt weakening factor");
        extraction_depth           = prm.get_double ("Melt extraction depth");
        melt_compressibility       = prm.get_double ("Melt compressibility");
        fractional_melting         = prm.get_bool ("Use fractional melting");
        freezing_rate              = prm.get_double ("Freezing rate");
        melting_time_scale         = prm.get_double ("Melting time scale for operator splitting");
        melt_bulk_modulus_derivative = prm.get_double ("Melt bulk modulus derivative");
        depletion_solidus_change   = prm.get_double ("Depletion solidus change");
        reference_permeability     = prm.get_double ("Reference permeability");


        if (this->convert_output_to_years() == true)
          {
            melting_time_scale *= year_in_seconds;
            freezing_rate /= year_in_seconds;
          }

        AssertThrow(melting_time_scale > 0,
                    ExcMessage("The Melting time scale for operator splitting must be larger than 0!"));

        if (this->get_parameters().reaction_solver_type == Parameters<dim>::ReactionSolverType::fixed_step)
          {
            AssertThrow(melting_time_scale >= this->get_parameters().reaction_time_step,
                        ExcMessage("The reaction time step " + Utilities::to_string(this->get_parameters().reaction_time_step)
                                   + " in the operator splitting scheme is too large to compute melting rates! "
                                   "You have to choose it in such a way that it is smaller than the 'Melting time scale for "
                                   "operator splitting' chosen in the material model, which is currently "
                                   + Utilities::to_string(melting_time_scale) + "."));

            AssertThrow(freezing_rate * this->get_parameters().reaction_time_step <= 1.0,
                        ExcMessage("The reaction time step " + Utilities::to_string(this->get_parameters().reaction_time_step)
                                   + " in the operator splitting scheme is too large to compute freezing rates! "
                                   "You have to choose it in such a way that it is smaller than the inverse of the "
                                   "'Freezing rate' chosen in the material model, which is currently "
                                   + Utilities::to_string(1.0/freezing_rate) + "."));
          }
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
#define INSTANTIATE(dim) \
  namespace ReactionModel \
  { \
    template class Katz2003MantleMelting<dim>; \
  }

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
