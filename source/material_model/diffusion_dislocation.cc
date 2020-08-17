/*
  Copyright (C) 2011 - 2020 by the authors of the ASPECT code.

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

#include <aspect/material_model/diffusion_dislocation.h>
#include <aspect/utilities.h>
#include <aspect/adiabatic_conditions/interface.h>

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    std::vector<double>
    DiffusionDislocation<dim>::
    calculate_isostrain_viscosities ( const std::vector<double> &volume_fractions,
                                      const double &pressure,
                                      const double &temperature,
                                      const SymmetricTensor<2,dim> &strain_rate) const
    {
      // This function calculates viscosities assuming that all the compositional fields
      // experience the same strain rate (isostrain).

      // If strain rate is zero (like during the first time step) set it to some very small number
      // to prevent a division-by-zero, and a floating point exception.
      // Otherwise, calculate the square-root of the norm of the second invariant of the deviatoric-
      // strain rate (often simplified as epsilondot_ii)
      const double edot_ii = std::max(std::sqrt(std::fabs(second_invariant(deviator(strain_rate)))),
                                      min_strain_rate);


      // Find effective viscosities for each of the individual phases
      // Viscosities should have same number of entries as compositional fields
      std::vector<double> composition_viscosities(volume_fractions.size());
      for (unsigned int j=0; j < volume_fractions.size(); ++j)
        {
          // Power law creep equation
          // edot_ii_i = A_i * stress_ii_i^{n_i} * d^{-m} \exp\left(-\frac{E_i^\ast + PV_i^\ast}{n_iRT}\right)
          // where ii indicates the square root of the second invariant and
          // i corresponds to diffusion or dislocation creep

          // For diffusion creep, viscosity is grain size dependent
          const Rheology::DiffusionCreepParameters diffusion_creep_parameters = diffusion_creep.compute_creep_parameters(j);

          // For dislocation creep, viscosity is grain size independent (m=0)
          const Rheology::DislocationCreepParameters dislocation_creep_parameters = dislocation_creep.compute_creep_parameters(j);

          // For diffusion creep, viscosity is grain size dependent
          const double prefactor_stress_diffusion = diffusion_creep_parameters.prefactor *
                                                    std::pow(grain_size, -diffusion_creep_parameters.grain_size_exponent) *
                                                    std::exp(-(std::max(diffusion_creep_parameters.activation_energy + pressure*diffusion_creep_parameters.activation_volume,0.0))/
                                                             (constants::gas_constant*temperature));

          // Because the ratios of the diffusion and dislocation strain rates are not known, stress is also unknown
          // We use Newton's method to find the second invariant of the stress tensor.
          // Start with the assumption that all strain is accommodated by diffusion creep:
          // If the diffusion creep prefactor is very small, that means that the diffusion viscosity is very large.
          // In this case, use the maximum viscosity instead to compute the starting guess.
          double stress_ii = (prefactor_stress_diffusion > (0.5 / max_visc)
                              ?
                              edot_ii/prefactor_stress_diffusion
                              :
                              0.5 / max_visc);
          double strain_rate_residual = 2*strain_rate_residual_threshold;
          double strain_rate_deriv = 0;
          unsigned int stress_iteration = 0;
          while (std::abs(strain_rate_residual) > strain_rate_residual_threshold
                 && stress_iteration < stress_max_iteration_number)
            {

              const std::pair<double, double> diff_edot_and_deriv = diffusion_creep.compute_strain_rate_and_derivative(stress_ii, pressure, temperature, diffusion_creep_parameters);
              const std::pair<double, double> disl_edot_and_deriv = dislocation_creep.compute_strain_rate_and_derivative(stress_ii, pressure, temperature, dislocation_creep_parameters);

              strain_rate_residual = diff_edot_and_deriv.first + disl_edot_and_deriv.first - edot_ii;
              strain_rate_deriv = diff_edot_and_deriv.second + disl_edot_and_deriv.second ;

              // If the strain rate derivative is zero, we catch it below.
              if (strain_rate_deriv>std::numeric_limits<double>::min())
                stress_ii -= strain_rate_residual/strain_rate_deriv;
              stress_iteration += 1;

              // In case the Newton iteration does not succeed, we do a fixpoint iteration.
              // This allows us to bound both the diffusion and dislocation viscosity
              // between a minimum and maximum value, so that we can compute the correct
              // viscosity values even if the parameters lead to one or both of the
              // viscosities being essentially zero or infinity.
              // If anything that would be used in the next iteration is not finite, the
              // Newton iteration would trigger an exception and we want to do the fixpoint
              // iteration instead.
              const bool abort_newton_iteration = !numbers::is_finite(stress_ii)
                                                  || !numbers::is_finite(strain_rate_residual)
                                                  || !numbers::is_finite(strain_rate_deriv)
                                                  || strain_rate_deriv < std::numeric_limits<double>::min()
                                                  || !numbers::is_finite(std::pow(stress_ii, diffusion_creep_parameters.stress_exponent-1))
                                                  || !numbers::is_finite(std::pow(stress_ii, dislocation_creep_parameters.stress_exponent-1))
                                                  || stress_iteration == stress_max_iteration_number;
              if (abort_newton_iteration)
                {
                  double diffusion_strain_rate = edot_ii;
                  double dislocation_strain_rate = min_strain_rate;
                  stress_iteration = 0;

                  do
                    {
                      const double old_diffusion_strain_rate = diffusion_strain_rate;

                      const double diffusion_prefactor = 0.5 * std::pow(diffusion_creep_parameters.prefactor,-1.0/diffusion_creep_parameters.stress_exponent);
                      const double diffusion_grain_size_dependence = std::pow(grain_size, diffusion_creep_parameters.grain_size_exponent/diffusion_creep_parameters.stress_exponent);
                      const double diffusion_strain_rate_dependence = std::pow(diffusion_strain_rate, (1.-diffusion_creep_parameters.stress_exponent)/diffusion_creep_parameters.stress_exponent);
                      const double diffusion_T_and_P_dependence = std::exp(std::max(diffusion_creep_parameters.activation_energy + pressure*diffusion_creep_parameters.activation_volume,0.0)/
                                                                           (constants::gas_constant*temperature));

                      const double diffusion_viscosity = std::min(std::max(diffusion_prefactor * diffusion_grain_size_dependence
                                                                           * diffusion_strain_rate_dependence * diffusion_T_and_P_dependence,
                                                                           min_visc), max_visc);

                      const double dislocation_prefactor = 0.5 * std::pow(dislocation_creep_parameters.prefactor,-1.0/dislocation_creep_parameters.stress_exponent);
                      const double dislocation_strain_rate_dependence = std::pow(dislocation_strain_rate, (1.-dislocation_creep_parameters.stress_exponent)/dislocation_creep_parameters.stress_exponent);
                      const double dislocation_T_and_P_dependence = std::exp(std::max(dislocation_creep_parameters.activation_energy + pressure*dislocation_creep_parameters.activation_volume,0.0)/
                                                                             (dislocation_creep_parameters.stress_exponent*constants::gas_constant*temperature));

                      const double dislocation_viscosity = std::min(std::max(dislocation_prefactor * dislocation_strain_rate_dependence
                                                                             * dislocation_T_and_P_dependence,
                                                                             min_visc), max_visc);

                      diffusion_strain_rate = dislocation_viscosity / (diffusion_viscosity + dislocation_viscosity) * edot_ii;
                      dislocation_strain_rate = diffusion_viscosity / (diffusion_viscosity + dislocation_viscosity) * edot_ii;

                      stress_iteration++;
                      AssertThrow(stress_iteration < stress_max_iteration_number,
                                  ExcMessage("No convergence has been reached in the loop that determines "
                                             "the ratio of diffusion/dislocation viscosity. Aborting! "
                                             "Residual is " + Utilities::to_string(strain_rate_residual) +
                                             " after " + Utilities::to_string(stress_iteration) + " iterations. "
                                             "You can increase the number of iterations by adapting the "
                                             "parameter 'Maximum strain rate ratio iterations'."));

                      strain_rate_residual = std::abs((diffusion_strain_rate-old_diffusion_strain_rate) / diffusion_strain_rate);
                      stress_ii = 2.0 * edot_ii * 1./(1./diffusion_viscosity + 1./dislocation_viscosity);
                    }
                  while (strain_rate_residual > strain_rate_residual_threshold);

                  break;
                }
            }

          // The effective viscosity, with minimum and maximum bounds
          composition_viscosities[j] = std::min(std::max(stress_ii/edot_ii/2, min_visc), max_visc);
        }
      return composition_viscosities;
    }

    template <int dim>
    void
    DiffusionDislocation<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
        {
          // const Point<dim> position = in.position[i];
          const double temperature = in.temperature[i];
          const double pressure= in.pressure[i];
          const std::vector<double> composition = in.composition[i];
          const std::vector<double> volume_fractions = MaterialUtilities::compute_volume_fractions(composition);

          // Averaging composition-field dependent properties

          // densities
          double density = 0.0;
          for (unsigned int j=0; j < volume_fractions.size(); ++j)
            {
              // not strictly correct if thermal expansivities are different, since we are interpreting
              // these compositions as volume fractions, but the error introduced should not be too bad.
              const double temperature_factor= (1.0 - thermal_expansivities[j] * (temperature - reference_T));
              density += volume_fractions[j] * densities[j] * temperature_factor;
            }

          // thermal expansivities
          double thermal_expansivity = 0.0;
          for (unsigned int j=0; j < volume_fractions.size(); ++j)
            thermal_expansivity += volume_fractions[j] * thermal_expansivities[j];

          // calculate effective viscosity
          if (in.requests_property(MaterialProperties::viscosity))
            {
              // Currently, the viscosities for each of the compositional fields are calculated assuming
              // isostrain amongst all compositions, allowing calculation of the viscosity ratio.
              // TODO: This is only consistent with viscosity averaging if the arithmetic averaging
              // scheme is chosen. It would be useful to have a function to calculate isostress viscosities.
              const std::vector<double> composition_viscosities =
                calculate_isostrain_viscosities(volume_fractions, pressure, temperature, in.strain_rate[i]);

              // The isostrain condition implies that the viscosity averaging should be arithmetic (see above).
              // We have given the user freedom to apply alternative bounds, because in diffusion-dominated
              // creep (where n_diff=1) viscosities are stress and strain-rate independent, so the calculation
              // of compositional field viscosities is consistent with any averaging scheme.
              out.viscosities[i] = MaterialUtilities::average_value(volume_fractions, composition_viscosities, viscosity_averaging);
            }

          out.densities[i] = density;
          out.thermal_expansion_coefficients[i] = thermal_expansivity;
          // Specific heat at the given positions.
          out.specific_heat[i] = heat_capacity;
          // Thermal conductivity at the given positions. If the temperature equation uses
          // the reference density profile formulation, use the reference density to
          // calculate thermal conductivity. Otherwise, use the real density. If the adiabatic
          // conditions are not yet initialized, the real density will still be used.
          if (this->get_parameters().formulation_temperature_equation ==
              Parameters<dim>::Formulation::TemperatureEquation::reference_density_profile &&
              this->get_adiabatic_conditions().is_initialized())
            out.thermal_conductivities[i] = thermal_diffusivity * heat_capacity *
                                            this->get_adiabatic_conditions().density(in.position[i]);
          else
            out.thermal_conductivities[i] = thermal_diffusivity * heat_capacity * density;
          // Compressibility at the given positions.
          // The compressibility is given as
          // $\frac 1\rho \frac{\partial\rho}{\partial p}$.
          out.compressibilities[i] = 0.0;
          // Pressure derivative of entropy at the given positions.
          out.entropy_derivative_pressure[i] = 0.0;
          // Temperature derivative of entropy at the given positions.
          out.entropy_derivative_temperature[i] = 0.0;
          // Change in composition due to chemical reactions at the
          // given positions. The term reaction_terms[i][c] is the
          // change in compositional field c at point i.
          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            out.reaction_terms[i][c] = 0.0;
        }
    }

    template <int dim>
    double
    DiffusionDislocation<dim>::
    reference_viscosity () const
    {
      return ref_visc;
    }

    template <int dim>
    bool
    DiffusionDislocation<dim>::
    is_compressible () const
    {
      return false;
    }

    template <int dim>
    void
    DiffusionDislocation<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection ("Diffusion dislocation");
        {
          // Reference and minimum/maximum values
          prm.declare_entry ("Reference temperature", "293.", Patterns::Double(0.),
                             "For calculating density by thermal expansivity. Units: \\si{\\kelvin}.");
          prm.declare_entry ("Minimum strain rate", "1.4e-20", Patterns::Double(0.),
                             "Stabilizes strain dependent viscosity. Units: \\si{\\per\\second}.");
          prm.declare_entry ("Minimum viscosity", "1e17", Patterns::Double(0.),
                             "Lower cutoff for effective viscosity. Units: \\si{\\pascal\\second}.");
          prm.declare_entry ("Maximum viscosity", "1e28", Patterns::Double(0.),
                             "Upper cutoff for effective viscosity. Units: \\si{\\pascal\\second}.");
          prm.declare_entry ("Effective viscosity coefficient", "1.0", Patterns::Double(0.),
                             "Scaling coefficient for effective viscosity.");
          prm.declare_entry ("Reference viscosity", "1e22", Patterns::Double(0.),
                             "The reference viscosity that is used for pressure scaling. "
                             "To understand how pressure scaling works, take a look at "
                             "\\cite{KHB12}. In particular, the value of this parameter "
                             "would not affect the solution computed by \\aspect{} if "
                             "we could do arithmetic exactly; however, computers do "
                             "arithmetic in finite precision, and consequently we need to "
                             "scale quantities in ways so that their magnitudes are "
                             "roughly the same. As explained in \\cite{KHB12}, we scale "
                             "the pressure during some computations (never visible by "
                             "users) by a factor that involves a reference viscosity. This "
                             "parameter describes this reference viscosity."
                             "\n\n"
                             "For problems with a constant viscosity, you will generally want "
                             "to choose the reference viscosity equal to the actual viscosity. "
                             "For problems with a variable viscosity, the reference viscosity "
                             "should be a value that adequately represents the order of "
                             "magnitude of the viscosities that appear, such as an average "
                             "value or the value one would use to compute a Rayleigh number."
                             "\n\n"
                             "Units: \\si{\\pascal\\second}.");

          // Viscosity iteration parameters
          prm.declare_entry ("Strain rate residual tolerance", "1e-22", Patterns::Double(0.),
                             "Tolerance for correct diffusion/dislocation strain rate ratio.");
          prm.declare_entry ("Maximum strain rate ratio iterations", "40", Patterns::Integer(0),
                             "Maximum number of iterations to find the correct "
                             "diffusion/dislocation strain rate ratio.");

          // Equation of state parameters
          prm.declare_entry ("Thermal diffusivity", "0.8e-6", Patterns::Double(0.),
                             "Units: \\si{\\meter\\squared\\per\\second}.");
          prm.declare_entry ("Heat capacity", "1.25e3",
                             Patterns::Double(0.),
                             "The value of the specific heat $C_p$. "
                             "Units: \\si{\\joule\\per\\kelvin\\per\\kilogram}.");
          prm.declare_entry ("Densities", "3300.",
                             Patterns::List(Patterns::Double(0.)),
                             "List of densities, $\\rho$, for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value. "
                             "Units: \\si{\\kilogram\\per\\meter\\cubed}.");
          prm.declare_entry ("Thermal expansivities", "3.5e-5",
                             Patterns::List(Patterns::Double(0.)),
                             "List of thermal expansivities for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value.  Units: \\si{\\per\\kelvin}.");

          // Rheological parameters
          prm.declare_entry ("Viscosity averaging scheme", "harmonic",
                             Patterns::Selection("arithmetic|harmonic|geometric|maximum composition"),
                             "When more than one compositional field is present at a point "
                             "with different viscosities, we need to come up with an average "
                             "viscosity at that point.  Select a weighted harmonic, arithmetic, "
                             "geometric, or maximum composition.");
          // Diffusion creep parameters
          prm.declare_entry ("Prefactors for diffusion creep", "1.5e-15",
                             Patterns::List(Patterns::Double(0.)),
                             "List of viscosity prefactors, $A$, for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value. "
                             "Units: \\si{\\per\\pascal} \\si{\\meter}$^{m_{\\text{diffusion}}}$ \\si{\\per\\second}.");
          prm.declare_entry ("Stress exponents for diffusion creep", "1.",
                             Patterns::List(Patterns::Double(0.)),
                             "List of stress exponents, $n_{\\text{diffusion}}$, for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value.  Units: None.");
          prm.declare_entry ("Grain size exponents for diffusion creep", "3.",
                             Patterns::List(Patterns::Double(0.)),
                             "List of grain size exponents, $m_{\\text{diffusion}}$, for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value.  Units: None.");
          prm.declare_entry ("Activation energies for diffusion creep", "375e3",
                             Patterns::List(Patterns::Double(0.)),
                             "List of activation energies, $E_a$, for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value. "
                             "Units: \\si{\\joule\\per\\mole}.");
          prm.declare_entry ("Activation volumes for diffusion creep", "6e-6",
                             Patterns::List(Patterns::Double(0.)),
                             "List of activation volumes, $V_a$, for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value. "
                             "Units: \\si{\\meter\\cubed\\per\\mole}.");

          // Dislocation creep parameters
          prm.declare_entry ("Prefactors for dislocation creep", "1.1e-16",
                             Patterns::List(Patterns::Double(0.)),
                             "List of viscosity prefactors, $A$, for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value. "
                             "Units: \\si{\\pascal}$^{-n_{\\text{dislocation}}}$\\si{\\per\\second}.");
          prm.declare_entry ("Stress exponents for dislocation creep", "3.5",
                             Patterns::List(Patterns::Double(0.)),
                             "List of stress exponents, $n_{\\text{dislocation}}$, for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value.  Units: None.");
          prm.declare_entry ("Activation energies for dislocation creep", "530e3",
                             Patterns::List(Patterns::Double(0.)),
                             "List of activation energies, $E_a$, for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value. "
                             "Units: \\si{\\joule\\per\\mole}.");
          prm.declare_entry ("Activation volumes for dislocation creep", "1.4e-5",
                             Patterns::List(Patterns::Double(0.)),
                             "List of activation volumes, $V_a$, for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value. "
                             "Units: \\si{\\meter\\cubed\\per\\mole}.");

          // Diffusion creep parameters
          Rheology::DiffusionCreep<dim>::declare_parameters(prm);

          // Dislocation creep parameters
          Rheology::DislocationCreep<dim>::declare_parameters(prm);

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    DiffusionDislocation<dim>::parse_parameters (ParameterHandler &prm)
    {
      // increment by one for background:
      const unsigned int n_fields = this->n_compositional_fields() + 1;

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection ("Diffusion dislocation");
        {
          // Initialise empty vector for compositional field variables
          std::vector<double> x_values;

          // Reference and minimum/maximum values
          reference_T = prm.get_double("Reference temperature");
          min_strain_rate = prm.get_double("Minimum strain rate");
          min_visc = prm.get_double ("Minimum viscosity");
          max_visc = prm.get_double ("Maximum viscosity");
          veff_coefficient = prm.get_double ("Effective viscosity coefficient");
          ref_visc = prm.get_double ("Reference viscosity");

          // Iteration parameters
          strain_rate_residual_threshold = prm.get_double ("Strain rate residual tolerance");
          stress_max_iteration_number = prm.get_integer ("Maximum strain rate ratio iterations");

          // Equation of state parameters
          thermal_diffusivity = prm.get_double("Thermal diffusivity");
          heat_capacity = prm.get_double("Heat capacity");

          // ---- Compositional parameters
          grain_size = prm.get_double("Grain size");
          densities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Densities"))),
                                                              n_fields,
                                                              "Densities");
          thermal_expansivities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Thermal expansivities"))),
                                                                          n_fields,
                                                                          "Thermal expansivities");

          viscosity_averaging = MaterialUtilities::parse_compositional_averaging_operation ("Viscosity averaging scheme",
                                prm);

          // Rheological parameters
          // Diffusion creep parameters
          diffusion_creep.initialize_simulator (this->get_simulator());
          diffusion_creep.parse_parameters(prm, std::make_shared<std::vector<unsigned int>>(n_fields));

          // Dislocation creep parameters
          dislocation_creep.initialize_simulator (this->get_simulator());
          dislocation_creep.parse_parameters(prm, std::make_shared<std::vector<unsigned int>>(n_fields));


        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Declare dependencies on solution variables
      this->model_dependence.viscosity = NonlinearDependence::temperature | NonlinearDependence::pressure | NonlinearDependence::strain_rate | NonlinearDependence::compositional_fields;
      this->model_dependence.density = NonlinearDependence::temperature | NonlinearDependence::pressure | NonlinearDependence::compositional_fields;
      this->model_dependence.compressibility = NonlinearDependence::none;
      this->model_dependence.specific_heat = NonlinearDependence::none;
      this->model_dependence.thermal_conductivity = NonlinearDependence::temperature | NonlinearDependence::pressure | NonlinearDependence::compositional_fields;
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(DiffusionDislocation,
                                   "diffusion dislocation",
                                   "An implementation of a viscous rheology including diffusion "
                                   "and dislocation creep. "
                                   "Compositional fields can each be assigned individual "
                                   "activation energies, reference densities, thermal expansivities, "
                                   "and stress exponents. The effective viscosity is defined as "
                                   "\n\n"
                                   "\\[\\eta_{\\text{eff}} = \\left(\\frac{1}{\\eta_{\\text{eff}}^\\text{diff}}+ "
                                   "\\frac{1}{\\eta_{\\text{eff}}^\\text{dis}}\\right)^{-1}\\] "
                                   "where "
                                   "\\[\\eta_{\\text{i}} = \\frac{1}{2} A^{-\\frac{1}{n_i}} d^\\frac{m_i}{n_i} "
                                   "\\dot{\\varepsilon_i}^{\\frac{1-n_i}{n_i}} "
                                   "\\exp\\left(\\frac{E_i^\\ast + PV_i^\\ast}{n_iRT}\\right)\\] "
                                   "\n\n"
                                   "where $d$ is grain size, $i$ corresponds to diffusion or dislocation creep, "
                                   "$\\dot{\\varepsilon}$ is the square root of the second invariant of the "
                                   "strain rate tensor, $R$ is the gas constant, $T$ is temperature, "
                                   "and $P$ is pressure. "
                                   "$A_i$ are prefactors, $n_i$ and $m_i$ are stress and grain size exponents "
                                   "$E_i$ are the activation energies and $V_i$ are the activation volumes. "
                                   "\n\n"
                                   "This form of the viscosity equation is commonly used in geodynamic simulations "
                                   "See, for example, Billen and Hirth (2007), G3, 8, Q08012. Significantly, "
                                   "other studies may use slightly different forms of the viscosity equation "
                                   "leading to variations in how specific terms are defined or combined. For "
                                   "example, the grain size exponent should always be positive in the diffusion "
                                   "viscosity equation used here, while other studies place the grain size term "
                                   "in the denominator and invert the sign of the grain size exponent. When "
                                   "examining previous work, one should carefully check how the viscous "
                                   "prefactor and grain size terms are defined. "
                                   " \n\n"
                                   "The ratio of diffusion to dislocation strain rate is found by Newton's "
                                   "method, iterating to find the stress which satisfies the above equations. "
                                   "The value for the components of this formula and additional "
                                   "parameters are read from the parameter file in subsection "
                                   "'Material model/DiffusionDislocation'.")
  }
}
