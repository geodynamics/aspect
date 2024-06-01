/*
  Copyright (C) 2020 - 2024 by the authors of the ASPECT code.

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


#include <aspect/material_model/rheology/diffusion_dislocation.h>
#include <aspect/material_model/utilities.h>
#include <aspect/utilities.h>

#include <deal.II/sundials/kinsol.h>
#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/parameter_handler.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace Rheology
    {
      template <int dim>
      DiffusionDislocation<dim>::DiffusionDislocation ()
        = default;

      template <int dim>
      std::vector<double>
      DiffusionDislocation<dim>::
      calculate_isostrain_viscosities ( const double pressure,
                                        const double temperature,
                                        const SymmetricTensor<2,dim> &strain_rate) const
      {
        // This function calculates viscosities assuming that all the compositional fields
        // experience the same strain rate (isostrain).

        // If strain rate is zero (like during the first time step) set it to some very small number
        // to prevent a division-by-zero, and a floating point exception.
        // Otherwise, calculate the square-root of the norm of the second invariant of the deviatoric-
        // strain rate (often simplified as epsilondot_ii)
        const double edot_ii = std::max(std::sqrt(std::max(-second_invariant(deviator(strain_rate)), 0.)),
                                        min_strain_rate);
        const double log_edot_ii = std::log(edot_ii);
        double log_strain_rate_deriv;

        // Find effective viscosities for each of the individual phases
        // Viscosities should have the same number of entries as compositional fields
        std::vector<double> composition_viscosities(n_chemical_composition_fields);
        for (unsigned int j=0; j < n_chemical_composition_fields; ++j)
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
            // We use KINSOL to find the second invariant of the stress tensor.
            // Start with the assumption that all strain is accommodated by diffusion creep:
            // If the diffusion creep prefactor is very small, that means that the diffusion viscosity is very large.
            // In this case, use the maximum viscosity instead to compute the starting guess.
            double stress_ii = (prefactor_stress_diffusion > (0.5 / maximum_viscosity)
                                ?
                                edot_ii/prefactor_stress_diffusion
                                :
                                0.5 / maximum_viscosity);
            Vector<double> log_stress_ii(1);
            log_stress_ii[0] = std::log(stress_ii);

            typename SUNDIALS::KINSOL<Vector<double>>::AdditionalData additional_data;
            additional_data.function_tolerance = log_strain_rate_residual_threshold;

            SUNDIALS::KINSOL<Vector<double>> nonlinear_solver(additional_data);

            nonlinear_solver.reinit_vector = [&](Vector<double> &x)
            {
              x.reinit(1);
            };

            nonlinear_solver.residual =
              [&](const Vector<double> &evaluation_point,
                  Vector<double>       &residual)
            {
              compute_log_strain_rate_residual(evaluation_point, residual, pressure, temperature, diffusion_creep_parameters, dislocation_creep_parameters, log_edot_ii);
            };

            nonlinear_solver.setup_jacobian =
              [&](const Vector<double> &current_u,
                  const Vector<double> & /*current_f*/)
            {
              compute_log_strain_rate_deriv(current_u, pressure, temperature, diffusion_creep_parameters, dislocation_creep_parameters, log_strain_rate_deriv);
            };

            nonlinear_solver.solve_with_jacobian = [&](const Vector<double> &rhs,
                                                       Vector<double> &solution,
                                                       const double /*tolerance*/)
            {
              solution(0) -= rhs(0)/log_strain_rate_deriv;
            };

            nonlinear_solver.solve(log_stress_ii);

            stress_ii = std::exp(log_stress_ii[0]);
            composition_viscosities[j] = std::min(std::max(stress_ii/edot_ii/2, minimum_viscosity), maximum_viscosity);
          }
        return composition_viscosities;
      }


      template <int dim>
      void
      DiffusionDislocation<dim>::compute_log_strain_rate_deriv(
        const Vector<double> &evaluation_point,
        double pressure,
        double temperature,
        const Rheology::DiffusionCreepParameters diffusion_creep_parameters,
        const Rheology::DislocationCreepParameters dislocation_creep_parameters,
        double &log_strain_rate_deriv) const
      {
        const std::pair<double, double> log_diff_edot_and_deriv = diffusion_creep.compute_log_strain_rate_and_derivative(evaluation_point[0], pressure, temperature, diffusion_creep_parameters);
        const std::pair<double, double> log_disl_edot_and_deriv = dislocation_creep.compute_log_strain_rate_and_derivative(evaluation_point[0], pressure, temperature, dislocation_creep_parameters);

        const double strain_rate_diffusion = std::exp(log_diff_edot_and_deriv.first);
        const double strain_rate_dislocation = std::exp(log_disl_edot_and_deriv.first);
        log_strain_rate_deriv = (strain_rate_diffusion * log_diff_edot_and_deriv.second + strain_rate_dislocation * log_disl_edot_and_deriv.second)/
                                (strain_rate_diffusion + strain_rate_dislocation);
      }



      template <int dim>
      void
      DiffusionDislocation<dim>::compute_log_strain_rate_residual(
        const Vector<double> &evaluation_point,
        Vector<double>       &residual,
        double pressure,
        double temperature,
        const Rheology::DiffusionCreepParameters diffusion_creep_parameters,
        const Rheology::DislocationCreepParameters dislocation_creep_parameters,
        double log_edot_ii) const
      {
        const std::pair<double, double> log_diff_edot_and_deriv = diffusion_creep.compute_log_strain_rate_and_derivative(evaluation_point[0], pressure, temperature, diffusion_creep_parameters);
        const std::pair<double, double> log_disl_edot_and_deriv = dislocation_creep.compute_log_strain_rate_and_derivative(evaluation_point[0], pressure, temperature, dislocation_creep_parameters);

        const double strain_rate_diffusion = std::exp(log_diff_edot_and_deriv.first);
        const double strain_rate_dislocation = std::exp(log_disl_edot_and_deriv.first);
        const double log_strain_rate_iterate = std::log(strain_rate_diffusion + strain_rate_dislocation);
        residual(0) = log_edot_ii - log_strain_rate_iterate;
      }

      template <int dim>
      double
      DiffusionDislocation<dim>::compute_viscosity (const double pressure,
                                                    const double temperature,
                                                    const std::vector<double> &volume_fractions,
                                                    const SymmetricTensor<2,dim> &strain_rate) const
      {
        // Currently, the viscosities for each of the compositional fields are calculated assuming
        // isostrain amongst all compositions, allowing calculation of the viscosity ratio.
        // TODO: This is only consistent with viscosity averaging if the arithmetic averaging
        // scheme is chosen. It would be useful to have a function to calculate isostress viscosities.
        const std::vector<double> composition_viscosities =
          calculate_isostrain_viscosities(pressure, temperature, strain_rate);

        // The isostrain condition implies that the viscosity averaging should be arithmetic (see above).
        // We have given the user freedom to apply alternative bounds, because in diffusion-dominated
        // creep (where n_diff=1) viscosities are stress and strain-rate independent, so the calculation
        // of compositional field viscosities is consistent with any averaging scheme.
        return MaterialUtilities::average_value(volume_fractions, composition_viscosities, viscosity_averaging);
      }


      template <int dim>
      void
      DiffusionDislocation<dim>::declare_parameters (ParameterHandler &prm)
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

        // Viscosity iteration parameters
        prm.declare_entry ("Strain rate residual tolerance", "1e-10", Patterns::Double(0.),
                           "Tolerance for determining the correct stress and viscosity from the "
                           "strain rate by internal iteration. The tolerance is expressed as the "
                           "difference between the natural logarithm of the input strain rate and "
                           "the strain rate at the current iteration. This determines that strain "
                           "rate is correctly partitioned between diffusion and dislocation creep "
                           "assuming that both mechanisms experience the same stress.");
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
                           "for a total of N+1 values, where N is the number of all compositional fields or only "
                           "those corresponding to chemical compositions. "
                           "If only one value is given, then all use the same value. "
                           "Units: \\si{\\kilogram\\per\\meter\\cubed}.");
        prm.declare_entry ("Thermal expansivities", "3.5e-5",
                           Patterns::List(Patterns::Double(0.)),
                           "List of thermal expansivities for background mantle and compositional fields, "
                           "for a total of N+1 values, where N is the number of all compositional fields or only "
                           "those corresponding to chemical compositions. "
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
                           "for a total of N+1 values, where N is the number of all compositional fields or only "
                           "those corresponding to chemical compositions. "
                           "If only one value is given, then all use the same value. "
                           "Units: \\si{\\per\\pascal} \\si{\\meter}$^{m_{\\text{diffusion}}}$ \\si{\\per\\second}.");
        prm.declare_entry ("Stress exponents for diffusion creep", "1.",
                           Patterns::List(Patterns::Double(0.)),
                           "List of stress exponents, $n_{\\text{diffusion}}$, for background mantle and compositional fields, "
                           "for a total of N+1 values, where N is the number of all compositional fields or only "
                           "those corresponding to chemical compositions. "
                           "If only one value is given, then all use the same value.  Units: None.");
        prm.declare_entry ("Grain size exponents for diffusion creep", "3.",
                           Patterns::List(Patterns::Double(0.)),
                           "List of grain size exponents, $m_{\\text{diffusion}}$, for background mantle and compositional fields, "
                           "for a total of N+1 values, where N is the number of all compositional fields or only "
                           "those corresponding to chemical compositions. "
                           "If only one value is given, then all use the same value.  Units: None.");
        prm.declare_entry ("Activation energies for diffusion creep", "375e3",
                           Patterns::List(Patterns::Double(0.)),
                           "List of activation energies, $E_a$, for background mantle and compositional fields, "
                           "for a total of N+1 values, where N is the number of all compositional fields or only "
                           "those corresponding to chemical compositions. "
                           "If only one value is given, then all use the same value. "
                           "Units: \\si{\\joule\\per\\mole}.");
        prm.declare_entry ("Activation volumes for diffusion creep", "6e-6",
                           Patterns::List(Patterns::Double(0.)),
                           "List of activation volumes, $V_a$, for background mantle and compositional fields, "
                           "for a total of N+1 values, where N is the number of all compositional fields or only "
                           "those corresponding to chemical compositions. "
                           "If only one value is given, then all use the same value. "
                           "Units: \\si{\\meter\\cubed\\per\\mole}.");

        // Dislocation creep parameters
        prm.declare_entry ("Prefactors for dislocation creep", "1.1e-16",
                           Patterns::List(Patterns::Double(0.)),
                           "List of viscosity prefactors, $A$, for background mantle and compositional fields, "
                           "for a total of N+1 values, where N is the number of all compositional fields or only "
                           "those corresponding to chemical compositions. "
                           "If only one value is given, then all use the same value. "
                           "Units: \\si{\\pascal}$^{-n_{\\text{dislocation}}}$\\si{\\per\\second}.");
        prm.declare_entry ("Stress exponents for dislocation creep", "3.5",
                           Patterns::List(Patterns::Double(0.)),
                           "List of stress exponents, $n_{\\text{dislocation}}$, for background mantle and compositional fields, "
                           "for a total of N+1 values, where N is the number of all compositional fields or only "
                           "those corresponding to chemical compositions. "
                           "If only one value is given, then all use the same value.  Units: None.");
        prm.declare_entry ("Activation energies for dislocation creep", "530e3",
                           Patterns::List(Patterns::Double(0.)),
                           "List of activation energies, $E_a$, for background mantle and compositional fields, "
                           "for a total of N+1 values, where N is the number of all compositional fields or only "
                           "those corresponding to chemical compositions. "
                           "If only one value is given, then all use the same value. "
                           "Units: \\si{\\joule\\per\\mole}.");
        prm.declare_entry ("Activation volumes for dislocation creep", "1.4e-5",
                           Patterns::List(Patterns::Double(0.)),
                           "List of activation volumes, $V_a$, for background mantle and compositional fields, "
                           "for a total of N+1 values, where N is the number of all compositional fields or only "
                           "those corresponding to chemical compositions. "
                           "If only one value is given, then all use the same value. "
                           "Units: \\si{\\meter\\cubed\\per\\mole}.");

        // Diffusion creep parameters
        Rheology::DiffusionCreep<dim>::declare_parameters(prm);

        // Dislocation creep parameters
        Rheology::DislocationCreep<dim>::declare_parameters(prm);
      }



      template <int dim>
      void
      DiffusionDislocation<dim>::parse_parameters (ParameterHandler &prm)
      {
        // Reference and minimum/maximum values
        min_strain_rate = prm.get_double("Minimum strain rate");
        minimum_viscosity = prm.get_double ("Minimum viscosity");
        maximum_viscosity = prm.get_double ("Maximum viscosity");

        // Iteration parameters
        log_strain_rate_residual_threshold = prm.get_double ("Strain rate residual tolerance");
        stress_max_iteration_number = prm.get_integer ("Maximum strain rate ratio iterations");

        // ---- Compositional parameters
        grain_size = prm.get_double("Grain size");

        viscosity_averaging = MaterialUtilities::parse_compositional_averaging_operation ("Viscosity averaging scheme",
                              prm);

        // Rheological parameters for chemical compositions
        // increment by one for background:
        n_chemical_composition_fields = this->introspection().n_chemical_composition_fields() + 1;

        // Diffusion creep parameters
        diffusion_creep.initialize_simulator (this->get_simulator());
        diffusion_creep.parse_parameters(prm, std::make_unique<std::vector<unsigned int>>(n_chemical_composition_fields));

        // Dislocation creep parameters
        dislocation_creep.initialize_simulator (this->get_simulator());
        dislocation_creep.parse_parameters(prm, std::make_unique<std::vector<unsigned int>>(n_chemical_composition_fields));

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
  namespace Rheology \
  { \
    template class DiffusionDislocation<dim>; \
  }

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
