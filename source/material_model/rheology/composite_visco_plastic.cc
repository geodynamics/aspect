/*
  Copyright (C) 2020 by the authors of the ASPECT code.

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


#include <aspect/material_model/rheology/composite_visco_plastic.h>
#include <aspect/material_model/utilities.h>
#include <aspect/utilities.h>

#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/parameter_handler.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace Rheology
    {
      // The following function is not yet used by ASPECT, so
      // it is commented out to avoid tester issues.
      //
      // It will provide useful output when a material model which
      // uses the CompositeViscoPlastic rheology has been implemented.
      //
      // The composite visco plastic rheology calculates the decomposed strain
      // rates for each of the following deformation mechanisms:
      // diffusion creep, dislocation creep, Peierls creep,
      // Drucker-Prager plasticity and a constant (high) viscosity limiter.
      // The values are provided in this order as a vector of additional
      // outputs. If the user declares one or more mechanisms inactive
      // (by assigning use_mechanism = False) then the corresponding
      // strain rate output will be equal to zero.
      //namespace
      //{
      //  std::vector<std::string> make_strain_rate_additional_outputs_names()
      //  {
      //    std::vector<std::string> names;
      //    names.emplace_back("edot_diffusion");
      //    names.emplace_back("edot_dislocation");
      //    names.emplace_back("edot_peierls");
      //    names.emplace_back("edot_drucker_prager");
      //    names.emplace_back("edot_limiter");
      //    return names;
      //  }
      //}


      template <int dim>
      CompositeViscoPlastic<dim>::CompositeViscoPlastic ()
      {}



      template <int dim>
      double
      CompositeViscoPlastic<dim>::compute_viscosity (const double pressure,
                                                     const double temperature,
                                                     const unsigned int composition,
                                                     const SymmetricTensor<2,dim> &strain_rate,
                                                     std::vector<double> &partial_strain_rates,
                                                     const std::vector<double> &phase_function_values,
                                                     const std::vector<unsigned int> &n_phases_per_composition) const
      {
        // If strain rate is zero (like during the first time step) set it to some very small number
        // to prevent a division-by-zero, and a floating point exception.
        // Otherwise, calculate the square-root of the norm of the second invariant of the deviatoric-
        // strain rate (often simplified as epsilondot_ii)
        const double edot_ii = std::max(std::sqrt(std::fabs(second_invariant(deviator(strain_rate)))),
                                        min_strain_rate);

        Rheology::DiffusionCreepParameters diffusion_creep_parameters;
        Rheology::DislocationCreepParameters dislocation_creep_parameters;
        Rheology::PeierlsCreepParameters peierls_creep_parameters;
        double eta_diff = max_viscosities[composition];
        double eta_disl = max_viscosities[composition];
        double eta_prls = max_viscosities[composition];


        if (use_diffusion_creep)
          {
            diffusion_creep_parameters = diffusion_creep->compute_creep_parameters(composition, phase_function_values, n_phases_per_composition);
            eta_diff = diffusion_creep->compute_viscosity(pressure, temperature, composition, phase_function_values, n_phases_per_composition);
          }

        if (use_dislocation_creep)
          {
            dislocation_creep_parameters = dislocation_creep->compute_creep_parameters(composition, phase_function_values, n_phases_per_composition);
            eta_disl = dislocation_creep->compute_viscosity(edot_ii, pressure, temperature, composition, phase_function_values, n_phases_per_composition);
          }

        if (use_peierls_creep)
          {
            peierls_creep_parameters = peierls_creep->compute_creep_parameters(composition);
            eta_prls = peierls_creep->compute_approximate_viscosity(edot_ii, pressure, temperature, composition);
          }
        // First guess at a stress using diffusion, dislocation, and Peierls creep viscosities calculated with the total second strain rate invariant.
        const double eta_guess = std::min(std::max(min_viscosities[composition], eta_diff*eta_disl*eta_prls/(eta_diff*eta_disl + eta_diff*eta_prls + eta_disl*eta_prls)), max_viscosities[composition]);

        double creep_stress = 2.*eta_guess*edot_ii;

        // Crude modification of the creep stress to be no higher than the
        // Drucker-Prager yield stress. Probably fine for a first guess.
        if (use_drucker_prager)
          {
            const double yield_stress = drucker_prager->compute_yield_stress(drucker_prager_parameters.cohesions[composition],
                                                                             drucker_prager_parameters.angles_internal_friction[composition],
                                                                             pressure,
                                                                             drucker_prager_parameters.max_yield_stress);
            creep_stress = std::min(creep_stress, yield_stress);
          }

        // In this rheology model, the total strain rate is partitioned between
        // different flow laws. We do not know how the strain is partitioned
        // between these flow laws, nor do we know the creep stress, which is
        // required to calculate the partitioning.

        // The following while loop conducts a Newton iteration to obtain the
        // creep stress, which we need in order to calculate the viscosity.
        double strain_rate_residual = 2*strain_rate_residual_threshold;
        double strain_rate_deriv = 0;
        unsigned int stress_iteration = 0;
        while (std::abs(strain_rate_residual) > strain_rate_residual_threshold
               && stress_iteration < stress_max_iteration_number)
          {

            const std::pair<double, double> edot_and_deriv = compute_strain_rate_and_derivative (creep_stress,
                                                             pressure,
                                                             temperature,
                                                             composition,
                                                             diffusion_creep_parameters,
                                                             dislocation_creep_parameters,
                                                             peierls_creep_parameters,
                                                             drucker_prager_parameters,
                                                             min_viscosities[composition],
                                                             max_viscosities[composition]);

            strain_rate_residual = edot_and_deriv.first - edot_ii;
            strain_rate_deriv = edot_and_deriv.second;

            // If the strain rate derivative is zero, we catch it below.
            if (strain_rate_deriv>std::numeric_limits<double>::min())
              creep_stress -= strain_rate_residual/strain_rate_deriv;
            stress_iteration += 1;

            // If anything that would be used in the next iteration is not finite, the
            // Newton iteration would trigger an exception and we want to abort the
            // iteration instead.
            // Currently, we still throw an exception, but if this exception is thrown,
            // another more robust iterative scheme should be implemented
            // (similar to that seen in the diffusion-dislocation material model).
            const bool abort_newton_iteration = !numbers::is_finite(creep_stress)
                                                || !numbers::is_finite(strain_rate_residual)
                                                || !numbers::is_finite(strain_rate_deriv)
                                                || strain_rate_deriv < std::numeric_limits<double>::min()
                                                || stress_iteration == stress_max_iteration_number;
            AssertThrow(!abort_newton_iteration,
                        ExcMessage("No convergence has been reached in the loop that determines "
                                   "the composite viscous creep stress. Aborting! "
                                   "Residual is " + Utilities::to_string(strain_rate_residual) +
                                   " after " + Utilities::to_string(stress_iteration) + " iterations. "
                                   "You can increase the number of iterations by adapting the "
                                   "parameter 'Maximum creep strain rate iterations'."));
          }

        // The creep stress is not the total stress, so we still need to do a little work to obtain the effective viscosity.
        // First, we compute the stress running through the strain rate limiter, and then add that to the creep stress
        // NOTE: The viscosity of the strain rate limiter is equal to (min_visc*max_visc)/(max_visc - min_visc)
        const double lim_stress = 2.*min_viscosities[composition]*(edot_ii - creep_stress/(2.*max_viscosities[composition]));
        const double total_stress = creep_stress + lim_stress;

        // Compute the strain rate experienced by the different mechanisms
        // These should sum to the total strain rate
        partial_strain_rates.resize(5, 0.);

        // The components of partial_strain_rates must be provided in the order
        // dictated by make_strain_rate_additional_outputs_names
        if (use_diffusion_creep)
          {
            const std::pair<double, double> diff_edot_and_deriv = diffusion_creep->compute_strain_rate_and_derivative(creep_stress, pressure, temperature, diffusion_creep_parameters);
            partial_strain_rates[0] = diff_edot_and_deriv.first;
          }

        if (use_dislocation_creep)
          {
            const std::pair<double, double> disl_edot_and_deriv = dislocation_creep->compute_strain_rate_and_derivative(creep_stress, pressure, temperature, dislocation_creep_parameters);
            partial_strain_rates[1] = disl_edot_and_deriv.first;
          }

        if (use_peierls_creep)
          {
            const std::pair<double, double> prls_edot_and_deriv = peierls_creep->compute_strain_rate_and_derivative(creep_stress, pressure, temperature, peierls_creep_parameters);
            partial_strain_rates[2] = prls_edot_and_deriv.first;
          }

        if (use_drucker_prager)
          {
            const std::pair<double, double> drpr_edot_and_deriv = drucker_prager->compute_strain_rate_and_derivative(creep_stress, pressure, composition, drucker_prager_parameters);
            partial_strain_rates[3] = drpr_edot_and_deriv.first;
          }

        partial_strain_rates[4] = total_stress/(2.*max_viscosities[composition]);

        // Now we return the viscosity using the total stress
        return total_stress/(2.*edot_ii);
      }



      // Overload the + operator to act on two pairs of doubles.
      std::pair<double,double> operator+(const std::pair<double,double> &x, const std::pair<double,double> &y)
      {
        return std::make_pair(x.first+y.first, x.second+y.second);
      }



      template <int dim>
      std::pair<double, double>
      CompositeViscoPlastic<dim>::compute_strain_rate_and_derivative (const double creep_stress,
                                                                      const double pressure,
                                                                      const double temperature,
                                                                      const unsigned int composition,
                                                                      const DiffusionCreepParameters diffusion_creep_parameters,
                                                                      const DislocationCreepParameters dislocation_creep_parameters,
                                                                      const PeierlsCreepParameters peierls_creep_parameters,
                                                                      const DruckerPragerParameters drucker_prager_parameters,
                                                                      const double min_viscosity,
                                                                      const double max_viscosity) const
      {
        std::pair<double, double> creep_edot_and_deriv = std::make_pair(0., 0.);

        if (use_diffusion_creep)
          creep_edot_and_deriv = creep_edot_and_deriv + diffusion_creep->compute_strain_rate_and_derivative(creep_stress, pressure, temperature, diffusion_creep_parameters);

        if (use_dislocation_creep)
          creep_edot_and_deriv = creep_edot_and_deriv + dislocation_creep->compute_strain_rate_and_derivative(creep_stress, pressure, temperature, dislocation_creep_parameters);

        if (use_peierls_creep)
          creep_edot_and_deriv = creep_edot_and_deriv + peierls_creep->compute_strain_rate_and_derivative(creep_stress, pressure, temperature, peierls_creep_parameters);

        if (use_drucker_prager)
          creep_edot_and_deriv = creep_edot_and_deriv + drucker_prager->compute_strain_rate_and_derivative(creep_stress, pressure, composition, drucker_prager_parameters);

        const double strain_rate = creep_stress/(2.*max_viscosity) + (max_viscosity/(max_viscosity - min_viscosity))*creep_edot_and_deriv.first;
        const double strain_rate_deriv = 1./(2.*max_viscosity) + (max_viscosity/(max_viscosity - min_viscosity))*creep_edot_and_deriv.second;

        return std::make_pair(strain_rate, strain_rate_deriv);
      }



      template <int dim>
      void
      CompositeViscoPlastic<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Include diffusion creep", "true",
                           Patterns::Bool (),
                           "Whether to include diffusion creep in the rheological formulation.");

        prm.declare_entry ("Include dislocation creep", "true",
                           Patterns::Bool (),
                           "Whether to include dislocation creep in the rheological formulation.");

        prm.declare_entry ("Include Peierls creep", "true",
                           Patterns::Bool (),
                           "Whether to include Peierls creep in the rheological formulation.");

        prm.declare_entry ("Include Drucker Prager plasticity", "true",
                           Patterns::Bool (),
                           "Whether to include Drucker Prager plasticity in the rheological formulation.");


        // Diffusion creep parameters
        Rheology::DiffusionCreep<dim>::declare_parameters(prm);

        // Dislocation creep parameters
        Rheology::DislocationCreep<dim>::declare_parameters(prm);

        // Dislocation creep parameters
        Rheology::PeierlsCreep<dim>::declare_parameters(prm);

        // Drucker Prager parameters
        Rheology::DruckerPrager<dim>::declare_parameters(prm);

        // Some of the parameters below are shared with the subordinate
        // rheology models (diffusion, dislocation, ...),
        // and will already have been declared. This is fine, the deal.II
        // parameter handler allows repeated declarations. The default values
        // below will override any defaults given in the subordinate
        // rheology models. The new defaults will apply both to
        // this rheology model and to the subordinate rheology modules.
        prm.declare_entry ("Minimum strain rate", "1.4e-20", Patterns::Double(0.),
                           "Stabilizes strain dependent viscosity. Units: \\si{\\per\\second}.");

        // Viscosity iteration parameters
        prm.declare_entry ("Strain rate residual tolerance", "1e-22", Patterns::Double(0.),
                           "Tolerance for correct diffusion/dislocation strain rate ratio.");
        prm.declare_entry ("Maximum creep strain rate iterations", "40", Patterns::Integer(0),
                           "Maximum number of iterations to find the correct "
                           "creep strain rate.");

        // Strain rate and stress limiting parameters
        prm.declare_entry ("Minimum viscosities", "1.e17",
                           Patterns::List(Patterns::Double(0.)),
                           "List of minimum effective viscosities. Units: \\si{\\pascal\\second}.");

        prm.declare_entry ("Maximum viscosities", "1.e28",
                           Patterns::List(Patterns::Double(0.)),
                           "List of maximum effective viscosities. Units: \\si{\\pascal\\second}.");
      }



      template <int dim>
      void
      CompositeViscoPlastic<dim>::parse_parameters (ParameterHandler &prm,
                                                    const std::shared_ptr<std::vector<unsigned int>> &expected_n_phases_per_composition)
      {
        // Retrieve the list of composition names
        const std::vector<std::string> list_of_composition_names = this->introspection().get_composition_names();

        // Establish that a background field is required here
        const bool has_background_field = true;

        min_strain_rate = prm.get_double("Minimum strain rate");

        // Iteration parameters
        strain_rate_residual_threshold = prm.get_double ("Strain rate residual tolerance");
        stress_max_iteration_number = prm.get_integer ("Maximum creep strain rate iterations");

        // Read min and max viscosity parameters, each of size of number of composition + number of phases + 1
        min_viscosities = Utilities::parse_map_to_double_array(prm.get("Minimum viscosities"),
                                                               list_of_composition_names,
                                                               has_background_field,
                                                               "Minimum viscosities",
                                                               true,
                                                               expected_n_phases_per_composition);
        max_viscosities = Utilities::parse_map_to_double_array(prm.get("Maximum viscosities"),
                                                               list_of_composition_names,
                                                               has_background_field,
                                                               "Maximum viscosities",
                                                               true,
                                                               expected_n_phases_per_composition);

        // Rheological parameters

        // Diffusion creep parameters
        use_diffusion_creep = prm.get_bool ("Include diffusion creep");
        if (use_diffusion_creep)
          {
            diffusion_creep = std_cxx14::make_unique<Rheology::DiffusionCreep<dim>>();
            diffusion_creep->initialize_simulator (this->get_simulator());
            diffusion_creep->parse_parameters(prm, expected_n_phases_per_composition);
          }

        // Dislocation creep parameters
        use_dislocation_creep = prm.get_bool ("Include dislocation creep");
        if (use_dislocation_creep)
          {
            dislocation_creep = std_cxx14::make_unique<Rheology::DislocationCreep<dim>>();
            dislocation_creep->initialize_simulator (this->get_simulator());
            dislocation_creep->parse_parameters(prm, expected_n_phases_per_composition);
          }

        // Peierls creep parameters
        use_peierls_creep = prm.get_bool ("Include Peierls creep");
        if (use_peierls_creep)
          {
            peierls_creep = std_cxx14::make_unique<Rheology::PeierlsCreep<dim>>();
            peierls_creep->initialize_simulator (this->get_simulator());
            peierls_creep->parse_parameters(prm);
          }

        // Drucker Prager parameters
        use_drucker_prager = prm.get_bool ("Include Drucker Prager plasticity");
        if (use_drucker_prager)
          {
            drucker_prager = std_cxx14::make_unique<Rheology::DruckerPrager<dim>>();
            drucker_prager_parameters = drucker_prager->parse_parameters(list_of_composition_names.size() + 1, prm);

            AssertThrow(drucker_prager_parameters.use_plastic_damper == true && drucker_prager_parameters.damper_viscosity > 0.,
                        ExcMessage("If Drucker Prager is included in the rheological formulation, you must use a viscous damper with a positive viscosity."));
          }

        AssertThrow(use_diffusion_creep == true || use_dislocation_creep == true || use_peierls_creep == true || use_drucker_prager == true,
                    ExcMessage("You need to include at least one deformation mechanism."));

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
    template class CompositeViscoPlastic<dim>; \
  }

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
