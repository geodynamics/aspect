/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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

#include <aspect/material_model/visco_plastic.h>
#include <aspect/utilities.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/signaling_nan.h>
#include <aspect/newton.h>
#include <aspect/adiabatic_conditions/interface.h>

namespace aspect
{
  namespace MaterialModel
  {

    namespace
    {
      std::vector<std::string> make_plastic_additional_outputs_names()
      {
        std::vector<std::string> names;
        names.emplace_back("current_cohesions");
        names.emplace_back("current_friction_angles");
        names.emplace_back("plastic_yielding");
        return names;
      }
    }

    template <int dim>
    bool
    ViscoPlastic<dim>::
    is_yielding (const double &pressure,
                 const double &temperature,
                 const std::vector<double> &composition,
                 const SymmetricTensor<2,dim> &strain_rate) const
    {
      /* The following returns whether or not the material is plastically yielding
       * as documented in evaluate.
       */
      bool plastic_yielding = false;

      const std::vector<double> volume_fractions = MaterialUtilities::compute_volume_fractions(composition, get_volumetric_composition_mask());

      const std::pair<std::vector<double>, std::vector<bool> > calculate_viscosities =
        calculate_isostrain_viscosities(volume_fractions, pressure, temperature, composition, strain_rate, viscous_flow_law, yield_mechanism);

      std::vector<double>::const_iterator max_composition = std::max_element(volume_fractions.begin(),volume_fractions.end());
      plastic_yielding = calculate_viscosities.second[std::distance(volume_fractions.begin(),max_composition)];

      return plastic_yielding;
    }

    template <int dim>
    PlasticAdditionalOutputs<dim>::PlasticAdditionalOutputs (const unsigned int n_points)
      :
      NamedAdditionalMaterialOutputs<dim>(make_plastic_additional_outputs_names()),
      cohesions(n_points, numbers::signaling_nan<double>()),
      friction_angles(n_points, numbers::signaling_nan<double>()),
      yielding(n_points, numbers::signaling_nan<double>())
    {}

    template <int dim>
    std::vector<double>
    PlasticAdditionalOutputs<dim>::get_nth_output(const unsigned int idx) const
    {
      AssertIndexRange (idx, 3);
      switch (idx)
        {
          case 0:
            return cohesions;

          case 1:
            return friction_angles;

          case 2:
            return yielding;

          default:
            AssertThrow(false, ExcInternalError());
        }
      // We will never get here, so just return something
      return cohesions;
    }



    template <int dim>
    std::pair<std::vector<double>, std::vector<bool> >
    ViscoPlastic<dim>::
    calculate_isostrain_viscosities (const std::vector<double> &volume_fractions,
                                     const double &pressure,
                                     const double &temperature,
                                     const std::vector<double> &composition,
                                     const SymmetricTensor<2,dim> &strain_rate,
                                     const ViscosityScheme &viscous_type,
                                     const YieldScheme &yield_type) const
    {
      // This function calculates viscosities assuming that all the compositional fields
      // experience the same strain rate (isostrain).

      // Calculate the square root of the second moment invariant for the deviatoric strain rate tensor.
      // The first time this function is called (first iteration of first time step)
      // a specified "reference" strain rate is used as the returned value would
      // otherwise be zero.
      const double edot_ii = ( (this->get_timestep_number() == 0 && strain_rate.norm() <= std::numeric_limits<double>::min())
                               ?
                               ref_strain_rate
                               :
                               std::max(std::sqrt(std::fabs(second_invariant(deviator(strain_rate)))),
                                        min_strain_rate) );

      // Choice of activation volume depends on whether there is an adiabatic temperature
      // gradient used when calculating the viscosity. This allows the same activation volume
      // to be used in incompressible and compressible models.
      const double temperature_for_viscosity = temperature + adiabatic_temperature_gradient_for_viscosity*pressure;
      Assert(temperature_for_viscosity != 0, ExcMessage(
               "The temperature used in the calculation of the visco-plastic rheology is zero. "
               "This is not allowed, because this value is used to divide through. It is probably "
               "being caused by the temperature being zero somewhere in the model. The relevant "
               "values for debugging are: temperature (" + Utilities::to_string(temperature) +
               "), adiabatic_temperature_gradient_for_viscosity ("
               + Utilities::to_string(adiabatic_temperature_gradient_for_viscosity) + ") and pressure ("
               + Utilities::to_string(pressure) + ")."));


      // First step: viscous behavior
      // Calculate viscosities for each of the individual compositional phases
      std::vector<double> composition_viscosities(volume_fractions.size());
      std::vector<bool> composition_yielding(volume_fractions.size());
      for (unsigned int j=0; j < volume_fractions.size(); ++j)
        {
          // Power law creep equation
          //    viscosity = 0.5 * A^(-1/n) * edot_ii^((1-n)/n) * d^(m/n) * exp((E + P*V)/(nRT))
          // A: prefactor, edot_ii: square root of second invariant of deviatoric strain rate tensor,
          // d: grain size, m: grain size exponent, E: activation energy, P: pressure,
          // V; activation volume, n: stress exponent, R: gas constant, T: temperature.
          // Note: values of A, d, m, E, V and n are distinct for diffusion & dislocation creep

          // Diffusion creep: viscosity is grain size dependent (m!=0) and strain-rate independent (n=1)
          double viscosity_diffusion = diffusion_creep.compute_viscosity(pressure, temperature_for_viscosity, j);
          //double viscosity_diffusion = 0.5 / prefactors_diffusion[j] *
          //                             std::exp((activation_energies_diffusion[j] + pressure*activation_volumes_diffusion[j])/
          //                                      (constants::gas_constant*temperature_for_viscosity)) *
          //                             std::pow(grain_size, grain_size_exponents_diffusion[j]);

          // For dislocation creep, viscosity is grain size independent (m=0) and strain-rate dependent (n>1)
          double viscosity_dislocation = dislocation_creep.compute_viscosity(edot_ii, pressure, temperature_for_viscosity, j);
          //double viscosity_dislocation = 0.5 * std::pow(prefactors_dislocation[j],-1/stress_exponents_dislocation[j]) *
          //                               std::exp((activation_energies_dislocation[j] + pressure*activation_volumes_dislocation[j])/
          //                                       (constants::gas_constant*temperature_for_viscosity*stress_exponents_dislocation[j])) *
          //                               std::pow(edot_ii,((1. - stress_exponents_dislocation[j])/stress_exponents_dislocation[j]));

          // Select what form of viscosity to use (diffusion, dislocation or composite)
          double viscosity_pre_yield = 0.0;
          switch (viscous_type)
            {
              case diffusion:
              {
                viscosity_pre_yield = viscosity_diffusion;
                break;
              }
              case dislocation:
              {
                viscosity_pre_yield = viscosity_dislocation;
                break;
              }
              case composite:
              {
                viscosity_pre_yield = (viscosity_diffusion * viscosity_dislocation)/(viscosity_diffusion + viscosity_dislocation);
                break;
              }
              default:
              {
                AssertThrow(false, ExcNotImplemented());
                break;
              }
            }


          // Second step: strain weakening

          // Calculate the strain weakening factors for cohesion, friction and viscosity. If no brittle and/or viscous strain weakening is applied, the factors are 1.
          const std::array<double, 3> weakening_factors = strain_rheology.compute_strain_weakening_factors(j, composition);

          const double current_cohesion = cohesions[j] * weakening_factors[0];
          const double current_friction = angles_internal_friction[j] * weakening_factors[1];
          viscosity_pre_yield *= weakening_factors[2];

          // Weakened friction and cohesion values
          std::pair<double, double> yield_parameters (current_cohesion, current_friction);

          // Third step: plastic yielding

          // Calculate Drucker-Prager yield strength (i.e. yield stress)
          const MaterialUtilities::DruckerPragerInputs plastic_in(yield_parameters.first, yield_parameters.second, std::max(pressure,0.0), edot_ii, max_yield_strength);
          MaterialUtilities::DruckerPragerOutputs plastic_out;
          MaterialUtilities::compute_drucker_prager_yielding<dim> (plastic_in, plastic_out);

          // If the viscous stress is greater than the yield strength, indicate we are in the yielding regime.
          const double viscous_stress = 2. * viscosity_pre_yield * edot_ii;
          if (viscous_stress >= plastic_out.yield_strength)
            composition_yielding[j] = true;

          // Select if yield viscosity is based on Drucker Prager or stress limiter rheology
          double viscosity_yield = viscosity_pre_yield;
          switch (yield_type)
            {
              case stress_limiter:
              {
                const double viscosity_limiter = plastic_out.yield_strength / (2.0 * ref_strain_rate)
                                                 * std::pow((edot_ii/ref_strain_rate), 1./exponents_stress_limiter[j] - 1.0);
                viscosity_yield = 1. / ( 1./viscosity_limiter + 1./viscosity_pre_yield);
                break;
              }
              case drucker_prager:
              {
                // If the viscous stress is greater than the yield strength, rescale the viscosity back to yield surface
                if (viscous_stress >= plastic_out.yield_strength)
                  viscosity_yield = plastic_out.plastic_viscosity;
                break;
              }
              default:
              {
                AssertThrow(false, ExcNotImplemented());
                break;
              }
            }

          // Limit the viscosity with specified minimum and maximum bounds
          composition_viscosities[j] = std::min(std::max(viscosity_yield, min_visc), max_visc);

        }
      return std::make_pair (composition_viscosities, composition_yielding);
    }


    template <int dim>
    void
    ViscoPlastic<dim>::
    fill_plastic_outputs(const unsigned int i,
                         const std::vector<double> &volume_fractions,
                         const bool plastic_yielding,
                         const MaterialModel::MaterialModelInputs<dim> &in,
                         MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      PlasticAdditionalOutputs<dim> *plastic_out = out.template get_additional_output<PlasticAdditionalOutputs<dim> >();

      if (plastic_out != nullptr)
        {
          double C = 0.;
          double phi = 0.;

          // set to weakened values, or unweakened values when strain weakening is not used
          for (unsigned int j=0; j < volume_fractions.size(); ++j)
            {
              // Calculate the strain weakening factors and weakened values
              const std::array<double, 3> weakening_factors = strain_rheology.compute_strain_weakening_factors(j, in.composition[i]);
              C   += volume_fractions[j] * (cohesions[j] * weakening_factors[0]);
              phi += volume_fractions[j] * (angles_internal_friction[j] * weakening_factors[1]);
            }

          plastic_out->cohesions[i] = C;
          plastic_out->friction_angles[i] = phi * 180. / numbers::PI;
          plastic_out->yielding[i] = plastic_yielding ? 1 : 0;
        }
    }



    template <int dim>
    void
    ViscoPlastic<dim>::
    compute_viscosity_derivatives(const unsigned int i,
                                  const std::vector<double> &volume_fractions,
                                  const std::vector<double> &composition_viscosities,
                                  const MaterialModel::MaterialModelInputs<dim> &in,
                                  MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      MaterialModel::MaterialModelDerivatives<dim> *derivatives =
        out.template get_additional_output<MaterialModel::MaterialModelDerivatives<dim> >();

      if (derivatives != nullptr)
        {
          // compute derivatives if necessary
          std::vector<SymmetricTensor<2,dim> > composition_viscosities_derivatives(volume_fractions.size());
          std::vector<double> composition_dviscosities_dpressure(volume_fractions.size());

          const double finite_difference_accuracy = 1e-7;

          // For each independent component, compute the derivative.
          for (unsigned int component = 0; component < SymmetricTensor<2,dim>::n_independent_components; ++component)
            {
              const TableIndices<2> strain_rate_indices = SymmetricTensor<2,dim>::unrolled_to_component_indices (component);

              const SymmetricTensor<2,dim> strain_rate_difference = in.strain_rate[i]
                                                                    + std::max(std::fabs(in.strain_rate[i][strain_rate_indices]), min_strain_rate)
                                                                    * finite_difference_accuracy
                                                                    * Utilities::nth_basis_for_symmetric_tensors<dim>(component);
              std::vector<double> eta_component =
                calculate_isostrain_viscosities(volume_fractions, in.pressure[i],
                                                in.temperature[i], in.composition[i],
                                                strain_rate_difference,
                                                viscous_flow_law,yield_mechanism).first;

              // For each composition of the independent component, compute the derivative.
              for (unsigned int composition_index = 0; composition_index < eta_component.size(); ++composition_index)
                {
                  // compute the difference between the viscosity with and without the strain-rate difference.
                  double viscosity_derivative = eta_component[composition_index] - composition_viscosities[composition_index];
                  if (viscosity_derivative != 0)
                    {
                      // when the difference is non-zero, divide by the difference.
                      viscosity_derivative /= std::max(std::fabs(strain_rate_difference[strain_rate_indices]), min_strain_rate)
                                              * finite_difference_accuracy;
                    }
                  composition_viscosities_derivatives[composition_index][strain_rate_indices] = viscosity_derivative;
                }
            }

          /**
           * Now compute the derivative of the viscosity to the pressure
           */
          const double pressure_difference = in.pressure[i] + (std::fabs(in.pressure[i]) * finite_difference_accuracy);

          const std::vector<double> viscosity_difference =
            calculate_isostrain_viscosities(volume_fractions, pressure_difference,
                                            in.temperature[i], in.composition[i], in.strain_rate[i],
                                            viscous_flow_law, yield_mechanism).first;


          for (unsigned int composition_index = 0; composition_index < viscosity_difference.size(); ++composition_index)
            {
              double viscosity_derivative = viscosity_difference[composition_index] - composition_viscosities[composition_index];
              if (viscosity_difference[composition_index] != 0)
                {
                  if (in.pressure[i] != 0)
                    {
                      viscosity_derivative /= std::fabs(in.pressure[i]) * finite_difference_accuracy;
                    }
                  else
                    {
                      viscosity_derivative = 0;
                    }
                }
              composition_dviscosities_dpressure[composition_index] = viscosity_derivative;
            }

          double viscosity_averaging_p = 0; // Geometric
          if (viscosity_averaging == MaterialUtilities::harmonic)
            viscosity_averaging_p = -1;
          if (viscosity_averaging == MaterialUtilities::arithmetic)
            viscosity_averaging_p = 1;
          if (viscosity_averaging == MaterialUtilities::maximum_composition)
            viscosity_averaging_p = 1000;


          derivatives->viscosity_derivative_wrt_strain_rate[i] =
            Utilities::derivative_of_weighted_p_norm_average(out.viscosities[i],
                                                             volume_fractions,
                                                             composition_viscosities,
                                                             composition_viscosities_derivatives,
                                                             viscosity_averaging_p);
          derivatives->viscosity_derivative_wrt_pressure[i] =
            Utilities::derivative_of_weighted_p_norm_average(out.viscosities[i],
                                                             volume_fractions,
                                                             composition_viscosities,
                                                             composition_dviscosities_dpressure,
                                                             viscosity_averaging_p);
        }
    }


    template <int dim>
    ComponentMask
    ViscoPlastic<dim>::
    get_volumetric_composition_mask() const
    {
      // Store which components to exclude during the volume fraction computation.
      ComponentMask strain_mask = strain_rheology.get_strain_composition_mask();

      return strain_mask;
    }



    template <int dim>
    void
    ViscoPlastic<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      // Store which components do not represent volumetric compositions (e.g. strain components).
      const ComponentMask volumetric_compositions = get_volumetric_composition_mask();

      EquationOfStateOutputs<dim> eos_outputs (this->n_compositional_fields()+1);

      // Loop through all requested points
      for (unsigned int i=0; i < in.temperature.size(); ++i)
        {
          equation_of_state.evaluate(in, i, eos_outputs);

          // First compute the equation of state variables and thermodynamic properties
          const std::vector<double> volume_fractions = MaterialUtilities::compute_volume_fractions(in.composition[i], volumetric_compositions);

          // not strictly correct if thermal expansivities are different, since we are interpreting
          // these compositions as volume fractions, but the error introduced should not be too bad.
          out.densities[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.densities, MaterialUtilities::arithmetic);
          out.thermal_expansion_coefficients[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.thermal_expansion_coefficients, MaterialUtilities::arithmetic);
          out.specific_heat[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.specific_heat_capacities, MaterialUtilities::arithmetic);
          double thermal_diffusivity = 0.0;

          for (unsigned int j=0; j < volume_fractions.size(); ++j)
            thermal_diffusivity += volume_fractions[j] * thermal_diffusivities[j];

          // Thermal conductivity at the given positions. If the temperature equation uses
          // the reference density profile formulation, use the reference density to
          // calculate thermal conductivity. Otherwise, use the real density. If the adiabatic
          // conditions are not yet initialized, the real density will still be used.
          if (this->get_parameters().formulation_temperature_equation ==
              Parameters<dim>::Formulation::TemperatureEquation::reference_density_profile &&
              this->get_adiabatic_conditions().is_initialized())
            out.thermal_conductivities[i] = thermal_diffusivity * out.specific_heat[i] *
                                            this->get_adiabatic_conditions().density(in.position[i]);
          else
            out.thermal_conductivities[i] = thermal_diffusivity * out.specific_heat[i] * out.densities[i];

          out.compressibilities[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.compressibilities, MaterialUtilities::arithmetic);
          out.entropy_derivative_pressure[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.entropy_derivative_pressure, MaterialUtilities::arithmetic);
          out.entropy_derivative_temperature[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.entropy_derivative_temperature, MaterialUtilities::arithmetic);

          // Compute the effective viscosity if requested and retrieve whether the material is plastically yielding
          bool plastic_yielding = false;
          if (in.strain_rate.size())
            {
              // Currently, the viscosities for each of the compositional fields are calculated assuming
              // isostrain amongst all compositions, allowing calculation of the viscosity ratio.
              // TODO: This is only consistent with viscosity averaging if the arithmetic averaging
              // scheme is chosen. It would be useful to have a function to calculate isostress viscosities.
              const std::pair<std::vector<double>, std::vector<bool> > calculate_viscosities =
                calculate_isostrain_viscosities(volume_fractions, in.pressure[i], in.temperature[i], in.composition[i], in.strain_rate[i],viscous_flow_law,yield_mechanism);

              // The isostrain condition implies that the viscosity averaging should be arithmetic (see above).
              // We have given the user freedom to apply alternative bounds, because in diffusion-dominated
              // creep (where n_diff=1) viscosities are stress and strain-rate independent, so the calculation
              // of compositional field viscosities is consistent with any averaging scheme.
              out.viscosities[i] = MaterialUtilities::average_value(volume_fractions, calculate_viscosities.first, viscosity_averaging);

              // Decide based on the maximum composition if material is yielding.
              // This avoids for example division by zero for harmonic averaging (as plastic_yielding
              // holds values that are either 0 or 1), but might not be consistent with the viscosity
              // averaging chosen.
              std::vector<double>::const_iterator max_composition = std::max_element(volume_fractions.begin(),volume_fractions.end());
              plastic_yielding = calculate_viscosities.second[std::distance(volume_fractions.begin(),max_composition)];

              // Compute viscosity derivatives if they are requested
              if (MaterialModel::MaterialModelDerivatives<dim> *derivatives =
                    out.template get_additional_output<MaterialModel::MaterialModelDerivatives<dim> >())
                compute_viscosity_derivatives(i,volume_fractions, calculate_viscosities.first, in, out);
            }

          // Now compute changes in the compositional fields (i.e. the accumulated strain).
          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            out.reaction_terms[i][c] = 0.0;

          // Calculate changes in strain invariants and update the reaction terms
          strain_rheology.fill_reaction_outputs(in, i, min_strain_rate, plastic_yielding, out);

          // Fill plastic outputs if they exist.
          fill_plastic_outputs(i,volume_fractions,plastic_yielding,in,out);
        }

      // If we use the full strain tensor, compute the change in the individual tensor components.
      strain_rheology.compute_finite_strain_reaction_terms(in, out);
    }

    template <int dim>
    double
    ViscoPlastic<dim>::
    reference_viscosity () const
    {
      return ref_visc;
    }

    template <int dim>
    bool
    ViscoPlastic<dim>::
    is_compressible () const
    {
      return equation_of_state.is_compressible();
    }

    template <int dim>
    double ViscoPlastic<dim>::
    get_min_strain_rate () const
    {
      return min_strain_rate;
    }

    template <int dim>
    void
    ViscoPlastic<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection ("Visco Plastic");
        {
          EquationOfState::MulticomponentIncompressible<dim>::declare_parameters (prm);

          Rheology::StrainDependent<dim>::declare_parameters (prm);

          // Reference and minimum/maximum values
          prm.declare_entry ("Reference temperature", "293", Patterns::Double(0),
                             "For calculating density by thermal expansivity. Units: $\\si{K}$");
          prm.declare_entry ("Minimum strain rate", "1.0e-20", Patterns::Double(0),
                             "Stabilizes strain dependent viscosity. Units: $1 / s$");
          prm.declare_entry ("Reference strain rate","1.0e-15",Patterns::Double(0),
                             "Reference strain rate for first time step. Units: $1 / s$");
          prm.declare_entry ("Minimum viscosity", "1e17", Patterns::Double(0),
                             "Lower cutoff for effective viscosity. Units: $Pa \\, s$");
          prm.declare_entry ("Maximum viscosity", "1e28", Patterns::Double(0),
                             "Upper cutoff for effective viscosity. Units: $Pa \\, s$");
          prm.declare_entry ("Reference viscosity", "1e22", Patterns::Double(0),
                             "Reference viscosity for nondimensionalization. "
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
                             "Units: $Pa \\, s$");

          // Equation of state parameters
          prm.declare_entry ("Thermal diffusivities", "0.8e-6",
                             Patterns::List(Patterns::Double(0)),
                             "List of thermal diffusivities, for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value.  Units: $m^2/s$");

          // Rheological parameters
          prm.declare_entry ("Viscosity averaging scheme", "harmonic",
                             Patterns::Selection("arithmetic|harmonic|geometric|maximum composition"),
                             "When more than one compositional field is present at a point "
                             "with different viscosities, we need to come up with an average "
                             "viscosity at that point.  Select a weighted harmonic, arithmetic, "
                             "geometric, or maximum composition.");
          prm.declare_entry ("Viscous flow law", "composite",
                             Patterns::Selection("diffusion|dislocation|composite"),
                             "Select what type of viscosity law to use between diffusion, "
                             "dislocation and composite options. Soon there will be an option "
                             "to select a specific flow law for each assigned composition ");
          prm.declare_entry ("Yield mechanism", "drucker",
                             Patterns::Selection("drucker|limiter"),
                             "Select what type of yield mechanism to use between Drucker Prager "
                             "and stress limiter options.");

          // Diffusion creep parameters
          Rheology::DiffusionCreep<dim>::declare_parameters(prm);

          // Dislocation creep parameters
          Rheology::DislocationCreep<dim>::declare_parameters(prm);


          // Plasticity parameters
          prm.declare_entry ("Angles of internal friction", "0",
                             Patterns::List(Patterns::Double(0)),
                             "List of angles of internal friction, $\\phi$, for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "For a value of zero, in 2D the von Mises criterion is retrieved. "
                             "Angles higher than 30 degrees are harder to solve numerically. Units: degrees.");
          prm.declare_entry ("Cohesions", "1e20",
                             Patterns::List(Patterns::Double(0)),
                             "List of cohesions, $C$, for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "The extremely large default cohesion value (1e20 Pa) prevents the viscous stress from "
                             "exceeding the yield stress. Units: $Pa$.");

          // Stress limiter parameters
          prm.declare_entry ("Stress limiter exponents", "1.0",
                             Patterns::List(Patterns::Double(0)),
                             "List of stress limiter exponents, $n_{\\text{lim}}$, "
                             "for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "Units: none.");

          // Limit maximum value of the drucker-prager yield stress
          prm.declare_entry ("Maximum yield stress", "1e12", Patterns::Double(0),
                             "Limits the maximum value of the yield stress determined by the "
                             "drucker-prager plasticity parameters. Default value is chosen so this "
                             "is not automatically used. Values of 100e6--1000e6 $Pa$ have been used "
                             "in previous models. Units: $Pa$");

          // Temperature in viscosity laws to include an adiabat (note units of K/Pa)
          prm.declare_entry ("Adiabat temperature gradient for viscosity", "0.0", Patterns::Double(0),
                             "Add an adiabatic temperature gradient to the temperature used in the flow law "
                             "so that the activation volume is consistent with what one would use in a "
                             "earth-like (compressible) model. Default is set so this is off. "
                             "Note that this is a linear approximation of the real adiabatic gradient, which "
                             "is okay for the upper mantle, but is not really accurate for the lower mantle. "
                             "Using a pressure gradient of 32436 Pa/m, then a value of "
                             "0.3 $K/km$ = 0.0003 $K/m$ = 9.24e-09 $K/Pa$ gives an earth-like adiabat."
                             "Units: $K/Pa$");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    ViscoPlastic<dim>::parse_parameters (ParameterHandler &prm)
    {
      // increment by one for background:
      const unsigned int n_fields = this->n_compositional_fields() + 1;

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection ("Visco Plastic");
        {
          // Equation of state parameters
          equation_of_state.initialize_simulator (this->get_simulator());
          equation_of_state.parse_parameters (prm);

          strain_rheology.initialize_simulator (this->get_simulator());
          strain_rheology.parse_parameters(prm);

          // Reference and minimum/maximum values
          reference_T = prm.get_double("Reference temperature");
          min_strain_rate = prm.get_double("Minimum strain rate");
          ref_strain_rate = prm.get_double("Reference strain rate");
          min_visc = prm.get_double ("Minimum viscosity");
          max_visc = prm.get_double ("Maximum viscosity");
          ref_visc = prm.get_double ("Reference viscosity");

          thermal_diffusivities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Thermal diffusivities"))),
                                                                          n_fields,
                                                                          "Thermal diffusivities");


          viscosity_averaging = MaterialUtilities::parse_compositional_averaging_operation ("Viscosity averaging scheme",
                                prm);

          // Rheological parameters
          if (prm.get ("Viscous flow law") == "composite")
            viscous_flow_law = composite;
          else if (prm.get ("Viscous flow law") == "diffusion")
            viscous_flow_law = diffusion;
          else if (prm.get ("Viscous flow law") == "dislocation")
            viscous_flow_law = dislocation;
          else
            AssertThrow(false, ExcMessage("Not a valid viscous flow law"));

          // Rheological parameters
          if (prm.get ("Yield mechanism") == "drucker")
            yield_mechanism = drucker_prager;
          else if (prm.get ("Yield mechanism") == "limiter")
            yield_mechanism = stress_limiter;
          else
            AssertThrow(false, ExcMessage("Not a valid yield mechanism."));

          // Rheological parameters
          // Diffusion creep parameters
          diffusion_creep.initialize_simulator (this->get_simulator());
          diffusion_creep.parse_parameters(prm);

          // Dislocation creep parameters
          dislocation_creep.initialize_simulator (this->get_simulator());
          dislocation_creep.parse_parameters(prm);

          // Plasticity parameters
          angles_internal_friction = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Angles of internal friction"))),
                                                                             n_fields,
                                                                             "Angles of internal friction");
          // Convert angles from degrees to radians
          for (unsigned int i = 0; i<n_fields; ++i)
            angles_internal_friction[i] *= numbers::PI/180.0;
          cohesions = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Cohesions"))),
                                                              n_fields,
                                                              "Cohesions");
          // Stress limiter parameter
          exponents_stress_limiter  = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Stress limiter exponents"))),
                                                                              n_fields,
                                                                              "Stress limiter exponents");

          // Limit maximum value of the drucker-prager yield stress
          max_yield_strength = prm.get_double("Maximum yield stress");

          // Include an adiabat temperature gradient in flow laws
          adiabatic_temperature_gradient_for_viscosity = prm.get_double("Adiabat temperature gradient for viscosity");
          if (this->get_heating_model_manager().adiabatic_heating_enabled())
            AssertThrow (adiabatic_temperature_gradient_for_viscosity == 0.0,
                         ExcMessage("If adiabatic heating is enabled you should not add another adiabatic gradient"
                                    "to the temperature for computing the viscosity, because the ambient"
                                    "temperature profile already includes the adiabatic gradient."));


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

    template <int dim>
    void
    ViscoPlastic<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      if (out.template get_additional_output<PlasticAdditionalOutputs<dim> >() == nullptr)
        {
          const unsigned int n_points = out.viscosities.size();
          out.additional_outputs.push_back(
            std_cxx14::make_unique<MaterialModel::PlasticAdditionalOutputs<dim>> (n_points));
        }
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(ViscoPlastic,
                                   "visco plastic",
                                   "An implementation of a visco-plastic rheology with options for "
                                   "selecting dislocation creep, diffusion creep or composite "
                                   "viscous flow laws.  Plasticity limits viscous stresses through "
                                   "a Drucker Prager yield criterion. The model is incompressible. "
                                   "Note that this material model is based heavily on the "
                                   "DiffusionDislocation (Bob Myhill) and DruckerPrager "
                                   "(Anne Glerum) material models. "
                                   "\n\n "
                                   "The viscosity for dislocation or diffusion creep is defined as "
                                   "\\[v = \\frac 12 A^{-\\frac{1}{n}} d^{\\frac{m}{n}} "
                                   "\\dot{\\varepsilon}_{ii}^{\\frac{1-n}{n}} "
                                   "\\exp\\left(\\frac{E + PV}{nRT}\\right)\\] "
                                   "where $A$ is the prefactor, $n$ is the stress exponent, "
                                   "$\\dot{\\varepsilon}_{ii}$ is the square root of the deviatoric "
                                   "strain rate tensor second invariant, $d$ is grain size, "
                                   "$m$ is the grain size exponent, $E$ is activation energy, "
                                   "$V$ is activation volume, $P$ is pressure, $R$ is the gas "
                                   "exponent and $T$ is temperature. "
                                   "This form of the viscosity equation is commonly used in "
                                   "geodynamic simulations. See, for example, Billen and Hirth "
                                   "(2007), G3, 8, Q08012. Significantly, other studies may use "
                                   "slightly different forms of the viscosity equation leading to "
                                   "variations in how specific terms are defined or combined. For "
                                   "example, the grain size exponent should always be positive in "
                                   "the diffusion viscosity equation used here, while other studies "
                                   "place the grain size term in the denominator and invert the sign "
                                   "of the grain size exponent. When examining previous work, one "
                                   "should carefully check how the viscous prefactor and grain size "
                                   "terms are defined. "
                                   "\n\n "
                                   "One may select to use the diffusion ($v_{\\text{diff}}$; $n=1$, $m!=0$), "
                                   "dislocation ($v_{\\text{disl}}$, $n>1$, $m=0$) or composite "
                                   "$\\frac{v_{\\text{diff}} v_{\\text{disl}}}{v_{\\text{diff}}+v_{\\text{disl}}}$ equation form. "
                                   "\n\n "
                                   "The diffusion and dislocation prefactors can be weakened with a factor "
                                   "between 0 and 1 according to the total or the viscous strain only. "
                                   "\n\n "
                                   "Viscosity is limited through one of two different `yielding' mechanisms. "
                                   "\n\n"
                                   "Plasticity limits viscous stress through a Drucker Prager "
                                   "yield criterion, where the yield stress in 3D is  "
                                   "$\\sigma_y = \\frac{6C\\cos(\\phi) + 2P\\sin(\\phi)} "
                                   "{\\sqrt(3)(3+\\sin(\\phi))}$ "
                                   "and "
                                   "$\\sigma_y = C\\cos(\\phi) + P\\sin(\\phi)$ "
                                   "in 2D. Above, $C$ is cohesion and $\\phi$  is the angle of "
                                   "internal friction.  Note that the 2D form is equivalent to the "
                                   "Mohr Coulomb yield surface.  If $\\phi$ is 0, the yield stress "
                                   "is fixed and equal to the cohesion (Von Mises yield criterion). "
                                   "When the viscous stress ($2v{\\varepsilon}_{ii}$) exceeds "
                                   "the yield stress, the viscosity is rescaled back to the yield "
                                   "surface: $v_{y}=\\sigma_{y}/(2{\\varepsilon}_{ii})$. "
                                   "This form of plasticity is commonly used in geodynamic models. "
                                   "See, for example, Thieulot, C. (2011), PEPI 188, pp. 47-68. "
                                   "\n\n"
                                   "The user has the option to linearly reduce the cohesion and "
                                   "internal friction angle as a function of the finite strain magnitude. "
                                   "The finite strain invariant or full strain tensor is calculated through "
                                   "compositional fields within the material model. This implementation is "
                                   "identical to the compositional field finite strain plugin and cookbook "
                                   "described in the manual (author: Gassmoeller, Dannberg). If the user selects to track "
                                   "the finite strain invariant ($e_{ii}$), a single compositional field tracks "
                                   "the value derived from $e_{ii}^t = (e_{ii})^{(t-1)} + \\dot{e}_{ii}\\; dt$, where $t$ and $t-1$ "
                                   "are the current and prior time steps, $\\dot{e}_{ii}$ is the second invariant of the "
                                   "strain rate tensor and $dt$ is the time step size. In the case of the "
                                   "full strain tensor $F$, the finite strain magnitude is derived from the "
                                   "second invariant of the symmetric stretching tensor $L$, where "
                                   "$L = F [F]^T$. The user must specify a single compositional "
                                   "field for the finite strain invariant or multiple fields (4 in 2D, 9 in 3D) "
                                   "for the finite strain tensor. These field(s) must be the first listed "
                                   "compositional fields in the parameter file. Note that one or more of the finite strain "
                                   "tensor components must be assigned a non-zero value initially. This value can be "
                                   "be quite small (e.g., 1.e-8), but still non-zero. While the option to track and use "
                                   "the full finite strain tensor exists, tracking the associated compositional fields "
                                   "is computationally expensive in 3D. Similarly, the finite strain magnitudes "
                                   "may in fact decrease if the orientation of the deformation field switches "
                                   "through time. Consequently, the ideal solution is track the finite strain "
                                   "invariant (single compositional) field within the material and track "
                                   "the full finite strain tensor through particles."
                                   "When only the second invariant of the strain is tracked, one has the option to "
                                   "track the full strain or only the plastic strain. In the latter case, strain is only tracked "
                                   "in case the material is plastically yielding, i.e. the viscous stess > yield strength. "
                                   ""
                                   "\n\n"
                                   "Viscous stress may also be limited by a non-linear stress limiter "
                                   "that has a form similar to the Peierls creep mechanism. "
                                   "This stress limiter assigns an effective viscosity "
                                   "$\\sigma_{\\text{eff}} = \\frac{\\tau_y}{2\\varepsilon_y} "
                                   "{\\frac{\\varepsilon_{ii}}{\\varepsilon_y}}^{\\frac{1}{n_y}-1}$ "
                                   "Above $\\tau_y$ is a yield stress, $\\varepsilon_y$ is the "
                                   "reference strain rate, $\\varepsilon_{ii}$ is the strain rate "
                                   "and $n_y$ is the stress limiter exponent.  The yield stress, "
                                   "$\\tau_y$, is defined through the Drucker Prager yield criterion "
                                   "formulation. This method of limiting viscous stress has been used "
                                   "in various forms within the geodynamic literature \\cite{chri92,vavv02,cibi13,cibi15}."
                                   "When $n_y$ is 1, it essentially becomes a linear viscosity model, "
                                   "and in the limit $n_y\\rightarrow \\infty$ it converges to the "
                                   "standard viscosity rescaling method (concretely, values $n_y>20$ "
                                   "are large enough)."
                                   "\n\n "
                                   "Compositional fields can each be assigned individual values of "
                                   "thermal diffusivity, heat capacity, density, thermal "
                                   "expansivity and rheological parameters. "
                                   "\n\n "
                                   "If more than one compositional field is present at a given "
                                   "point, viscosities are averaged with an arithmetic, geometric "
                                   "harmonic (default) or maximum composition scheme. "
                                   "\n\n "
                                   "The value for the components of this formula and additional "
                                   "parameters are read from the parameter file in subsection "
                                   " 'Material model/Visco Plastic'.")
  }
}
