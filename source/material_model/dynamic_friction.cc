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


#include <aspect/material_model/dynamic_friction.h>
#include <aspect/simulator.h>
#include <aspect/utilities.h>

#include <numeric>


namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    const std::vector<double>
    DynamicFriction<dim>::
    compute_viscosities(const double pressure,
                        const SymmetricTensor<2,dim> &strain_rate) const
    {

      std::vector<double> viscosities(mu_s.size());

      // second invariant for strain tensor
      const double strain_rate_dev_inv2 = ( (this->get_timestep_number() == 0 && strain_rate.norm() <= std::numeric_limits<double>::min())
                                            ?
                                            reference_strain_rate * reference_strain_rate
                                            :
                                            std::fabs(second_invariant(deviator(strain_rate))));

      for (unsigned int i = 0; i < mu_s.size(); i++)
        {
          // Calculate effective steady-state friction coefficient. The formula below is equivalent to the
          // equation 13 in van Dinther et al., (2013, JGR). Although here the dynamic friction coefficient
          // is directly specified. In addition, we also use a reference strain rate in place of a characteristic
          // velocity divided by local element size.
          const double mu  = mu_d[i] + (mu_s[i] - mu_d[i]) / ( (1 + strain_rate_dev_inv2/reference_strain_rate) );

          // Convert effective steady-state friction coefficient to internal angle of friction.
          const double phi = std::atan (mu);

          // Compute the viscosity according to the Drucker-Prager yield criterion.
          const double plastic_viscosity = drucker_prager_plasticity.compute_viscosity(cohesions[i],
                                                                                       phi,
                                                                                       std::max(pressure,0.0),
                                                                                       std::sqrt(strain_rate_dev_inv2),
                                                                                       std::numeric_limits<double>::infinity());

          // Cut off the viscosity between a minimum and maximum value to avoid
          // a numerically unfavourable large viscosity range.
          viscosities[i] = 1.0 / ( ( 1.0 / (plastic_viscosity + minimum_viscosity) ) + (1.0 / maximum_viscosity) );

        }
      return viscosities;
    }



    template <int dim>
    void
    DynamicFriction<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      EquationOfStateOutputs<dim> eos_outputs (this->n_compositional_fields()+1);

      for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
        {
          const std::vector<double> volume_fractions = MaterialUtilities::compute_field_fractions(in.composition[i]);

          if (in.requests_property(MaterialProperties::viscosity))
            {
              const std::vector<double> viscosities = compute_viscosities(in.pressure[i], in.strain_rate[i]);
              out.viscosities[i] = MaterialUtilities::average_value (volume_fractions, viscosities, viscosity_averaging);
            }

          equation_of_state.evaluate(in, i, eos_outputs);

          // The averaging is not strictly correct if thermal expansivities are different, since we are interpreting
          // these compositions as volume fractions, but the error introduced should not be too bad.
          out.densities[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.densities, MaterialUtilities::arithmetic);
          out.thermal_expansion_coefficients[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.thermal_expansion_coefficients, MaterialUtilities::arithmetic);
          out.specific_heat[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.specific_heat_capacities, MaterialUtilities::arithmetic);

          // Arithmetic averaging of thermal conductivities
          // This may not be strictly the most reasonable thing, but for most Earth materials we hope
          // that they do not vary so much that it is a big problem.
          out.thermal_conductivities[i] = MaterialUtilities::average_value (volume_fractions, thermal_conductivities, MaterialUtilities::arithmetic);

          out.compressibilities[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.compressibilities, MaterialUtilities::arithmetic);
          out.entropy_derivative_pressure[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.entropy_derivative_pressure, MaterialUtilities::arithmetic);
          out.entropy_derivative_temperature[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.entropy_derivative_temperature, MaterialUtilities::arithmetic);

          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            out.reaction_terms[i][c] = 0.0;

        }
    }

    template <int dim>
    double
    DynamicFriction<dim>::
    reference_viscosity () const
    {
      return background_viscosities[0]; //background
    }

    template <int dim>
    bool
    DynamicFriction<dim>::
    is_compressible () const
    {
      return equation_of_state.is_compressible();
    }

    template <int dim>
    void
    DynamicFriction<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Dynamic Friction");
        {
          EquationOfState::MulticomponentIncompressible<dim>::declare_parameters (prm, 4.e-5);

          prm.declare_entry ("Reference temperature", "293.",
                             Patterns::Double (0.),
                             "The reference temperature $T_0$. Units: \\si{\\kelvin}.");
          prm.declare_entry ("Thermal conductivities", "4.7",
                             Patterns::List(Patterns::Double (0.)),
                             "List of thermal conductivities for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one value is given, then all use the same value. Units: \\si{\\watt\\per\\meter\\per\\kelvin}.");
          prm.declare_entry ("Viscosity averaging scheme", "harmonic",
                             Patterns::Selection("arithmetic|harmonic|geometric|maximum composition"),
                             "When more than one compositional field is present at a point "
                             "with different viscosities, we need to come up with an average "
                             "viscosity at that point.  Select a weighted harmonic, arithmetic, "
                             "geometric, or maximum composition.");
          prm.enter_subsection("Viscosities");
          {
            prm.declare_entry ("Minimum viscosity", "1e19",
                               Patterns::Double (0.),
                               "The value of the minimum viscosity cutoff $\\eta_min$. "
                               "Units: \\si{\\pascal\\second}.");
            prm.declare_entry ("Maximum viscosity", "1e24",
                               Patterns::Double (0.),
                               "The value of the maximum viscosity cutoff $\\eta_max$. "
                               "Units: \\si{\\pascal\\second}.");
            prm.declare_entry ("Reference strain rate", "1e-15",
                               Patterns::Double (0.),
                               "The value of the initial strain rate prescribed during the "
                               "first nonlinear iteration $\\dot{\\epsilon}_ref$. Units: \\si{\\per\\second}.");
            prm.declare_entry ("Coefficients of static friction", "0.5",
                               Patterns::List(Patterns::Double (0.)),
                               "List of coefficients of static friction for background mantle and compositional fields,"
                               "for a total of N+1 values, where N is the number of compositional fields."
                               "If only one value is given, then all use the same value. Units: dimensionless.");
            prm.declare_entry ("Coefficients of dynamic friction", "0.4",
                               Patterns::List(Patterns::Double (0.)),
                               "List of coefficients of dynamic friction for background mantle and compositional fields,"
                               "for a total of N+1 values, where N is the number of compositional fields."
                               "If only one value is given, then all use the same value. Units: dimensionless.");
            prm.declare_entry ("Cohesions", "4.e6",
                               Patterns::List(Patterns::Double (0.)),
                               "List of cohesions for background mantle and compositional fields,"
                               "for a total of N+1 values, where N is the number of compositional fields."
                               "If only one value is given, then all use the same value. Units: \\si{\\pascal}.");
            prm.declare_entry ("Background Viscosities", "1.e20",
                               Patterns::List(Patterns::Double (0.)),
                               "List of background viscosities for mantle and compositional fields,"
                               "for a total of N+1 values, where N is the number of compositional fields."
                               "If only one value is given, then all use the same value. "
                               "Units: \\si{\\pascal\\second}.");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    DynamicFriction<dim>::parse_parameters (ParameterHandler &prm)
    {
      //not pretty, but we need to get the number of compositional fields before
      //simulatoraccess has been initialized here...

      prm.enter_subsection ("Compositional fields");
      const unsigned int n_fields = this->n_compositional_fields() + 1;
      prm.leave_subsection();

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Dynamic Friction");
        {
          equation_of_state.initialize_simulator (this->get_simulator());
          equation_of_state.parse_parameters (prm);

          reference_T = prm.get_double ("Reference temperature");

          viscosity_averaging = MaterialUtilities::parse_compositional_averaging_operation ("Viscosity averaging scheme",
                                prm);

          thermal_conductivities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Thermal conductivities"))),
                                                                           n_fields,
                                                                           "Thermal conductivities");
          prm.enter_subsection("Viscosities");
          {
            minimum_viscosity  = prm.get_double ("Minimum viscosity");
            maximum_viscosity  = prm.get_double ("Maximum viscosity");
            reference_strain_rate    = prm.get_double ("Reference strain rate");

            mu_s = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Coefficients of static friction"))),
                                                           n_fields,
                                                           "Coefficients of static friction");
            mu_d = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Coefficients of dynamic friction"))),
                                                           n_fields,
                                                           "Coefficients of dynamic friction");
            cohesions = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Cohesions"))),
                                                                n_fields,
                                                                "Cohesions");
            background_viscosities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Background Viscosities"))),
                                                                             n_fields,
                                                                             "Background Viscosities");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Declare dependencies on solution variables
      this->model_dependence.viscosity = NonlinearDependence::compositional_fields | NonlinearDependence::strain_rate;
      this->model_dependence.density = NonlinearDependence::temperature | NonlinearDependence::compositional_fields;
      this->model_dependence.compressibility = NonlinearDependence::none;
      this->model_dependence.specific_heat = NonlinearDependence::compositional_fields;
      this->model_dependence.thermal_conductivity = NonlinearDependence::compositional_fields;
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(DynamicFriction,
                                   "dynamic friction",
                                   "This model is for use with an arbitrary number of compositional fields, where each field "
                                   "represents a rock type which can have completely different properties from the others."
                                   "Each rock type itself has constant material properties, with the exception of viscosity "
                                   "which is modified according to a Drucker-Prager yield criterion. Unlike the drucker prager "
                                   "or visco plastic material models, the angle of internal friction is a function of velocity. "
                                   "This relationship is similar to rate-and-state friction constitutive relationships, which "
                                   "are applicable to the strength of rocks during earthquakes. The formulation used here is "
                                   "derived from van Dinther et al. 2013, JGR. Each compositional field is interpreed as a volume fraction. "
                                   "If the sum of the fields is greater than one, they are renormalized. If it is less than one, material properties "
                                   "for ``background material'' make up the rest. When more than one field is present, the "
                                   "material properties are averaged arithmetically. An exception is the viscosity, "
                                   "where the averaging should make more of a difference. For this, the user selects"
                                   "between arithmetic, harmonic, geometric, or maximum composition averaging. ")
  }
}
