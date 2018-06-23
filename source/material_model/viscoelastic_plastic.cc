/*
  Copyright (C) 2011 - 2018 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/material_model/viscoelastic_plastic.h>
#include <aspect/utilities.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <aspect/global.h>
#include <numeric>
#include <deal.II/base/signaling_nan.h>

namespace aspect
{
  namespace MaterialModel
  {

    namespace
    {
      std::vector<std::string> make_elastic_additional_outputs_names()
      {
        std::vector<std::string> names;
        names.push_back("elastic_shear_modulus");
        return names;
      }
    }

    template <int dim>
    ElasticAdditionalOutputs<dim>::ElasticAdditionalOutputs (const unsigned int n_points)
      :
      NamedAdditionalMaterialOutputs<dim>(make_elastic_additional_outputs_names()),
      elastic_shear_moduli(n_points, numbers::signaling_nan<double>())
    {}

    template <int dim>
    std::vector<double>
    ElasticAdditionalOutputs<dim>::get_nth_output(const unsigned int idx) const
    {
      AssertIndexRange (idx, 1);
      switch (idx)
        {
          case 0:
            return elastic_shear_moduli;

          default:
            AssertThrow(false, ExcInternalError());
        }
      // we will never get here, so just return something
      return elastic_shear_moduli;
    }

    template <int dim>
    double
    ViscoelasticPlastic<dim>::
    average_value (const std::vector<double> &volume_fractions,
                   const std::vector<double> &parameter_values,
                   const enum AveragingScheme &average_type) const
    {
      double averaged_parameter = 0.0;

      switch (average_type)
        {
          case arithmetic:
          {
            for (unsigned int i=0; i< volume_fractions.size(); ++i)
              averaged_parameter += volume_fractions[i]*parameter_values[i];
            break;
          }
          case harmonic:
          {
            for (unsigned int i=0; i< volume_fractions.size(); ++i)
              {
                AssertThrow(parameter_values[i] != 0,
                            ExcMessage ("All values must be greater than 0 during harmonic averaging"));
                averaged_parameter += volume_fractions[i]/(parameter_values[i]);
              }
            averaged_parameter = 1.0/averaged_parameter;
            break;
          }
          case geometric:
          {
            for (unsigned int i=0; i < volume_fractions.size(); ++i)
              averaged_parameter += volume_fractions[i]*std::log(parameter_values[i]);
            averaged_parameter = std::exp(averaged_parameter);
            break;
          }
          case maximum_composition:
          {
            const unsigned int i = (unsigned int)(std::max_element( volume_fractions.begin(),
                                                                    volume_fractions.end() )
                                                  - volume_fractions.begin());
            averaged_parameter = parameter_values[i];
            break;
          }
          default:
          {
            AssertThrow( false, ExcNotImplemented() );
            break;
          }
        }
      return averaged_parameter;
    }

    template <int dim>
    void
    ViscoelasticPlastic<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {

      // Store which components to exclude during volume fraction computation.
      ComponentMask composition_mask(this->n_compositional_fields(),true);
      // Elastic stress fields
      for (unsigned int i = 0; i < SymmetricTensor<2,dim>::n_independent_components ; ++i)
        composition_mask.set(i,false);

      // Create the structure for the elastic force terms that are needed to compute the
      // right-hand side of the Stokes system
      MaterialModel::ElasticOutputs<dim>
      *force_out = out.template get_additional_output<MaterialModel::ElasticOutputs<dim> >();

      // The elastic time step (dte) is equal to the numerical time step if the time step number
      // is greater than 0 and the parameter 'use_fixed_elastic_time_step' is set to false.
      // On the first (0) time step the elastic time step is always equal to the value
      // specified in 'fixed_elastic_time_step', which is also used in all subsequent time
      // steps if 'use_fixed_elastic_time_step' is set to true.
      const double dte = ( ( this->get_timestep_number() > 0 && use_fixed_elastic_time_step == false )
                           ?
                           this->get_timestep()
                           :
                           fixed_elastic_time_step);

      for (unsigned int i=0; i < in.temperature.size(); ++i)
        {

          const double temperature = in.temperature[i];
          const std::vector<double> composition = in.composition[i];
          const double pressure = in.pressure[i];
          const SymmetricTensor<2,dim> strain_rate = in.strain_rate[i];
          const std::vector<double> volume_fractions = compute_volume_fractions(composition, composition_mask);

          // Calculate the square root of the second moment invariant for the deviatoric strain rate tensor.
          // The first time this function is called (first iteration of first time step)
          // a specified "reference" strain rate is used as the returned value would
          // otherwise be zero.
          const double edot_ii = ( (this->get_timestep_number() == 0 && strain_rate.norm() <= std::numeric_limits<double>::min())
                                   ?
                                   reference_strain_rate
                                   :
                                   std::max(std::sqrt(std::fabs(second_invariant(deviator(strain_rate)))),
                                            minimum_strain_rate) );

          out.specific_heat[i] = average_value(volume_fractions, specific_heats, arithmetic);

          // Arithmetic averaging of thermal conductivities
          // This may not be strictly the most reasonable thing, but for most Earth materials we hope
          // that they do not vary so much that it is a big problem.
          out.thermal_conductivities[i] = average_value(volume_fractions, thermal_conductivities, arithmetic);

          double density = 0.0;
          for (unsigned int j=0; j < volume_fractions.size(); ++j)
            {
              // not strictly correct if thermal expansivities are different, since we are interpreting
              // these compositions as volume fractions, but the error introduced should not be too bad.
              const double temperature_factor= (1.0 - thermal_expansivities[j] * (temperature - reference_temperature));
              density += volume_fractions[j] * densities[j] * temperature_factor;
            }
          out.densities[i] = density;

          out.thermal_expansion_coefficients[i] = average_value(volume_fractions, thermal_expansivities, arithmetic);

          // Compressibility at the given positions.
          // The compressibility is given as
          // $\frac 1\rho \frac{\partial\rho}{\partial p}$.
          // (here we use an incompressible medium)
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

          const std::vector<double> viscosities_pre_yield = linear_viscosities;

          // TODO: Add strain-weakening of cohesion and friction
          const std::vector<double> coh = cohesions;
          const std::vector<double> phi = angles_internal_friction;

          // Initialize variables
          std::vector<double> stresses_viscous(volume_fractions.size());
          std::vector<double> stresses_yield(volume_fractions.size());
          std::vector<double> viscosities_drucker_prager(volume_fractions.size());
          std::vector<double> viscosities_yield(volume_fractions.size());
          std::vector<double> viscosities_viscoelastic(volume_fractions.size());

          // Loop through all compositions
          for (unsigned int j=0; j < volume_fractions.size(); ++j)
            {

              stresses_viscous[j] = 2. * viscosities_pre_yield[j] * edot_ii;

              stresses_yield[j] = ( (dim==3)
                                    ?
                                    ( 6.0 * coh[j] * std::cos(phi[j]) + 6.0 * std::max(pressure,0.0) * std::sin(phi[j]) )
                                    / ( std::sqrt(3.0) * (3.0 + std::sin(phi[j]) ) )
                                    :
                                    coh[j] * std::cos(phi[j]) + std::max(pressure,0.0) * std::sin(phi[j]) );

              if ( stresses_viscous[j] >= stresses_yield[j]  )
                {
                  viscosities_drucker_prager[j] = stresses_yield[j] / (2.0 * edot_ii);
                }
              else
                {
                  viscosities_drucker_prager[j] = viscosities_pre_yield[j];
                }

              viscosities_yield[j] = viscosities_drucker_prager[j];

              viscosities_viscoelastic[j] = ( viscosities_yield[j] * dte ) / ( dte + ( viscosities_yield[j] / elastic_shear_moduli[j] ) );

            }

          out.viscosities[i] = average_value(volume_fractions,viscosities_viscoelastic,viscosity_averaging);

          if (ElasticAdditionalOutputs<dim> *elastic_out = out.template get_additional_output<ElasticAdditionalOutputs<dim> >())
            {
              elastic_out->elastic_shear_moduli[i] = average_value(volume_fractions,elastic_shear_moduli,viscosity_averaging);
            }

          // Fill elastic force outputs (assumed to be zero during intial time step)
          if (force_out)
            {
              force_out->elastic_force[i] = 0.;
            }
        }

      // Viscoelasticity section
      if (in.current_cell.state() == IteratorState::valid && this->get_timestep_number() > 0 && in.strain_rate.size() > 0)
        {
          // Get old (previous time step) velocity gradients
          std::vector<Point<dim> > quadrature_positions(in.position.size());
          for (unsigned int i=0; i < in.position.size(); ++i)
            quadrature_positions[i] = this->get_mapping().transform_real_to_unit_cell(in.current_cell, in.position[i]);

          // FEValues requires a quadrature and we provide the default quadrature
          // as we only need to evaluate the solution and gradients.
          FEValues<dim> fe_values (this->get_mapping(),
                                   this->get_fe(),
                                   Quadrature<dim>(quadrature_positions),
                                   update_gradients);

          fe_values.reinit (in.current_cell);
          std::vector<Tensor<2,dim> > old_velocity_gradients (quadrature_positions.size(), Tensor<2,dim>());
          fe_values[this->introspection().extractors.velocities].get_function_gradients (this->get_old_solution(),
                                                                                         old_velocity_gradients);

          MaterialModel::ElasticOutputs<dim>
          *force_out = out.template get_additional_output<MaterialModel::ElasticOutputs<dim> >();

          for (unsigned int i=0; i < in.position.size(); ++i)
            {

              const std::vector<double> composition = in.composition[i];
              const std::vector<double> volume_fractions = compute_volume_fractions(composition, composition_mask);

              // Get old stresses from compositional fields
              SymmetricTensor<2,dim> stress_old;
              for (unsigned int j=0; j < SymmetricTensor<2,dim>::n_independent_components; ++j)
                stress_old[SymmetricTensor<2,dim>::unrolled_to_component_indices(j)] = in.composition[i][j];

              // Calculate the rotated stresses
              // Rotation (vorticity) tensor (equation 25 in Moresi et al., 2003, J. Comp. Phys.)
              const Tensor<2,dim> rotation = 0.5 * ( old_velocity_gradients[i] - transpose(old_velocity_gradients[i]) );

              // Recalculate average elastic shear modulus
              const double elastic_shear_modulus = average_value(volume_fractions,
                                                                 elastic_shear_moduli,
                                                                 viscosity_averaging);

              // Average viscoelastic viscosity
              const double viscoelastic_viscosity = out.viscosities[i];

              // Calculate the current (new) viscoelastic stress, which is a function of the material
              // properties (viscoelastic viscosity, shear modulus), elastic time step size, strain rate,
              // vorticity and prior (inherited) viscoelastic stresses (see equation 29 in Moresi et al.,
              // 2003, J. Comp. Phys.)
              SymmetricTensor<2,dim> stress_new = ( 2. * viscoelastic_viscosity * deviator(in.strain_rate[i]) ) +
                                                  ( ( viscoelastic_viscosity / ( elastic_shear_modulus * dte ) ) * stress_old ) +
                                                  ( ( viscoelastic_viscosity / elastic_shear_modulus ) *
                                                    ( symmetrize(rotation * Tensor<2,dim>(stress_old) ) - symmetrize(Tensor<2,dim>(stress_old) * rotation) ) );

              // Stress averaging scheme to account for difference betweed fixed elastic time step
              // and numerical time step (see equation 32 in Moresi et al., 2003, J. Comp. Phys.)
              const double dt = this->get_timestep();
              if (use_fixed_elastic_time_step == true && use_stress_averaging == true)
                {
                  stress_new = ( ( 1. - ( dt / dte ) ) * stress_old ) + ( ( dt / dte ) * stress_new ) ;
                }

              // Fill reaction terms
              for (unsigned int j = 0; j < SymmetricTensor<2,dim>::n_independent_components ; ++j)
                out.reaction_terms[i][j] = -stress_old[SymmetricTensor<2,dim>::unrolled_to_component_indices(j)]
                                           + stress_new[SymmetricTensor<2,dim>::unrolled_to_component_indices(j)];

              // Fill elastic force outputs (See equation 30 in Moresi et al., 2003, J. Comp. Phys.)
              if (force_out)
                {
                  force_out->elastic_force[i] = -1. * ( ( viscoelastic_viscosity / ( elastic_shear_modulus * dte  ) ) * stress_old );
                }

            }
        }

    }

    template <int dim>
    double
    ViscoelasticPlastic<dim>::
    reference_viscosity () const
    {
      return input_reference_viscosity;
    }

    template <int dim>
    bool
    ViscoelasticPlastic<dim>::
    is_compressible () const
    {
      return false;
    }

    template <int dim>
    void
    ViscoelasticPlastic<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Viscoelastic Plastic");
        {

          // Reference and minimum/maximum values
          prm.declare_entry ("Reference temperature", "293",
                             Patterns::Double (0),
                             "The reference temperature $T_0$. Units: $K$.");
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
          prm.declare_entry ("Densities", "3300.",
                             Patterns::List(Patterns::Double(0)),
                             "List of densities for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value.  Units: $kg / m^3$");
          prm.declare_entry ("Thermal expansivities", "4.e-5",
                             Patterns::List(Patterns::Double(0)),
                             "List of thermal expansivities for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value. Units: $1/K$");
          prm.declare_entry ("Specific heats", "1250.",
                             Patterns::List(Patterns::Double(0)),
                             "List of specific heats $C_p$ for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value. Units: $J /kg /K$");
          prm.declare_entry ("Thermal conductivities", "4.7",
                             Patterns::List(Patterns::Double(0)),
                             "List of thermal conductivities for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value. Units: $W/m/K$ ");
          prm.declare_entry ("Viscosity averaging scheme", "harmonic",
                             Patterns::Selection("arithmetic|harmonic|geometric|maximum composition "),
                             "When more than one compositional field is present at a point "
                             "with different viscosities, we need to come up with an average "
                             "viscosity at that point.  Select a weighted harmonic, arithmetic, "
                             "geometric, or maximum composition.");

          prm.declare_entry ("Linear viscosities", "1.e21",
                             Patterns::List(Patterns::Double(0)),
                             "List of linear (fixed) viscosities for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value. Units: $Pa s$");

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

          // Viscoelasticity parameters
          prm.declare_entry ("Elastic shear moduli", "75.0e9",
                             Patterns::List(Patterns::Double(0)),
                             "List of elastic shear moduli, $G$, "
                             "for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "The default value of 75 GPa is representative of mantle rocks. Units: Pa.");
          prm.declare_entry ("Use fixed elastic time step", "unspecified",
                             Patterns::Selection("true|false|unspecified"),
                             "Select whether the material time scale in the viscoelastic constitutive "
                             "relationship uses the regular numerical time step or a separate fixed "
                             "elastic time step throughout the model run. The fixed elastic time step "
                             "is always used during the initial time step. If a fixed elastic time "
                             "step is used throughout the model run, a stress averaging scheme can be "
                             "applied to account for differences with the numerical time step. An "
                             "alternative approach is to limit the maximum time step size so that it "
                             "is equal to the elastic time step. The default value of this parameter is "
                             "'unspecified', which throws an exception during runtime. In order for "
                             "the model to run the user must select 'true' or 'false'.");
          prm.declare_entry ("Fixed elastic time step", "1.e3",
                             Patterns::Double (0),
                             "The fixed elastic time step $dte$. Units: years if the "
                             "'Use years in output instead of seconds' parameter is set; "
                             "seconds otherwise.");
          prm.declare_entry ("Use stress averaging","false",
                             Patterns::Bool (),
                             "Whether to apply a stress averaging scheme to account for differences "
                             "between the fixed elastic time step and numerical time step. ");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    ViscoelasticPlastic<dim>::parse_parameters (ParameterHandler &prm)
    {

      // Get the number of fields for composition-dependent material properties
      const unsigned int n_fields = this->n_compositional_fields() + 1;

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Viscoelastic Plastic");
        {

          reference_temperature = prm.get_double ("Reference temperature");
          minimum_strain_rate = prm.get_double("Minimum strain rate");
          reference_strain_rate = prm.get_double("Reference strain rate");
          minimum_viscosity = prm.get_double ("Minimum viscosity");
          maximum_viscosity = prm.get_double ("Maximum viscosity");
          input_reference_viscosity = prm.get_double ("Reference viscosity");

          if (prm.get ("Viscosity averaging scheme") == "harmonic")
            viscosity_averaging = harmonic;
          else if (prm.get ("Viscosity averaging scheme") == "arithmetic")
            viscosity_averaging = arithmetic;
          else if (prm.get ("Viscosity averaging scheme") == "geometric")
            viscosity_averaging = geometric;
          else if (prm.get ("Viscosity averaging scheme") == "maximum composition")
            viscosity_averaging = maximum_composition;
          else
            AssertThrow(false, ExcMessage("Not a valid viscosity averaging scheme"));


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

          // Parse viscoelastic properties
          densities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Densities"))),
                                                              n_fields,
                                                              "Densities");
          linear_viscosities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Linear viscosities"))),
                                                                       n_fields,
                                                                       "Viscosities");
          thermal_conductivities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Thermal conductivities"))),
                                                                           n_fields,
                                                                           "Thermal conductivities");
          thermal_expansivities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Thermal expansivities"))),
                                                                          n_fields,
                                                                          "Thermal expansivities");
          specific_heats = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Specific heats"))),
                                                                   n_fields,
                                                                   "Specific heats");
          elastic_shear_moduli = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Elastic shear moduli"))),
                                                                         n_fields,
                                                                         "Elastic shear moduli");

          if (prm.get ("Use fixed elastic time step") == "true")
            use_fixed_elastic_time_step = true;
          else if (prm.get ("Use fixed elastic time step") == "false")
            use_fixed_elastic_time_step = false;
          else
            AssertThrow(false, ExcMessage("'Use fixed elastic time step' must be set to 'true' or 'false'"));

          use_stress_averaging = prm.get_bool ("Use stress averaging");
          if (use_stress_averaging)
            AssertThrow(use_fixed_elastic_time_step == true,
                        ExcMessage("Stress averaging can only be used if 'Use fixed elastic time step' is set to true'"));

          fixed_elastic_time_step = prm.get_double ("Fixed elastic time step");
          AssertThrow(fixed_elastic_time_step > 0,
                      ExcMessage("The fixed elastic time step must be greater than zero"));

          if (this->convert_output_to_years())
            fixed_elastic_time_step *= year_in_seconds;

          AssertThrow(this->get_parameters().enable_elasticity == true,
                      ExcMessage ("Material model Viscoelastic only works if 'Enable elasticity' is set to true"));

          // Check whether the compositional fields representing the viscoelastic
          // stress tensor are both named correctly and listed in the right order.
          if (dim == 2)
            {
              AssertThrow(this->introspection().compositional_index_for_name("stress_xx") == 0,
                          ExcMessage("Material model Viscoelastic only works if the first "
                                     "compositional field is called stress_xx."));
              AssertThrow(this->introspection().compositional_index_for_name("stress_yy") == 1,
                          ExcMessage("Material model Viscoelastic only works if the second "
                                     "compositional field is called stress_yy."));
              AssertThrow(this->introspection().compositional_index_for_name("stress_xy") == 2,
                          ExcMessage("Material model Viscoelastic only works if the third "
                                     "compositional field is called stress_xy."));
            }
          else if (dim == 3)
            {
              AssertThrow(this->introspection().compositional_index_for_name("stress_xx") == 0,
                          ExcMessage("Material model Viscoelastic Plastic only works if the first "
                                     "compositional field is called stress_xx."));
              AssertThrow(this->introspection().compositional_index_for_name("stress_yy") == 1,
                          ExcMessage("Material model Viscoelastic Plastic only works if the second "
                                     "compositional field is called stress_yy."));
              AssertThrow(this->introspection().compositional_index_for_name("stress_zz") == 2,
                          ExcMessage("Material model Viscoelastic Plastic only works if the third "
                                     "compositional field is called stress_zz."));
              AssertThrow(this->introspection().compositional_index_for_name("stress_xy") == 3,
                          ExcMessage("Material model Viscoelastic Plastic only works if the fourth "
                                     "compositional field is called stress_xy."));
              AssertThrow(this->introspection().compositional_index_for_name("stress_xz") == 4,
                          ExcMessage("Material model Viscoelastic Plastic only works if the fifth "
                                     "compositional field is called stress_xz."));
              AssertThrow(this->introspection().compositional_index_for_name("stress_yz") == 5,
                          ExcMessage("Material model Viscoelastic Plastic only works if the sixth "
                                     "compositional field is called stress_yz."));
            }
          else
            ExcNotImplemented();

          // It does not make sense to use this material with the 'no Advection, iterated Stokes'
          // or 'single Advection, no Stokes' solver schemes.
          AssertThrow((this->get_parameters().nonlinear_solver ==
                       Parameters<dim>::NonlinearSolver::single_Advection_single_Stokes
                       ||
                       this->get_parameters().nonlinear_solver ==
                       Parameters<dim>::NonlinearSolver::single_Advection_iterated_Stokes
                       ||
                       this->get_parameters().nonlinear_solver ==
                       Parameters<dim>::NonlinearSolver::iterated_Advection_and_Stokes
                       ||
                       this->get_parameters().nonlinear_solver ==
                       Parameters<dim>::NonlinearSolver::iterated_Advection_and_Newton_Stokes),
                      ExcMessage("The material model will now work with the nonlinear "
                                 "solver schemes 'no Advection, iterated Stokes' and "
                                 "'single Advection, no Stokes'."));

          // Functionality to average the additional RHS terms over the cell is not implemented.
          // This enforces that the variable 'Material averaging' is set to 'none'.
          AssertThrow((this->get_parameters().material_averaging ==
                       MaterialModel::MaterialAveraging::none),
                      ExcMessage("The viscoelastic_plastic material model cannot be used with "
                                 "material averaging. The variable 'Material averaging' "
                                 "in the 'Material model' subsection must be set to 'none'."));
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();



      // Declare dependencies on solution variables
      this->model_dependence.viscosity = NonlinearDependence::compositional_fields;
      this->model_dependence.density = NonlinearDependence::temperature | NonlinearDependence::compositional_fields;
      this->model_dependence.compressibility = NonlinearDependence::none;
      this->model_dependence.specific_heat = NonlinearDependence::compositional_fields;
      this->model_dependence.thermal_conductivity = NonlinearDependence::compositional_fields;
    }



    template <int dim>
    void
    ViscoelasticPlastic<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      if (out.template get_additional_output<ElasticAdditionalOutputs<dim> >() == NULL)
        {
          const unsigned int n_points = out.viscosities.size();
          out.additional_outputs.push_back(
            std_cxx11::shared_ptr<MaterialModel::AdditionalMaterialOutputs<dim> >
            (new MaterialModel::ElasticAdditionalOutputs<dim> (n_points)));
        }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(ViscoelasticPlastic,
                                   "viscoelastic plastic",
                                   "A material model that combines non-linear plasticity with a simple "
                                   "linear viscoelastic material behavior. The model is incompressible. "
                                   "Note that this material model is based heavily on and combines "
                                   "functionality from the following material models: "
                                   "DiffusionDislocation, DruckerPrager, ViscoPlastic and Viscoelastic. "
                                   "\n\n"
                                   "Plasticity limits viscous stress through a Drucker Prager "
                                   "yield criterion, where the yield stress in 3D is  "
                                   "$\\sigma_y = \\frac{6*C*\\cos(\\phi) + 2*P*\\sin(\\phi)} "
                                   "{\\sqrt(3)*(3+\\sin(\\phi))}$ "
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
                                   "The viscoelastic rheology behavior takes into account the elastic shear "
                                   "strength (e.g., shear modulus), while the tensile and volumetric "
                                   "strength (e.g., Young's and bulk modulus) are not considered. The "
                                   "model is incompressible and allows specifying an arbitrary number "
                                   "of compositional fields, where each field represents a different "
                                   "rock type or component of the viscoelastic stress tensor. The stress "
                                   "tensor in 2D and 3D, respectively, contains 3 or 6 components. The "
                                   "compositional fields representing these components must be named "
                                   "and listed in a very specific format, which is designed to minimize"
                                   "mislabeling stress tensor components as distinct 'compositional "
                                   "rock types' (or vice versa). For 2D models, the first three "
                                   "compositional fields must be labeled 'stress_xx', 'stress_yy' and 'stress_xy'. "
                                   "In 3D, the first six compositional fields must be labeled 'stress_xx', "
                                   "'stress_yy', 'stress_zz', 'stress_xy', 'stress_xz', 'stress_yz'. "
                                   "\n\n "
                                   "Combining this viscoelasticity implementation with non-linear viscous flow "
                                   "and plasticity produces a constitutive relationship commonly referred to "
                                   "as partial elastoviscoplastic (e.g., pEVP) in the geodynamics community. "
                                   "While extensively discussed and applied within the geodynamics "
                                   "literature, notable references include: "
                                   "Moresi et al. (2003), J. Comp. Phys., v. 184, p. 476-497. "
                                   "Gerya and Yuen (2007), Phys. Earth. Planet. Inter., v. 163, p. 83-105. "
                                   "Gerya (2010), Introduction to Numerical Geodynamic Modeling. "
                                   "Kaus (2010), Tectonophysics, v. 484, p. 36-47. "
                                   "Choi et al. (2013), J. Geophys. Res., v. 118, p. 2429-2444. "
                                   "Keller et al. (2013), Geophys. J. Int., v. 195, p. 1406-1442. "
                                   "\n\n "
                                   "The overview below directly follows Moresi et al. (2003) eqns. 23-32. "
                                   "However, an important distinction between this material model and "
                                   "the studies above is the use of compositional fields, rather than "
                                   "tracers, to track individual components of the viscoelastic stress "
                                   "tensor. The material model will be udpated when an option to track "
                                   "and calculate viscoelastic stresses with tracers is implemented. "
                                   "\n\n "
                                   "Moresi et al. (2003) begins (eqn. 23) by writing the deviatoric "
                                   "rate of deformation ($\\hat{D}$) as the sum of elastic "
                                   "($\\hat{D_{e}}$) and viscous ($\\hat{D_{v}}$) components: "
                                   "$\\hat{D} = \\hat{D_{e}} + \\hat{D_{v}}$.  "
                                   "These terms further decompose into "
                                   "$\\hat{D_{v}} = \\frac{\\tau}{2\\eta}$ and "
                                   "$\\hat{D_{e}} = \\frac{\\overset{\\triangledown}{\\tau}}{2\\mu}$, where "
                                   "$\\tau$ is the viscous deviatoric stress, $\\eta$ is the shear viscosity, "
                                   "$\\mu$ is the shear modulus and $\\overset{\\triangledown}{\\tau}$ is the "
                                   "Jaumann corotational stress rate. This later term (eqn. 24) contains the "
                                   "time derivative of the deviatoric stress ($\\dot{\\tau}$) and terms that "
                                   "account for material spin (e.g., rotation) due to advection: "
                                   "$\\overset{\\triangledown}{\\tau} = \\dot{\\tau} + {\\tau}W -W\\tau$. "
                                   "Above, $W$ is the material spin tensor (eqn. 25): "
                                   "$W_{ij} = \\frac{1}{2} \\left (\\frac{\\partial V_{i}}{\\partial x_{j}} - "
                                   "\\frac{\\partial V_{j}}{\\partial x_{i}} \\right )$. "
                                   "\n\n "
                                   "The Jaumann stress-rate can also be approximated using terms from the time "
                                   "at the previous time step ($t$) and current time step ($t + \\Delta t_^{e}$): "
                                   "$\\smash[t]{\\overset{\\triangledown}{\\tau}}^{t + \\Delta t^{e}} \\approx "
                                   "\\frac{\\tau^{t + \\Delta t^{e} - \\tau^{t}}}{\\Delta t^{e}} - "
                                   "W^{t}\\tau^{t} + \\tau^{t}W^{t}$. "
                                   "In this material model, the size of the time step above ($\\Delta t^{e}$) "
                                   "can be specified as the numerical time step size or an independent fixed time "
                                   "step. If the latter case is a selected, the user has an option to apply a "
                                   "stress averaging scheme to account for the differences between the numerical "
                                   "and fixed elastic time step (eqn. 32). If one selects to use a fixed elastic time "
                                   "step throughout the model run, this can still be achieved by using CFL and "
                                   "maximum time step values that restrict the numerical time step to a specific time. "
                                   "\n\n "
                                   "The formulation above allows rewriting the total rate of deformation (eqn. 29) as "
                                   "$\\tau^{t + \\Delta t^{e}} = \\eta_{eff} \\left ( "
                                   "2\\hat{D}^{t + \\triangle t^{e}} + \\frac{\\tau^{t}}{\\mu \\Delta t^{e}} + "
                                   "\\frac{W^{t}\\tau^{t} - \\tau^{t}W^{t}}{\\mu}  \\right )$. "
                                   "\n\n "
                                   "The effective viscosity (eqn. 28) is a function of the viscosity ($\\eta$), "
                                   "elastic time step size ($\\Delta t^{e}$) and shear relaxation time "
                                   "($ \\alpha = \\frac{\\eta}{\\mu} $): "
                                   "$\\eta_{eff} = \\eta \\frac{\\Delta t^{e}}{\\Delta t^{e} + \\alpha}$ "
                                   "The magnitude of the shear modulus thus controls how much the effective "
                                   "viscosity is reduced relative to the initial viscosity. "
                                   "\n\n "
                                   "Elastic effects are introduced into the governing stokes equations through "
                                   "an elastic force term (eqn. 30) using stresses from the previous time step: "
                                   "$F^{e,t} = -\\frac{\\eta_{eff}}{\\mu \\Delta t^{e}} \\tau^{t}$. "
                                   "This force term is added onto the right-hand side force vector in the "
                                   "system of equations. "
                                   "\n\n "
                                   "The value of each compositional field representing distinck rock types at a "
                                   "point is interpreted to be a volume fraction of that rock type. If the sum of "
                                   "the compositional field volume fractions is less than one, then the remainder "
                                   "of the volume is assumed to be 'background material'."
                                   "\n\n "
                                   "Several model parameters (densities, elastic shear moduli, thermal expansivities, "
                                   "plasticity parameters, viscosity terms, etc) can be defined per-compositional field. "
                                   "For each material parameter the user supplies a comma delimited list of length "
                                   "N+1, where N is the number of compositional fields. The additional field corresponds "
                                   "to the value for background material. They should be ordered ''background, "
                                   "composition1, composition2...''. However, the first 3 (2D) or 6 (3D) composition "
                                   "fields correspond to components of the elastic stress tensor and their material "
                                   "values will not contribute to the volume fractions. If a single value is given, then "
                                   "all the compositional fields are given that value. Other lengths of lists are not "
                                   "allowed. For a given compositional field the material parameters are treated as "
                                   "constant, except density, which varies linearly with temperature according to the "
                                   "thermal expansivity. "
                                   "\n\n "
                                   "When more than one compositional field is present at a point, they are averaged "
                                   "arithmetically. An exception is viscosity, which may be averaged arithmetically, "
                                   "harmonically, geometrically, or by selecting the viscosity of the composition field "
                                   "with the greatest volume fraction.")
  }
}
