/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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

#include <aspect/material_model/simple_nonlinear.h>

using namespace dealii;

namespace aspect
{
  namespace
  {
    std::vector<double>
    get_vector_double (const std::string &parameter, const unsigned int n_fields, ParameterHandler &prm)
    {
      std::vector<double> parameter_list;
      parameter_list = Utilities::string_to_double(Utilities::split_string_list(prm.get (parameter)));
      if (parameter_list.size() == 1)
        parameter_list.resize(n_fields, parameter_list[0]);

      AssertThrow(parameter_list.size() == n_fields,
                  ExcMessage("Length of " + parameter + " list (size "+ std::to_string(parameter_list.size()) +") must be either one,"
                             " or n_compositional_fields+1 (= " + std::to_string(n_fields) + ")."));

      return parameter_list;
    }
  }

  namespace MaterialModel
  {
    template <int dim>
    std::vector<double>
    SimpleNonlinear<dim>::
    compute_volume_fractions( const std::vector<double> &compositional_fields) const
    {
      /*std::cout << "vfFlag 1 " << std::endl;
         std::vector<double> volume_fractions(compositional_fields.size());
         std::cout << "vfFlag 2 " << std::endl;
         //clip the compositional fields so they are between zero and one
         std::vector<double> x_comp = compositional_fields;
         for ( unsigned int i=0; i < x_comp.size(); ++i)
           x_comp[i] = std::min(std::max(x_comp[i], 0.0), 1.0);
         std::cout << "vfFlag 5 " << std::endl;
         //sum the compositional fields for normalization purposes
         double sum_composition = 0.0;
         for ( unsigned int i=0; i < x_comp.size(); ++i)
           sum_composition += x_comp[i];
         std::cout << "vfFlag 8 " << std::endl;
         // If the sum of the compositions is smaller
         // than 0.1 something is going really wrong
         // and we abort.
         AssertThrow(sum_composition >= 0.1,
             ExcMessage("The sum of the compositions is smaller than 0.1. There are " + std::to_string(x_comp.size()) + " compositions. It is " +  std::to_string(sum_composition) + ", xcomp 0: "  +  std::to_string(x_comp[0]) + ", comp 0: "  +  std::to_string(compositional_fields[0]) ))
         std::cout << "vfFlag 10 " << std::endl;
         //volume_fractions[0] = 0.0;  //background mantle
         for ( unsigned int i=0; i < x_comp.size(); ++i)
           volume_fractions[i] = x_comp[i]/sum_composition;std::cout << "vfFlag 15 " << std::endl;
         //std::cout << "Multiple compositions, return them." << std::endl;
         return volume_fractions;*/

      std::vector<double> volume_fractions( compositional_fields.size()+1);

      //clip the compositional fields so they are between zero and one
      std::vector<double> x_comp = compositional_fields;
      for ( unsigned int i=0; i < x_comp.size(); ++i)
        x_comp[i] = std::min(std::max(x_comp[i], 0.0), 1.0);

      //sum the compositional fields for normalization purposes
      double sum_composition = 0.0;
      for ( unsigned int i=0; i < x_comp.size(); ++i)
        sum_composition += x_comp[i];

      if (sum_composition >= 1.0)
        {
          volume_fractions[0] = 0.0;  //background mantle
          for ( unsigned int i=1; i <= x_comp.size(); ++i)
            volume_fractions[i] = x_comp[i-1]/sum_composition;
        }
      else
        {
          volume_fractions[0] = 1.0 - sum_composition; //background mantle
          for ( unsigned int i=1; i <= x_comp.size(); ++i)
            volume_fractions[i] = x_comp[i-1];
        } //std::cout << "Background mantle = " << volume_fractions[0] << std::endl;
      return volume_fractions;
    }

    template <int dim>
    double
    SimpleNonlinear<dim>::
    p_norm_average ( const std::vector<double> &composition,
                     const std::vector<double> &parameter_values,
                     const double p) const
    {
      double averaged_parameter = 0.0;
      const std::vector<double> volume_fractions = compute_volume_fractions(composition);

      // first look at the specal case which can be done faster
      // Minimum
      if (p <= -1000)
        {
          const unsigned int i = (unsigned int)(std::min_element( volume_fractions.begin(),
                                                                  volume_fractions.end() )
                                                - volume_fractions.begin());
          averaged_parameter = parameter_values[i];
          return averaged_parameter;
        }
      // Harmonic average (TODO: Check if we need to use a epsilon)
      if (p == -1)
        {
          for (unsigned int i=0; i< volume_fractions.size(); ++i)
            averaged_parameter += volume_fractions[i]/(parameter_values[i]);
          averaged_parameter = 1.0/averaged_parameter;
          return averaged_parameter;
        }

      // Geometric average (TODO: Check if we need to use a epsilon)
      if (p == 0)
        {
          for (unsigned int i=0; i < volume_fractions.size(); ++i)
            averaged_parameter += volume_fractions[i]*std::log(parameter_values[i]);
          averaged_parameter = std::exp(averaged_parameter);
          return averaged_parameter;
        }

      // Arithmetic average (TODO: Check if we need to use a epsilon)
      if (p == 1)
        {
          for (unsigned int i=0; i< volume_fractions.size(); ++i)
            averaged_parameter += volume_fractions[i]*parameter_values[i];
          return averaged_parameter;
        }

      // Quadratic average (RMS) (TODO: Check if we need to use a epsilon)
      if (p == 2)
        {
          for (unsigned int i=0; i< volume_fractions.size(); ++i)
            averaged_parameter += volume_fractions[i]*parameter_values[i]*parameter_values[i];
          averaged_parameter = std::sqrt(averaged_parameter);
          return averaged_parameter;
        }

      // Cubic average (TODO: Check if we need to use a epsilon)
      if (p == 3)
        {
          for (unsigned int i=0; i< volume_fractions.size(); ++i)
            averaged_parameter += volume_fractions[i]*parameter_values[i]*parameter_values[i]*parameter_values[i];
          averaged_parameter = std::cbrt(averaged_parameter);
          return averaged_parameter;
        }

      // Maximum
      if (p >= 1000)
        {
          const unsigned int i = (unsigned int)(std::max_element( volume_fractions.begin(),
                                                                  volume_fractions.end() )
                                                - volume_fractions.begin());
          averaged_parameter = parameter_values[i];
          return averaged_parameter;
        }

      for (unsigned int i=0; i< volume_fractions.size(); ++i)//{
        averaged_parameter += volume_fractions[i]*std::pow(parameter_values[i],p);
      //std::cout << "ap = " << averaged_parameter << " = " << volume_fractions[i] << "*std::pow(" << parameter_values[i] << "," << p << "), i = " << i << std::endl;}
      averaged_parameter = std::pow(averaged_parameter,1/p); //std::cout << "ap = " << averaged_parameter << std::endl;

      return averaged_parameter;
    }


    template <int dim>
    SymmetricTensor<2,dim>
    SimpleNonlinear<dim>::
    p_norm_average_derivatives (const double averaged_parameter,
                                const std::vector<double> &composition,
                                const std::vector<double> &parameter_values,
                                const std::vector<SymmetricTensor<2,dim> > &parameter_derivatives,
                                const double p) const
    {
      const std::vector<double> volume_fractions = compute_volume_fractions(composition);

      SymmetricTensor<2,dim> averaged_parameter_derivative;
      double averaged_parameter_derivative_part_1 = 0.0;
      SymmetricTensor<2,dim> averaged_parameter_derivative_part_2;
      std::vector<double> parameter_values_p(volume_fractions.size());

      double averaged_parameter_check = 0;

      for (unsigned int i=0; i< volume_fractions.size(); ++i)
        {
          parameter_values_p[i] = std::pow(parameter_values[i],p);
          //std::cout << parameter_values_p[i] << " = std::pow("<< parameter_values[i] << "," << p << ")" << std::endl;
        }

      //const std::vector<double> volume_fractions = compute_volume_fractions(composition);

      // first look at the specal case which can be done faster
      // Minimum
      if (p <= -1000)
        {
          const unsigned int i = (unsigned int)(std::min_element( volume_fractions.begin(),
                                                                  volume_fractions.end() )
                                                - volume_fractions.begin());
          averaged_parameter_derivative = parameter_derivatives[i];
          return averaged_parameter_derivative;
        }
      // Harmonic average
      // Seems to work without a epsilon. TODO: Might be possible to optimalize further.
      if (p == -1)
        {
          for (unsigned int i=0; i< volume_fractions.size(); ++i)
            {
              averaged_parameter_derivative_part_1 += volume_fractions[i] * parameter_values_p[i];
              averaged_parameter_derivative_part_2 += volume_fractions[i]*(1/(parameter_values[i]*parameter_values[i]))*parameter_derivatives[i];
            }
          averaged_parameter_derivative = std::pow(averaged_parameter_derivative_part_1,-2) * averaged_parameter_derivative_part_2;
          return averaged_parameter_derivative;
        }

      // Geometric average
      // Seems to work without a epsilon. TODO: Might be possible to optimalize further.
      if (p == 0)
        {
          for (unsigned int i=0; i < volume_fractions.size(); ++i)
            {
              averaged_parameter_derivative_part_1 += volume_fractions[i]*std::log(parameter_values[i]);
              averaged_parameter_derivative_part_2 += volume_fractions[i]*(1/parameter_values[i])*parameter_derivatives[i];
            }
          averaged_parameter_derivative = std::exp(averaged_parameter_derivative_part_1) * averaged_parameter_derivative_part_2;
          return averaged_parameter_derivative;
        }

      // Arithmetic average
      // Seems to work without a epsilon.
      if (p == 1)
        {
          for (unsigned int i=0; i< volume_fractions.size(); ++i)
            averaged_parameter_derivative_part_2 += volume_fractions[i]*parameter_derivatives[i];
          return averaged_parameter_derivative_part_2;
        }

      // Maximum
      if (p >= 1000)
        {
          const unsigned int i = (unsigned int)(std::max_element( volume_fractions.begin(),
                                                                  volume_fractions.end() )
                                                - volume_fractions.begin());
          averaged_parameter_derivative = parameter_derivatives[i];
          return averaged_parameter_derivative;
        }

      //TODO: This can probably be optimized by using:
      //averaged_parameter_derivative_part_2 += volume_fractions[i]*parameter_values_p[i]*(1/parameter_values[i])*parameter_derivatives[i]; and
      //averaged_parameter_derivative = averaged_parameter * (1/averaged_parameter_derivative_part_1) * averaged_parameter_derivative_part_2;
      for (unsigned int i=0; i< volume_fractions.size(); ++i)
        {
          averaged_parameter_derivative_part_1 += volume_fractions[i] * parameter_values_p[i];
          averaged_parameter_derivative_part_2 += volume_fractions[i]*std::pow(parameter_values[i],p-1)*parameter_derivatives[i];
        }
      averaged_parameter_derivative = std::pow(averaged_parameter_derivative_part_1,(1/p)-1) * averaged_parameter_derivative_part_2;

      return averaged_parameter_derivative;
    }

    template <int dim>
    void
    SimpleNonlinear<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      //set up additional output for the derivatives
      MaterialModelDerivatives<dim> *derivatives;
      derivatives = out.template get_additional_output<MaterialModelDerivatives<dim> >();

      for (unsigned int i=0; i < in.temperature.size(); ++i)
        {
          // const Point<dim> position = in.position[i];
          const double temperature = in.temperature[i];
          const double pressure= in.pressure[i];

          // Averaging composition-field dependent properties
          // Compositions
          //This assert may be help when designing tests for this material model
          AssertThrow(in.composition[i].size()+1 == n_fields,
                      ExcMessage("Number of compositional fields + 1 not equal to number of fields given in input file."));

          const std::vector<double> volume_fractions = compute_volume_fractions(in.composition[i]);
          double density = 0.0;
          for (unsigned int c=0; c < volume_fractions.size(); ++c)
            {
              //not strictly correct if thermal expansivities are different, since we are interpreting
              //these compositions as volume fractions, but the error introduced should not be too bad.
              const double temperature_factor = (1.0 - thermal_expansivities[c] * (temperature - reference_T));
              density += volume_fractions[c] * densities[c] * temperature_factor;
            }

          // thermal expansivities
          double thermal_expansivity = 0.0;
          for (unsigned int c=0; c < volume_fractions.size(); ++c)
            thermal_expansivity += volume_fractions[c] * thermal_expansivities[c];

          // Specific heat at the given positions.
          double specific_heat = 0.0;
          for (unsigned int c=0; c < volume_fractions.size(); ++c)
            specific_heat += volume_fractions[c] * heat_capacity[c];

          // Thermal conductivity at the given positions.
          double thermal_conductivities = 0.0;
          for (unsigned int c=0; c < volume_fractions.size(); ++c)
            thermal_conductivities += volume_fractions[c] * thermal_diffusivity[c] * heat_capacity[c] * densities[c];

          // calculate effective viscosity
          if (in.strain_rate.size())
            {
              // This function calculates viscosities assuming that all the compositional fields
              // experience the same strain rate (isostrain). Since there is only one process in
              // this material model (a general powerlaw) we do not need to worry about how to
              // distribute the strain-rate and stress over the processes.
              std::vector<double> composition_viscosities(volume_fractions.size());
              std::vector<SymmetricTensor<2,dim> > composition_viscosities_derivatives(volume_fractions.size());

              for (unsigned int c=0; c < volume_fractions.size(); ++c)
                {
                  // If strain rate is zero (like during the first time step) set it to some very small number
                  // to prevent a division-by-zero, and a floating point exception.
                  // Otherwise, calculate the square-root of the norm of the second invariant of the deviatoric-
                  // strain rate (often simplified as epsilondot_ii)
                  const double edot_ii_strict = std::sqrt(0.5*deviator(in.strain_rate[i])*deviator(in.strain_rate[i]));
                  const double edot_ii = 2.0 * std::max(edot_ii_strict, min_strain_rate[c] * min_strain_rate[c]);

                  // Find effective viscosities for each of the individual phases
                  // Viscosities should have same number of entries as compositional fields

                  // Power law equation
                  // edot_ii_i = A_i * stress_ii_i^{n_i} * d^{-m} \exp\left(-\frac{E_i^* + PV_i^*}{n_iRT}\right)
                  // where ii indicates the square root of the second invariant and
                  // i corresponds to diffusion or dislocation creep

                  // The isostrain condition implies that the viscosity averaging should be arithmetic (see above).
                  // We have given the user freedom to apply alternative bounds, because in diffusion-dominated
                  // creep (where n_diff=1) viscosities are stress and strain-rate independent, so the calculation
                  // of compositional field viscosities is consistent with any averaging scheme.
                  // TODO: Mailed Bob and asked why the averaging should be arithmetic. Bob replyed that this has to
                  //       do with effective medium theory. Have to look into this a bit more.
                  const double stress_exponent_inv = 1/stress_exponent[c];
                  composition_viscosities[c] = std::max(std::min(std::pow(prefactor[c],-stress_exponent_inv) * std::pow(edot_ii,stress_exponent_inv-1), max_visc[c]), min_visc[c]);
                  Assert(dealii::numbers::is_finite(composition_viscosities[c]),ExcMessage ("Error: Viscosity is not finite."));
                  if (derivatives != NULL)
                    {
                      if (edot_ii_strict > min_strain_rate[c] * min_strain_rate[c] && composition_viscosities[c] < max_visc[c] && composition_viscosities[c] > min_visc[c])
                        {
                          //strictly speaking the derivative is this: 0.5 * ((1/stress_exponent)-1) * std::pow(2,2) * out.viscosities[i] * (1/(edot_ii*edot_ii)) * deviator(in.strain_rate[i])
                          composition_viscosities_derivatives[c] = 2 * (stress_exponent_inv-1) * composition_viscosities[c] * (1/(edot_ii*edot_ii)) * deviator(in.strain_rate[i]);
                          //std::cout << composition_viscosities_derivatives[i] << ", dev strt:" << deviator(in.strain_rate[i]) << std::endl;
                        }
                      else
                        {
                          composition_viscosities_derivatives[c] = 0;
                          //std::cout << "Set " << i << " to 0" << std::endl;
                        }
                    }
                }
              out.viscosities[i] = p_norm_average(in.composition[i], composition_viscosities, viscosity_averaging_p);
              Assert(dealii::numbers::is_finite(out.viscosities[i]),ExcMessage ("Error: Averaged viscosity is not finite."));

              if (derivatives != NULL)
                {
                  derivatives->dviscosities_dstrain_rate[i] = p_norm_average_derivatives(out.viscosities[i],in.composition[i], composition_viscosities, composition_viscosities_derivatives, viscosity_averaging_p);

#ifdef DEBUG
                  for (int x = 0; x < dim; x++)
                    for (int y = 0; y < dim; y++)
                      if (!dealii::numbers::is_finite(derivatives->dviscosities_dstrain_rate[i][x][y]))
                        std::cout << "Error: Averaged viscosity to strain-rate devrivative is not finite." << std::endl;

                  derivatives->dviscosities_dpressure[i]    = 0;
                  if (!dealii::numbers::is_finite(derivatives->dviscosities_dpressure[i]))
                    std::cout << "Error: Averaged viscosity to pressure devrivative is not finite." << std::endl;
#endif
                }
            }
          out.densities[i] = density;
          out.thermal_expansion_coefficients[i] = thermal_expansivity;
          // Specific heat at the given positions.
          out.specific_heat[i] = specific_heat;
          // Thermal conductivity at the given positions.
          out.thermal_conductivities[i] = thermal_conductivities;
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
          for (unsigned int c=0; c < in.composition[i].size(); ++c)
            out.reaction_terms[i][c] = 0.0;
        }
    }
    template <int dim>
    double
    SimpleNonlinear<dim>::
    reference_viscosity () const
    {
      return ref_visc;
    }

    template <int dim>
    double
    SimpleNonlinear<dim>::
    reference_density () const
    {
      return densities[0];
    }

    template <int dim>
    bool
    SimpleNonlinear<dim>::
    is_compressible () const
    {
      return false;
    }

    template <int dim>
    void
    SimpleNonlinear<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Compositional fields");
      {
        prm.declare_entry ("Number of fields", "0",
                           Patterns::Integer (0),
                           "The number of fields that will be advected along with the flow field, excluding "
                           "velocity, pressure and temperature.");
      }
      prm.leave_subsection();
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection ("Simple nonlinear");
        {
          // Reference and minimum/maximum values
          prm.declare_entry ("Reference temperature", "293", Patterns::Double(0),
                             "For calculating density by thermal expansivity. Units: $K$");
          prm.declare_entry ("Minimum strain rate", "1.4e-20", Patterns::List(Patterns::Double(0)),
                             "Stabilizes strain dependent viscosity. Units: $1 / s$");
          prm.declare_entry ("Minimum viscosity", "1e10", Patterns::List(Patterns::Double(0)),
                             "Lower cutoff for effective viscosity. Units: $Pa s$");
          prm.declare_entry ("Maximum viscosity", "1e28", Patterns::List(Patterns::Double(0)),
                             "Upper cutoff for effective viscosity. Units: $Pa s$");
          prm.declare_entry ("Effective viscosity coefficient", "1.0", Patterns::List(Patterns::Double(0)),
                             "Scaling coefficient for effective viscosity.");
          prm.declare_entry ("Reference viscosity", "1e22", Patterns::List(Patterns::Double(0)),
                             "Reference viscosity for nondimensionalization. Units $Pa s$");

          // Equation of state parameters
          prm.declare_entry ("Thermal diffusivity", "0.8e-6", Patterns::List(Patterns::Double(0)), "Units: $m^2/s$");
          prm.declare_entry ("Heat capacity", "1.25e3", Patterns::List(Patterns::Double(0)), "Units: $J / (K * kg)$");
          prm.declare_entry ("Densities", "3300.",
                             Patterns::List(Patterns::Double(0)),
                             "List of densities, $\\rho$, for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one values is given, then all use the same value.  Units: $kg / m^3$");
          prm.declare_entry ("Thermal expansivities", "3.5e-5",
                             Patterns::List(Patterns::Double(0)),
                             "List of thermal expansivities for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one values is given, then all use the same value.  Units: $1 / K$");


          // SimpleNonlinear creep parameters
          prm.declare_entry ("Prefactor", "1e-37",
                             Patterns::List(Patterns::Double(0)),
                             "List of viscosity prefactors, $A$, for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one values is given, then all use the same value. "
                             "Units: $Pa^{-n_{dislocation}} m^{n_{dislocation}/m_{dislocation}} s^{-1}$");
          prm.declare_entry ("Stress exponent", "3",
                             Patterns::List(Patterns::Double(0)),
                             "List of stress exponents, $n_dislocation$, for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one values is given, then all use the same value.  Units: None");

          // averaging parameters
          prm.declare_entry ("Viscosity averaging p", "-1",
                             Patterns::Double(),
                             "This is the p value in the generalized weighed average eqation: "
                             " mean = \\frac{1}{k}(\\sum_{i=1}^k \\big(c_i \\eta_{\\text{eff}_i}^p)\\big)^{\\frac{1}{p}}. "
                             " Units: $Pa s$");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    SimpleNonlinear<dim>::parse_parameters (ParameterHandler &prm)
    {
      // can't use this->n_compositional_fields(), because some
      // tests never initiate the simulator, but uses the material
      // model directly.
      prm.enter_subsection ("Compositional fields");
      {
        n_fields = prm.get_integer ("Number of fields")+1;
        //AssertThrow(n_fields > 0, ExcMessage("This material model needs at least one compositional field."))
      }
      prm.leave_subsection();

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection ("Simple nonlinear");
        {

          // Reference and minimum/maximum values
          reference_T = prm.get_double("Reference temperature");
          ref_visc = prm.get_double ("Reference viscosity");
          min_strain_rate = get_vector_double("Minimum strain rate",n_fields,prm);
          min_visc = get_vector_double ("Minimum viscosity",n_fields,prm);
          max_visc = get_vector_double ("Maximum viscosity",n_fields,prm);
          veff_coefficient = get_vector_double ("Effective viscosity coefficient",n_fields,prm);

          // Equation of state parameters
          thermal_diffusivity = get_vector_double("Thermal diffusivity",n_fields,prm);
          heat_capacity = get_vector_double("Heat capacity",n_fields,prm);

          // ---- Compositional parameters
          densities = get_vector_double("Densities",n_fields,prm);
          thermal_expansivities = get_vector_double("Thermal expansivities",n_fields,prm);

          // Rheological parameters
          prefactor = get_vector_double("Prefactor",n_fields,prm);
          stress_exponent = get_vector_double("Stress exponent",n_fields,prm);

          // averaging parameters
          viscosity_averaging_p = prm.get_double("Viscosity averaging p");

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
    ASPECT_REGISTER_MATERIAL_MODEL(SimpleNonlinear,
                                   "simple nonlinear",
                                   " An implementation of a viscous rheology including diffusion"
                                   " and dislocation creep."
                                   " Compositional fields can each be assigned individual"
                                   " activation energies, reference densities, thermal expansivities,"
                                   " and stress exponents. The effective viscosity is defined as"
                                   " \n\n"
                                   " \\[v_\\text{eff} = \\left(\\frac{1}{v_\\text{eff}^\\text{diff}}+"
                                   " \\frac{1}{v_\\text{eff}^\\text{dis}}\\right)^{-1}\\]"
                                   " where"
                                   " \\[v_\\text{i} = 0.5 * A^{-\\frac{1}{n_i}} d^\\frac{m_i}{n_i}"
                                   " \\dot{\\varepsilon_i}^{\\frac{1-n_i}{n_i}}"
                                   " \\exp\\left(\\frac{E_i^* + PV_i^*}{n_iRT}\\right)\\]"
                                   " \n\n"
                                   " where $d$ is grain size, $i$ corresponds to diffusion or dislocation creep,"
                                   " $\\dot{\\varepsilon}$ is the square root of the second invariant of the"
                                   " strain rate tensor, $R$ is the gas constant, $T$ is temperature, "
                                   " and $P$ is pressure."
                                   " $A_i$ are prefactors, $n_i$ and $m_i$ are stress and grain size exponents"
                                   " $E_i$ are the activation energies and $V_i$ are the activation volumes."
                                   " \n\n"
                                   " The ratio of diffusion to dislocation strain rate is found by Newton's"
                                   " method, iterating to find the stress which satisfies the above equations."
                                   " The value for the components of this formula and additional"
                                   " parameters are read from the parameter file in subsection"
                                   " 'Material model/SimpleNonlinear'.")
  }
}
