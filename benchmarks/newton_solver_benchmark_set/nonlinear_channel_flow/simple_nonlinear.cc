/*
  Copyright (C) 2011 - 2022 by the authors of the ASPECT code.

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
#include <aspect/simulator.h>
#include <deal.II/grid/tria.h>
#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/newton.h>

#include <iostream>

#include <aspect/utilities.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * A material model based on a simple power law rheology and
     * implementing the derivatives needed for the Newton method.
     *
     * Power law equation to compute the viscosity $\eta$ per composition:
     * $\eta = A^{\frac{-1}{n}} * \left(2.0 * \dot\varepsilon_{II}\right)
     * ^{\frac{1}{n}-1}$ where $A$ is the prefactor, $\dot\varepsion$ is
     * the strain-rate, II indicates the square root of the second invariant
     * defined as $\frac{1}{2} \dot\varepsilon_{ij} \dot\varepsilon_{ij}$,
     * and $n$ is the stress exponent.
     *
     * The viscosities per composition are averaged using the
     * utilities weighted p-norm function. The volume fractions are used
     * as weights for the averaging.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class SimpleNonlinear : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const;

        /**
         * This material model is incompressible.
         */
        virtual bool is_compressible () const;


        static
        void
        declare_parameters (ParameterHandler &prm);

        virtual
        void
        parse_parameters (ParameterHandler &prm);

      private:
        /**
         * For calculating density by thermal expansivity. Units: $\\si{\\kelvin}$
         */
        double reference_temperature;

        /**
         * Defining a minimum strain rate stabilizes the viscosity calculation,
         * which involves a division by the strain rate. Units: $\\si{\\per\\second}$.
         */
        std::vector<double> min_strain_rate;
        std::vector<double> min_viscosity;
        std::vector<double> max_viscosity;

        std::vector<double> thermal_diffusivity;
        std::vector<double> heat_capacity;


        std::vector<double> densities;
        std::vector<double> thermal_expansivities;


        std::vector<double> viscosity_prefactor;
        std::vector<double> stress_exponent;

        unsigned int n_fields;
        /**
         * averaging parameters used for the power exponent
         * in the Utilities::weighted_p_norm_average and
         * Utilities::derivative_of_weighted_p_norm_average.
         */
        double viscosity_averaging_p;

        bool use_deviator_of_strain_rate;
    };
  }
}



using namespace dealii;

namespace aspect
{

  namespace MaterialModel
  {
    template <int dim>
    void
    SimpleNonlinear<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      //set up additional output for the derivatives
      MaterialModelDerivatives<dim> *derivatives = out.template get_additional_output<MaterialModelDerivatives<dim>>();

      for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
        {
          const double temperature = in.temperature[i];

          // Averaging composition-field dependent properties
          // Compositions
          //This assert may be help when designing tests for this material model
          AssertThrow(in.composition[i].size()+1 == n_fields,
                      ExcMessage("Number of compositional fields + 1 not equal to number of fields given in input file."));

          const std::vector<double> volume_fractions = MaterialUtilities::compute_composition_fractions(in.composition[i]);
          double density = 0.0;
          for (unsigned int c=0; c < volume_fractions.size(); ++c)
            {
              //not strictly correct if thermal expansivities are different, since we are interpreting
              //these compositions as volume fractions, but the error introduced should not be too bad.
              const double temperature_factor = (1.0 - thermal_expansivities[c] * (temperature - reference_temperature));
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
          if (in.requests_property(MaterialProperties::viscosity))
            {
              // This function calculates viscosities assuming that all the compositional fields
              // experience the same strain rate (isostrain). Since there is only one process in
              // this material model (a general powerlaw) we do not need to worry about how to
              // distribute the strain-rate and stress over the processes.
              std::vector<double> composition_viscosities(volume_fractions.size());

              std::vector<SymmetricTensor<2,dim>> composition_viscosities_derivatives(volume_fractions.size());

              for (unsigned int c=0; c < volume_fractions.size(); ++c)
                {
                  // If strain rate is zero (like during the first time step) set it to some very small number
                  // to prevent a division-by-zero, and a floating point exception.
                  // Otherwise, calculate the square-root of the norm of the second invariant of the deviatoric-
                  // strain rate (often simplified as epsilondot_ii)
                  const SymmetricTensor<2,dim> edot = use_deviator_of_strain_rate ? deviator(in.strain_rate[i]) : in.strain_rate[i];
                  const double edot_ii_strict = std::sqrt(0.5*edot*edot);
                  const double edot_ii = 2.0 * std::max(edot_ii_strict, min_strain_rate[c]);

                  const double stress_exponent_inv = 1/stress_exponent[c];
                  composition_viscosities[c] = std::max(std::min(std::pow(viscosity_prefactor[c],-stress_exponent_inv) * std::pow(edot_ii,stress_exponent_inv-1), max_viscosity[c]), min_viscosity[c]);
                  Assert(dealii::numbers::is_finite(composition_viscosities[c]), ExcMessage ("Error: Viscosity is not finite."));
                  if (derivatives != NULL)
                    {
                      if (edot_ii_strict > min_strain_rate[c] && composition_viscosities[c] < max_viscosity[c] && composition_viscosities[c] > min_viscosity[c])
                        {
                          composition_viscosities_derivatives[c] = (stress_exponent_inv-1) * composition_viscosities[c] * (1.0/(edot*edot)) * edot;

                          if (use_deviator_of_strain_rate == true)
                            composition_viscosities_derivatives[c] = composition_viscosities_derivatives[c] * deviator_tensor<dim>();
                        }
                      else
                        {
                          composition_viscosities_derivatives[c] = 0;
                        }
                    }
                }

              out.viscosities[i] = Utilities::weighted_p_norm_average(volume_fractions, composition_viscosities, viscosity_averaging_p);
              Assert(dealii::numbers::is_finite(out.viscosities[i]), ExcMessage ("Error: Averaged viscosity is not finite."));

              if (derivatives != NULL)
                {
                  derivatives->viscosity_derivative_wrt_strain_rate[i] = Utilities::derivative_of_weighted_p_norm_average(out.viscosities[i],volume_fractions, composition_viscosities, composition_viscosities_derivatives, viscosity_averaging_p);

                  for (unsigned int x = 0; x < dim; x++)
                    for (unsigned int y = 0; y < dim; y++)
                      Assert (dealii::numbers::is_finite(derivatives->viscosity_derivative_wrt_strain_rate[i][x][y]),
                              ExcMessage ("Averaged viscosity to strain-rate devrivative is not finite."));

                  derivatives->viscosity_derivative_wrt_pressure[i]  = 0;

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
          prm.declare_entry ("Reference temperature", "293.", Patterns::Double(0.),
                             "For calculating density by thermal expansivity. Units: \\si{\\kelvin}");
          prm.declare_entry ("Minimum strain rate", "1.96e-40", Patterns::List(Patterns::Double(0.)),
                             "Stabilizes strain dependent viscosity. Units: \\si{\\per\\second}");
          prm.declare_entry ("Minimum viscosity", "1e10", Patterns::List(Patterns::Double(0.)),
                             "Lower cutoff for effective viscosity. Units: \\si{\\pascal\\second}");
          prm.declare_entry ("Maximum viscosity", "1e28", Patterns::List(Patterns::Double(0.)),
                             "Upper cutoff for effective viscosity. Units: \\si{\\pascal\\second}");
          prm.declare_entry ("Effective viscosity coefficient", "1.0", Patterns::List(Patterns::Double(0.)),
                             "Scaling coefficient for effective viscosity.");

          // Equation of state parameters
          prm.declare_entry ("Thermal diffusivity", "0.8e-6", Patterns::List(Patterns::Double(0.)),
                             "Units: \\si{\\meter\\squared\\per\\second}");
          prm.declare_entry ("Heat capacity", "1.25e3", Patterns::List(Patterns::Double(0.)),
                             "Units: \\si{\\joule\\per\\kelvin\\per\\kilogram}");
          prm.declare_entry ("Densities", "3300.",
                             Patterns::List(Patterns::Double(0.)),
                             "List of densities, $\\rho$, for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value. "
                             "Units: \\si{\\kilogram\\per\\meter\\cubed}");
          prm.declare_entry ("Thermal expansivities", "3.5e-5",
                             Patterns::List(Patterns::Double(0.)),
                             "List of thermal expansivities for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value.  Units: \\si{\\per\\kelvin}");


          // SimpleNonlinear creep parameters
          prm.declare_entry ("Viscosity prefactor", "1e-37",
                             Patterns::List(Patterns::Double(0)),
                             "List of viscosity prefactors, $A$, for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value. "
                             "Units: \\si{\\pascal}$^{-n_{\\text{dislocation}}}$ \\si{\\meter}$^{n_{\\text{dislocation}}/m_{\\text{dislocation}}}$ \\si{\\per\\second}");
          prm.declare_entry ("Stress exponent", "3.",
                             Patterns::List(Patterns::Double(0.)),
                             "List of stress exponents, $n_dislocation$, for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value.  Units: None");

          // averaging parameters
          prm.declare_entry ("Viscosity averaging p", "-1.",
                             Patterns::Double(),
                             "This is the p value in the generalized weighed average equation: "
                             " $\\text{mean} = \\frac{1}{k}(\\sum_{i=1}^k \\big(c_i \\eta_{\\text{eff}_i}^p)\\big)^{\\frac{1}{p}}$. "
                             " Units: \\si{\\pascal\\second}");

          // strain-rate deviator parameter
          prm.declare_entry ("Use deviator of strain-rate", "true",
                             Patterns::Bool(),
                             "This value determines wheter to use the deviator of the strain-rate in computing the viscosity, "
                             "or simply the strain rate $\\varepsilon(\\mathbf u)$.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    SimpleNonlinear<dim>::parse_parameters (ParameterHandler &prm)
    {
      // Can't use this->n_compositional_fields(), because some
      // tests never initialize the simulator, but uses the material
      // model directly.
      prm.enter_subsection ("Compositional fields");
      {
        n_fields = prm.get_integer ("Number of fields")+1;
      }
      prm.leave_subsection();

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection ("Simple nonlinear");
        {

          // Reference and minimum/maximum values
          reference_temperature = prm.get_double("Reference temperature");
          min_strain_rate =  Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Minimum strain rate"))),n_fields,"Minimum strain rate");
          min_viscosity =  Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get ("Minimum viscosity"))),n_fields,"Minimum viscosity");
          max_viscosity =  Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get ("Maximum viscosity"))),n_fields,"Maximum viscosity");

          // Equation of state parameters
          thermal_diffusivity =  Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Thermal diffusivity"))),n_fields,"Thermal diffusivity");
          heat_capacity =  Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Heat capacity"))),n_fields,"Heat capacity");

          // Compositional parameters
          densities =  Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Densities"))),n_fields,"Densities");
          thermal_expansivities =  Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Thermal expansivities"))),n_fields,"Thermal expansivities");

          // Rheological parameters
          viscosity_prefactor =  Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Viscosity prefactor"))),n_fields,"Viscosity prefactor");
          stress_exponent =  Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Stress exponent"))),n_fields,"Stress exponent");

          // averaging parameters
          viscosity_averaging_p = prm.get_double("Viscosity averaging p");

          use_deviator_of_strain_rate = prm.get_bool ("Use deviator of strain-rate");

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
                                   "A material model based on a simple power law rheology and implementing the derivatives "
                                   "needed for the Newton method. "
                                   "Power law equation to compute the viscosity $\\eta$ per composition comes from the book Introduction "
                                   "to Geodynamic Modelling by Taras Gerya (2010): "
                                   "$\\eta = A^{\\frac{-1}{n}} * \\left(2.0 * \\dot\\varepsilon_{II}\\right)^{\\frac{1}{n}-1}$ "
                                   "where $A$ is the prefactor, $\\dot\\varepsion$ is the strain-rate, II indicates the square "
                                   "root of the second invariant defined as $\\frac{1}{2} \\dot\\varepsilon_{ij} \\dot\\varepsilon_{ij}$, "
                                   "and $n$ is the stress exponent. "
                                   "The viscosities per composition are averaged using the utilities weighted "
                                   "p-norm function. The volume fractions are used as weights for the averaging. ")
  }
}
