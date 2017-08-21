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
     * $\eta = A * \dot\varepsilon_{II}^{\frac{1}{n}-1}$ where $A$ is the
     * prefactor, $\dot\varepsion$ is the strain-rate, II indicates
     * the square root of the second invariant defined as
     * $\frac{1}{2} \dot\varepsilon_{ij} \dot\varepsilon_{ij}$, and
     * $n$ is the stress exponent.
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

        virtual double reference_viscosity () const;


        static
        void
        declare_parameters (ParameterHandler &prm);

        virtual
        void
        parse_parameters (ParameterHandler &prm);

      private:
        std::vector<double> compute_volume_fractions( const std::vector<double> &compositional_fields) const;


        /**
         * For calculating density by thermal expansivity. Units: $K$
         */
        double reference_temperature;

        /**
         * Defining a minimum strain rate stabilizes the viscosity calculation,
         * which involves a division by the strain rate. Units: $1/s$.
         */
        std::vector<double> min_strain_rate;
        std::vector<double> min_viscosity;
        std::vector<double> max_viscosity;
        // the reference viscosity returned by the reference_viscosity function.
        double ref_viscosity;

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
    };
  }
}



using namespace dealii;

namespace aspect
{

  namespace MaterialModel
  {
    template <int dim>
    std::vector<double>
    SimpleNonlinear<dim>::
    compute_volume_fractions( const std::vector<double> &compositional_fields) const
    {
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
        }
      return volume_fractions;
    }

    template <int dim>
    void
    SimpleNonlinear<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      //set up additional output for the derivatives
      MaterialModelDerivatives<dim> *derivatives = out.template get_additional_output<MaterialModelDerivatives<dim> >();

      for (unsigned int i=0; i < in.temperature.size(); ++i)
        {
          const double temperature = in.temperature[i];

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

                  const double stress_exponent_inv = 1/stress_exponent[c];
                  composition_viscosities[c] = std::max(std::min(std::pow(viscosity_prefactor[c],-stress_exponent_inv) * std::pow(edot_ii,stress_exponent_inv-1), max_viscosity[c]), min_viscosity[c]);
                  Assert(dealii::numbers::is_finite(composition_viscosities[c]), ExcMessage ("Error: Viscosity is not finite."));
                  if (derivatives != NULL)
                    {
                      if (edot_ii_strict > min_strain_rate[c] * min_strain_rate[c] && composition_viscosities[c] < max_viscosity[c] && composition_viscosities[c] > min_viscosity[c])
                        {
                          /**
                           * strictly speaking the derivative is this:
                           * 0.5 * ((1/stress_exponent)-1) * std::pow(2,2) * out.viscosities[i] * (1/(edot_ii*edot_ii)) * deviator(in.strain_rate[i]),
                           * but we can (and do) simply it as:
                           */
                          composition_viscosities_derivatives[c] = 2.0 * (stress_exponent_inv-1) * composition_viscosities[c] * (1.0/(edot_ii*edot_ii)) * deviator(in.strain_rate[i]);
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
    double
    SimpleNonlinear<dim>::
    reference_viscosity () const
    {
      return ref_viscosity;
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
          prm.declare_entry ("Viscosity prefactor", "1e-37",
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
                             "This is the p value in the generalized weighed average equation: "
                             " $\\text{mean} = \\frac{1}{k}(\\sum_{i=1}^k \\big(c_i \\eta_{\\text{eff}_i}^p)\\big)^{\\frac{1}{p}}$. "
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
          ref_viscosity = prm.get_double ("Reference viscosity");
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
                                   "Power law equation to compute the viscosity $\\eta$ per composition: "
                                   "$\\eta = A * \\dot\\varepsilon_{II}^{\\frac{1}{n}-1}$ where $A$ is the prefactor, "
                                   "$\\dot\\varepsion$ is the strain-rate, II indicates the square root of the second "
                                   "invariant defined as $\\frac{1}{2} \\dot\\varepsilon_{ij} \\dot\\varepsilon_{ij}$, and "
                                   "$n$ is the stress exponent. "
                                   "The viscosities per composition are averaged using the utilities weighted "
                                   "p-norm function. The volume fractions are used as weights for the averaging. ")
  }
}



int f(double parameter)
{

  std::cout << std::endl << "Test for p = " << parameter << std::endl;

  const int dim=2;
  using namespace aspect::MaterialModel;
  MaterialModelInputs<dim> in_base(5,3);
  in_base.composition[0][0] = 0;
  in_base.composition[0][1] = 0;
  in_base.composition[0][2] = 0;
  in_base.composition[1][0] = 0.75;
  in_base.composition[1][1] = 0.15;
  in_base.composition[1][2] = 0.10;
  in_base.composition[2][0] = 0;
  in_base.composition[2][1] = 0.2;
  in_base.composition[2][2] = 0.4;
  in_base.composition[3][0] = 0;
  in_base.composition[3][1] = 0.2;
  in_base.composition[3][2] = 0.4;
  in_base.composition[4][0] = 1;
  in_base.composition[4][1] = 0;
  in_base.composition[4][2] = 0;

  in_base.pressure[0] = 1e9;
  in_base.pressure[1] = 5e9;
  in_base.pressure[2] = 2e10;
  in_base.pressure[3] = 2e11;
  in_base.pressure[4] = 2e12;

  /**
   * We can't take to small strain-rates, because then the difference in the
   * viscosity will be too small for the double accuracy which stores
   * the viscosity solutions and the finite difference solution.
   */
  in_base.strain_rate[0] = SymmetricTensor<2,dim>();
  in_base.strain_rate[0][0][0] = 1e-12;
  in_base.strain_rate[0][0][1] = 1e-12;
  in_base.strain_rate[0][1][1] = 1e-11;
  in_base.strain_rate[1] = SymmetricTensor<2,dim>(in_base.strain_rate[0]);
  in_base.strain_rate[1][0][0] = -1.71266e-13;
  in_base.strain_rate[1][0][1] = -5.82647e-12;
  in_base.strain_rate[1][1][1] = 4.21668e-14;
  in_base.strain_rate[2] = SymmetricTensor<2,dim>(in_base.strain_rate[0]);
  in_base.strain_rate[2][1][1] = 1e-13;
  in_base.strain_rate[2][0][1] = 1e-11;
  in_base.strain_rate[2][0][0] = -1e-12;
  in_base.strain_rate[3] = SymmetricTensor<2,dim>(in_base.strain_rate[0]);
  in_base.strain_rate[3][1][1] = 4.9e-21;
  in_base.strain_rate[3][0][1] = 4.9e-21;
  in_base.strain_rate[3][0][0] = 4.9e-21;
  in_base.strain_rate[4] = SymmetricTensor<2,dim>(in_base.strain_rate[0]);
  in_base.strain_rate[4][1][1] = 1e-11;
  in_base.strain_rate[4][0][1] = 1e-11;
  in_base.strain_rate[4][0][0] = 1e-11;

  in_base.temperature[0] = 293;
  in_base.temperature[1] = 1600;
  in_base.temperature[2] = 2000;
  in_base.temperature[3] = 2100;
  in_base.temperature[4] = 2200;

  SymmetricTensor<2,dim> zerozero = SymmetricTensor<2,dim>();
  SymmetricTensor<2,dim> onezero = SymmetricTensor<2,dim>();
  SymmetricTensor<2,dim> oneone = SymmetricTensor<2,dim>();

  zerozero[0][0] = 1;
  onezero[1][0]  = 0.5; // because symmetry doubles this entry
  oneone[1][1]   = 1;

  double finite_difference_accuracy = 1e-7;
  double finite_difference_factor = 1+finite_difference_accuracy;

  bool Error = false;

  MaterialModelInputs<dim> in_dviscositydpressure(in_base);
  in_dviscositydpressure.pressure[0] *= finite_difference_factor;
  in_dviscositydpressure.pressure[1] *= finite_difference_factor;
  in_dviscositydpressure.pressure[2] *= finite_difference_factor;
  in_dviscositydpressure.pressure[3] *= finite_difference_factor;
  in_dviscositydpressure.pressure[4] *= finite_difference_factor;

  MaterialModelInputs<dim> in_dviscositydstrainrate_zerozero(in_base);
  MaterialModelInputs<dim> in_dviscositydstrainrate_onezero(in_base);
  MaterialModelInputs<dim> in_dviscositydstrainrate_oneone(in_base);

  in_dviscositydstrainrate_zerozero.strain_rate[0] += std::fabs(in_dviscositydstrainrate_zerozero.strain_rate[0][0][0]) * finite_difference_accuracy * zerozero;
  in_dviscositydstrainrate_zerozero.strain_rate[1] += std::fabs(in_dviscositydstrainrate_zerozero.strain_rate[1][0][0]) * finite_difference_accuracy * zerozero;
  in_dviscositydstrainrate_zerozero.strain_rate[2] += std::fabs(in_dviscositydstrainrate_zerozero.strain_rate[2][0][0]) * finite_difference_accuracy * zerozero;
  in_dviscositydstrainrate_zerozero.strain_rate[3] += std::fabs(in_dviscositydstrainrate_zerozero.strain_rate[3][0][0]) * finite_difference_accuracy * zerozero;
  in_dviscositydstrainrate_zerozero.strain_rate[4] += std::fabs(in_dviscositydstrainrate_zerozero.strain_rate[4][0][0]) * finite_difference_accuracy * zerozero;
  in_dviscositydstrainrate_onezero.strain_rate[0]  += std::fabs(in_dviscositydstrainrate_onezero.strain_rate[0][1][0]) * finite_difference_accuracy * onezero;
  in_dviscositydstrainrate_onezero.strain_rate[1]  += std::fabs(in_dviscositydstrainrate_onezero.strain_rate[1][1][0]) * finite_difference_accuracy * onezero;
  in_dviscositydstrainrate_onezero.strain_rate[2]  += std::fabs(in_dviscositydstrainrate_onezero.strain_rate[2][1][0]) * finite_difference_accuracy * onezero;
  in_dviscositydstrainrate_onezero.strain_rate[3]  += std::fabs(in_dviscositydstrainrate_onezero.strain_rate[3][1][0]) * finite_difference_accuracy * onezero;
  in_dviscositydstrainrate_onezero.strain_rate[4]  += std::fabs(in_dviscositydstrainrate_onezero.strain_rate[4][1][0]) * finite_difference_accuracy * onezero;
  in_dviscositydstrainrate_oneone.strain_rate[0]   += std::fabs(in_dviscositydstrainrate_oneone.strain_rate[0][1][1]) * finite_difference_accuracy * oneone;
  in_dviscositydstrainrate_oneone.strain_rate[1]   += std::fabs(in_dviscositydstrainrate_oneone.strain_rate[1][1][1]) * finite_difference_accuracy * oneone;
  in_dviscositydstrainrate_oneone.strain_rate[2]   += std::fabs(in_dviscositydstrainrate_oneone.strain_rate[2][1][1]) * finite_difference_accuracy * oneone;
  in_dviscositydstrainrate_oneone.strain_rate[3]   += std::fabs(in_dviscositydstrainrate_oneone.strain_rate[3][1][1]) * finite_difference_accuracy * oneone;
  in_dviscositydstrainrate_oneone.strain_rate[4]   += std::fabs(in_dviscositydstrainrate_oneone.strain_rate[4][1][1]) * finite_difference_accuracy * oneone;

  MaterialModelInputs<dim> in_dviscositydtemperature(in_base);
  in_dviscositydtemperature.temperature[0] *= 1.0000000001;
  in_dviscositydtemperature.temperature[1] *= 1.0000000001;
  in_dviscositydtemperature.temperature[2] *= 1.0000000001;
  in_dviscositydtemperature.temperature[3] *= 1.0000000001;
  in_dviscositydtemperature.temperature[4] *= 1.0000000001;


  MaterialModelOutputs<dim> out_base(5,3);

  MaterialModelOutputs<dim> out_dviscositydpressure(5,3);
  MaterialModelOutputs<dim> out_dviscositydstrainrate_zerozero(5,3);
  MaterialModelOutputs<dim> out_dviscositydstrainrate_onezero(5,3);
  MaterialModelOutputs<dim> out_dviscositydstrainrate_oneone(5,3);
  MaterialModelOutputs<dim> out_dviscositydtemperature(5,3);

  if (out_base.get_additional_output<MaterialModelDerivatives<dim> >() != NULL)
    throw "error";

  out_base.additional_outputs.push_back(std::make_shared<MaterialModelDerivatives<dim> > (5));

  SimpleNonlinear<dim> mat;
  ParameterHandler prm;
  mat.declare_parameters(prm);

  prm.enter_subsection("Compositional fields");
  {
    prm.set("Number of fields","3");
  }
  prm.leave_subsection();
  prm.enter_subsection("Material model");
  {
    prm.enter_subsection ("Simple nonlinear");
    {
      prm.set ("Viscosity prefactor", "1e-37,1e-36,1e-35,5e-36");
      prm.set ("Viscosity averaging p", std::to_string(parameter));
      prm.set ("Minimum strain rate", 1.4e-20);
    }
    prm.leave_subsection();
  }
  prm.leave_subsection();

  mat.parse_parameters(prm);

  mat.evaluate(in_base, out_base);
  mat.evaluate(in_dviscositydpressure, out_dviscositydpressure);
  mat.evaluate(in_dviscositydstrainrate_zerozero, out_dviscositydstrainrate_zerozero);
  mat.evaluate(in_dviscositydstrainrate_onezero, out_dviscositydstrainrate_onezero);
  mat.evaluate(in_dviscositydstrainrate_oneone, out_dviscositydstrainrate_oneone);
  mat.evaluate(in_dviscositydtemperature, out_dviscositydtemperature);

  //set up additional output for the derivatives
  MaterialModelDerivatives<dim> *derivatives;
  derivatives = out_base.get_additional_output<MaterialModelDerivatives<dim> >();

  double temp;
  for (unsigned int i = 0; i < 5; i++)
    {
      // prevent division by zero. If it is zero, the test has passed, because or
      // the finite difference and the analytical result match perfectly, or (more
      // likely) the material model in independent of this variable.
      temp = (out_dviscositydpressure.viscosities[i] - out_base.viscosities[i]);
      if (in_base.pressure[i] != 0)
        {
          temp /= (in_base.pressure[i] * finite_difference_accuracy);
        }

      if (temp > derivatives->viscosity_derivative_wrt_pressure[i] * finite_difference_factor || temp < derivatives->viscosity_derivative_wrt_pressure[i] * (2-finite_difference_factor))
        {
          std::cout << "Error: The derivative of the viscosity to the pressure is too different from the analytical value." << std::endl;
          Error = true;
        }

    }

  for (unsigned int i = 0; i < 5; i++)
    {
      // prevent division by zero. If it is zero, the test has passed, because or
      // the finite difference and the analytical result match perfectly, or (more
      // likely) the material model in independent of this variable.
      temp = out_dviscositydstrainrate_zerozero.viscosities[i] - out_base.viscosities[i];
      if (temp != 0)
        {
          temp /= std::fabs(in_dviscositydstrainrate_zerozero.strain_rate[i][0][0]) * finite_difference_accuracy;
        }
      std::cout << "zerozero " << i << ": Finite difference = " << temp << ". Analytical derivative = " << derivatives->viscosity_derivative_wrt_strain_rate[i][0][0]  << std::endl;
      if (std::fabs(temp - derivatives->viscosity_derivative_wrt_strain_rate[i][0][0]) > 1e-3 * (std::fabs(temp) + std::fabs(derivatives->viscosity_derivative_wrt_strain_rate[i][0][0])))
        {
          std::cout << "   Error: The derivative of the viscosity to the strain rate is too different from the analytical value." << std::endl;
          Error = true;
        }



    }

  for (unsigned int i = 0; i < 5; i++)
    {
      // prevent division by zero. If it is zero, the test has passed, because or
      // the finite difference and the analytical result match perfectly, or (more
      // likely) the material model in independent of this variable.
      temp = out_dviscositydstrainrate_onezero.viscosities[i] - out_base.viscosities[i];
      if (temp != 0)
        {
          temp /= std::fabs(in_dviscositydstrainrate_onezero.strain_rate[i][1][0]) * finite_difference_accuracy;
        }
      std::cout << "onezero " << i << ": Finite difference = " << temp << ". Analytical derivative = " << derivatives->viscosity_derivative_wrt_strain_rate[i][1][0]   << std::endl;
      if (std::fabs(temp - derivatives->viscosity_derivative_wrt_strain_rate[i][1][0]) > 1e-3 * (std::fabs(temp) + std::fabs(derivatives->viscosity_derivative_wrt_strain_rate[i][1][0])) )
        {
          std::cout << "   Error: The derivative of the viscosity to the strain rate is too different from the analytical value." << std::endl;
          Error = true;
        }
    }

  for (unsigned int i = 0; i < 5; i++)
    {
      // prevent division by zero. If it is zero, the test has passed, because or
      // the finite difference and the analytical result match perfectly, or (more
      // likely) the material model in independent of this variable.
      temp = out_dviscositydstrainrate_oneone.viscosities[i] - out_base.viscosities[i];
      if (temp != 0)
        {
          temp /= std::fabs(in_dviscositydstrainrate_oneone.strain_rate[i][1][1]) * finite_difference_accuracy;
        }
      std::cout << "oneone " << i << ": Finite difference = " << temp << ". Analytical derivative = " << derivatives->viscosity_derivative_wrt_strain_rate[i][1][1]  << std::endl;
      if (std::fabs(temp - derivatives->viscosity_derivative_wrt_strain_rate[i][1][1]) > 1e-3 * (std::fabs(temp) + std::fabs(derivatives->viscosity_derivative_wrt_strain_rate[i][1][1])) )
        {
          std::cout << "   Error: The derivative of the viscosity to the strain rate is too different from the analytical value." << std::endl;
          Error = true;
        }

    }

  if (Error)
    {
      std::cout << "Some parts of the test where not succesfull." << std::endl;
    }
  else
    {
      std::cout << "OK" << std::endl;
    }

  return 42;
}

int exit_function()
{
  exit(0);
  return 42;
}
// run this function by initializing a global variable by it
int ii = f(-1000); // Testing min function
int iz = f(-2); // Testing generalized p norm mean with negative p
int ij = f(-1.5); // Testing generalized p norm mean with negative, non int p
int ik = f(-1); // Testing harmonic mean
int ji = f(0); // Testing geometric mean
int jj = f(1); // Testing arithmetic mean
int jk = f(2); // Testing generalized p norm mean with positive p
int kj = f(1000); // Testing max function
int kl = exit_function();


