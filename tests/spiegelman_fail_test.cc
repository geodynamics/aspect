#include <aspect/simulator.h>
#include <deal.II/grid/tria.h>
#include <aspect/simulator_access.h>
#include <aspect/newton.h>
#include <aspect/parameters.h>

#include <iostream>

#include <aspect/utilities.h>

#include <aspect/material_model/interface.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * The same material model as Drucker Prager, but this one supports multiple
     * compositions in a way to reproduce the Spiegelman 2016 paper.
     *
     * The model is considered incompressible, following the definition
     * described in Interface::is_compressible.
     *
     * The viscosity is computed according to the Drucker Prager frictional
     * plasticity criterion based on a user-defined internal angle of friction $\phi$
     * and cohesion $C$. In 3D:
     * $\sigma_y = \frac{6 C \cos(\phi)}{\sqrt(3) (3+\sin(\phi))} +
     * \frac{6 P \sin(\phi)}{\sqrt(3) (3+\sin(\phi))}$,
     * where $P$ is the pressure.
     * See for example Zienkiewicz, O. C., Humpheson, C. and Lewis, R. W. (1975),
     * G&eacute;otechnique 25, No. 4, 671-689.
     * With this formulation we circumscribe instead of inscribe the Mohr Coulomb
     * yield surface.
     * In 2D the Drucker Prager yield surface is the same
     * as the Mohr Coulomb surface:
     * $\sigma_y = P \sin(\phi) + C \cos(\phi)$.
     * Note that in 2D for $\phi=0$, these criteria
     * revert to the von Mises criterion (no pressure dependence).
     * See for example Thieulot, C. (2011), PEPI 188, 47-68.
     *
     * Note that we enforce the pressure to be positive in the computation of
     * the yield strength by replacing it with
     * a zero value whenever it is negative to prevent negative
     * yield strengths and viscosities.
     * We then use the computed yield strength to scale back the viscosity on
     * to the yield surface using the Viscosity Rescaling Method described in
     * Kachanov, L. M. (2004), Fundamentals of the Theory of Plasticity,
     * Dover Publications, Inc.
     *
     * To avoid numerically unfavourably large (or even negative) viscosity ranges,
     * we regularize the viscosity by $\eta_{\text{eff}} = \eta_{\text{ref}} * \eta
     * / (\eta_{\text{ref}} + \eta)$. Futhermore,
     * we cut off the viscosity with a user-defined minimum and maximum viscosity:
     * $\eta_eff = \frac{1}{\frac{1}{\eta_min + \eta}+
     * \frac{1}{\eta_max}}$.
     *
     * Note that this model uses the formulation that assumes an incompressible
     * medium despite the fact that the density follows the law
     * $\rho(T)=\rho_0(1-\beta(T-T_{\text{ref}}))$.
     *
     *
     * @ingroup MaterialModels
     */

    template <int dim>
    class SpiegelmanMaterial : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        double compute_second_invariant(const SymmetricTensor<2,dim> strain_rate, const double min_strain_rate) const;

        double compute_viscosity(const double edot_ii,const double pressure,const int comp, const double constant_viscosity,const bool compute_full_viscosity, const double min_visc, const double max_visc) const;

        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const;

        /**
         * Return whether the model is compressible or not.  Incompressibility
         * does not necessarily imply that the density is constant; rather, it
         * may still depend on temperature or pressure. In the current
         * context, compressibility means whether we should solve the continuity
         * equation as $\nabla \cdot (\rho \mathbf u)=0$ (compressible Stokes)
         * or as $\nabla \cdot \mathbf{u}=0$ (incompressible Stokes).
         *
         * This material model is incompressible.
         */
        virtual bool is_compressible () const;

        virtual double reference_viscosity () const;

        virtual double reference_density () const;

        static
        void
        declare_parameters (ParameterHandler &prm);

        virtual
        void
        parse_parameters (ParameterHandler &prm);

      private:

        double reference_T;


        /**
         * Whether or not to use the deviator of the strain-rate in the
         * viscosity calculation.
         */
        bool use_deviator_of_strain_rate;

        /**
         * Defining a minimum strain rate stabilizes the viscosity calculation,
         * which involves a division by the strain rate. Units: $\\si{\\per\\second}$.
         */
        std::vector<double> min_strain_rate;
        std::vector<double> min_visc;
        std::vector<double> max_visc;
        std::vector<double> veff_coefficient;
        double ref_visc;
        double reference_compressibility;
        std::vector<double> ref_visc_list;

        std::vector<double> thermal_diffusivity;
        std::vector<double> heat_capacity;


        std::vector<double> densities;
        std::vector<double> thermal_expansivities;

        std::vector<double> cohesion;
        std::vector<double> phi;

        /**
         * The friction angles phi for each of the compositional fields. Since we use their
         * sines and cosines a lot, we only compute these values (for each of the phis) once
         * and cache them.
         */
        std::vector<double> sin_phi;
        std::vector<double> cos_phi;

        std::vector<double> constant_viscosity;

        unsigned int n_fields;

        /**
         * This parameter is used in the averaging of the  of the different compositions.
         * It is implemented as a generalized mean, where p < -1000 is the minimum, p ==-1
         * is the harmonic average, p == 0 is the geometric average, p == 1 is the arithmetic
         * average and p > 1000 is the maximum.
         */
        double viscosity_averaging_p;

        /**
         * This value determines whether to use the analytical or the
         * finite difference derivative for the Newton solver. Although
         * the analytical derivative should be much faster to compute, the
         * finite difference derivative can be useful when changing the
         * rheology.
         */
        bool use_analytical_derivative;
    };

  }


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
    double
    SpiegelmanMaterial<dim>::
    compute_second_invariant(const SymmetricTensor<2,dim> strain_rate, const double min_strain_rate) const
    {
      const double edot_ii_strict = std::sqrt(strain_rate*strain_rate);
      const double edot_ii = std::max(edot_ii_strict, min_strain_rate);
      return edot_ii;
    }

    template <int dim>
    double
    SpiegelmanMaterial<dim>::
    compute_viscosity(const double edot_ii,
                      const double pressure,
                      const int comp,
                      const double constant_viscosity,
                      const bool compute_full_viscosity,
                      const double min_visc,
                      const double max_visc) const
    {
      double viscosity;

      if (constant_viscosity == 0)
        {
          const double strength = ( (dim==3)
                                    ?
                                    ( 6.0 * cohesion[comp] * cos_phi[comp] + 6.0 * std::max(pressure,0.0) * sin_phi[comp] )
                                    / ( std::sqrt(3) * ( 3.0 + sin_phi[comp] ) )
                                    :
                                    cohesion[comp] * cos_phi[comp] + std::max(pressure,0.0) * sin_phi[comp] );

          // Rescale the viscosity back onto the yield surface
          if (strength != 0 && edot_ii != 0)
            viscosity = strength / ( 2.0 * std::sqrt(0.5) * edot_ii );
          else
            viscosity = ref_visc;

          if (compute_full_viscosity == true)
            viscosity = ref_visc * viscosity / (ref_visc + viscosity);
        }
      else
        {
          viscosity = constant_viscosity;
        }

      if (compute_full_viscosity == true)
        return std::max(std::min(viscosity, max_visc), min_visc);
      else
        return viscosity;
    }

    template <int dim>
    void
    SpiegelmanMaterial<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      //set up additional output for the derivatives
      MaterialModelDerivatives<dim> *derivatives;
      derivatives = out.template get_additional_output<MaterialModelDerivatives<dim> >();

      for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
        {
          const double temperature = in.temperature[i];
          const double pressure = in.pressure[i];

          // Averaging composition-field dependent properties
          // Compositions
          //This assert may be help when designing tests for this material model
          AssertThrow(in.composition[i].size()+1 == n_fields,
                      ExcMessage("Number of compositional fields + 1 not equal to number of fields given in input file."));

          const std::vector<double> volume_fractions = MaterialUtilities::compute_field_fractions(in.composition[i]);
          double density = 0.0;
          for (unsigned int c=0; c < volume_fractions.size(); ++c)
            {
              const double temperature_factor = (1 - thermal_expansivities[c] * (temperature - reference_T));
              density += volume_fractions[c] * densities[c] * std::exp(reference_compressibility * pressure) * temperature_factor;
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
              std::vector<SymmetricTensor<2,dim> > composition_viscosities_derivatives(volume_fractions.size());
              std::vector<double> composition_dviscosities_dpressure(volume_fractions.size());

              const SymmetricTensor<2,dim> edot = use_deviator_of_strain_rate ? deviator(in.strain_rate[i]) : in.strain_rate[i];

              for (unsigned int c=0; c < volume_fractions.size(); ++c)
                {
                  // If strain rate is zero (like during the first time step) set it to some very small number
                  // to prevent a division-by-zero, and a floating point exception.
                  // Otherwise, calculate the square-root of the norm of the second invariant of the deviatoric-
                  // strain rate (often simplified as epsilondot_ii)

                  const double edot_ii = compute_second_invariant(edot, min_strain_rate[c]);

                  // Find effective viscosities for each of the individual phases
                  // Viscosities should have same number of entries as compositional fields

                  // Regularized Drucker Prager viscosity (see Spiegelman et al, 2016)
                  composition_viscosities[c] = compute_viscosity(edot_ii,pressure,c,constant_viscosity[c],true,min_visc[c],max_visc[c]);

                  Assert(dealii::numbers::is_finite(composition_viscosities[c]),ExcMessage ("Error: Viscosity is not finite."));

                  if (derivatives != nullptr)
                    {
                      if (use_analytical_derivative)
                        {
                          //analytic
                          if (constant_viscosity[c] == 0  && composition_viscosities[c] <= max_visc[c] && composition_viscosities[c] >= min_visc[c])
                            {
                              // we only want the pure durcker prager here, so the regularization should be off and the min and max values should not be used.
                              // Therefore we set the compute_full_viscosity flag to flase. The min and max viscosity values  won't be used with this flag,
                              // so we enter dummy values.
                              const double drucker_prager_viscosity = compute_viscosity(edot_ii,pressure,c,constant_viscosity[c],false,std::numeric_limits<double>::min(),std::numeric_limits<double>::max());
                              const double regularization_adjustment = (ref_visc * ref_visc)
                                                                       / (ref_visc * ref_visc + 2.0 * ref_visc * drucker_prager_viscosity
                                                                          + drucker_prager_viscosity * drucker_prager_viscosity);

                              composition_viscosities_derivatives[c] = -regularization_adjustment *
                                                                       (drucker_prager_viscosity / (edot_ii * edot_ii)) * edot;

                              if (use_deviator_of_strain_rate == true)
                                composition_viscosities_derivatives[c] = composition_viscosities_derivatives[c] * deviator_tensor<dim>();

                              composition_dviscosities_dpressure[c] = regularization_adjustment *
                                                                      ((dim == 3)
                                                                       ?
                                                                       6 * sin_phi[c] / (std::sqrt(3) * (3.0 + sin_phi[c]) * 2.0 * std::sqrt(0.5) * edot_ii)
                                                                       :
                                                                       sin_phi[c] / (2.0 * std::sqrt(0.5) * edot_ii));
                            }
                          else
                            {
                              composition_viscosities_derivatives[c] = 0;
                              composition_dviscosities_dpressure[c] = 0;
                            }
                        }
                      else
                        {
                          // finite difference
                          const double finite_difference_accuracy = 1e-7;
                          // For each independent component, compute the derivative.
                          for (unsigned int component = 0; component < SymmetricTensor<2,dim>::n_independent_components; ++component)
                            {
                              // compute which of the independent index of the strain-rate tensor we are now looking at.
                              const TableIndices<2> strain_rate_indices = SymmetricTensor<2,dim>::unrolled_to_component_indices (component);

                              // add a small difference to one independent component of hte strain-rate tensor
                              const SymmetricTensor<2,dim> strain_rate_difference_plus = edot + std::max(edot_ii, min_strain_rate[c]) * finite_difference_accuracy
                                                                                         * Utilities::nth_basis_for_symmetric_tensors<dim>(component);
                              const double second_invariant_strain_rate_difference_plus = compute_second_invariant(strain_rate_difference_plus, min_strain_rate[c]);
                              const double eta_component_plus = compute_viscosity(second_invariant_strain_rate_difference_plus,pressure,c,constant_viscosity[c],true,min_visc[c],max_visc[c]);

                              // compute the difference between the viscosity with and without the strain-rate difference.
                              double viscosity_derivative = eta_component_plus - composition_viscosities[c];
                              if (viscosity_derivative != 0)
                                {
                                  // when the difference is non-zero, divide by the difference.
                                  viscosity_derivative /= std::max(edot_ii, min_strain_rate[c]) * finite_difference_accuracy;
                                }

                              composition_viscosities_derivatives[c][strain_rate_indices] = viscosity_derivative;
                            }

                          /**
                           * Now compute the finite difference derivative of the viscoisty to the pressure
                           */
                          double pressure_difference = in.pressure[i] + (std::fabs(in.pressure[i]) * finite_difference_accuracy);

                          double pressure_difference_eta = compute_viscosity(edot_ii, pressure_difference,c,constant_viscosity[c],true,min_visc[c],max_visc[c]);
                          double deriv_pressure = pressure_difference_eta - composition_viscosities[c];


                          if (pressure_difference_eta != 0)
                            {
                              if (in.pressure[i] != 0)
                                {
                                  deriv_pressure /= std::fabs(in.pressure[i]) * finite_difference_accuracy;
                                }
                              else
                                {
                                  deriv_pressure = 0;
                                }
                            }

                          composition_dviscosities_dpressure[c] = deriv_pressure;
                        }
                    }
                }
              out.viscosities[i] = Utilities::weighted_p_norm_average(volume_fractions, composition_viscosities, viscosity_averaging_p);
              Assert(dealii::numbers::is_finite(out.viscosities[i]),ExcMessage ("Error: Averaged viscosity is not finite."));

              if (derivatives != nullptr)
                {
                  derivatives->viscosity_derivative_wrt_strain_rate[i] = Utilities::derivative_of_weighted_p_norm_average(out.viscosities[i],volume_fractions, composition_viscosities, composition_viscosities_derivatives, viscosity_averaging_p);
                  derivatives->viscosity_derivative_wrt_pressure[i] = Utilities::derivative_of_weighted_p_norm_average(out.viscosities[i],volume_fractions, composition_viscosities, composition_dviscosities_dpressure, viscosity_averaging_p);

                  Assert(dealii::numbers::is_finite(derivatives->viscosity_derivative_wrt_pressure[i]),
                         ExcMessage ("Error: Averaged dviscosities_dpressure is not finite."));

                  for (int x = 0; x < dim; ++x)
                    for (int y = 0; y < dim; ++y)
                      Assert(dealii::numbers::is_finite(derivatives->viscosity_derivative_wrt_strain_rate[i][x][y]),
                             ExcMessage ("Error: Averaged dviscosities_dstrain_rate is not finite."));
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
          out.compressibilities[i] = reference_compressibility;
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
    SpiegelmanMaterial<dim>::
    reference_viscosity () const
    {
      return ref_visc;
    }

    template <int dim>
    double
    SpiegelmanMaterial<dim>::
    reference_density () const
    {
      return densities[0];
    }

    template <int dim>
    bool
    SpiegelmanMaterial<dim>::
    is_compressible () const
    {
      return (reference_compressibility != 0);
    }

    template <int dim>
    void
    SpiegelmanMaterial<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Compositional fields");
      {
        prm.declare_entry ("Number of fields", "0",
                           Patterns::Integer (0),
                           "The number of fields that will be advected along with the flow field, excluding "
                           "velocity, pressure and temperature.");

        // temperature and density parameters
        prm.declare_entry ("List of conductivities", "2.25",
                           Patterns::List (Patterns::Double(0)),
                           "A list of thermal conductivities equal to the number of "
                           "compositional fields.");
        prm.declare_entry ("List of capacities", "1250",
                           Patterns::List (Patterns::Double(0)),
                           "A list of heat capacities equal to the number of "
                           "compositional fields.");
        prm.declare_entry ("Thermal expansivities", "3.5e-5",
                           Patterns::List(Patterns::Double(0)),
                           "List of thermal expansivities for background mantle and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "If only one value is given, then all use the same value.  Units: \\si{\\per\\kelvin}");
        prm.declare_entry ("List of reference densities", "2700",
                           Patterns::List (Patterns::Double(0)),
                           "A list of reference densities equal to the number of "
                           "compositional fields.");

        // visocosity parameters
        prm.declare_entry ("List of cohesions", "1e20",
                           Patterns::List (Patterns::Double(0)),
                           "A list of initial viscosities equal to the number of "
                           "compositional fields.");
        prm.declare_entry ("List of angles of internal friction", "30",
                           Patterns::List (Patterns::Double(0)),
                           "A list of initial viscosities equal to the number of "
                           "compositional fields.");
        prm.declare_entry ("List of initial viscosities", "1e21",
                           Patterns::List (Patterns::Double(0)),
                           "A list of initial viscosities equal to the number of "
                           "compositional fields.");
        prm.declare_entry ("List of constant viscosities", "0",
                           Patterns::List (Patterns::Double(0)),
                           "A list of viscous viscosities equal to the number of "
                           "compositional fields. When it is zero, the Drucker Prager "
                           "rheology is used, when it is nonzero, it sets the viscosity "
                           "with the set value.");

      }
      prm.leave_subsection();

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection ("Spiegelman 2016");
        {
          // Reference and minimum/maximum values
          prm.declare_entry ("Reference temperature", "293", Patterns::Double(0),
                             "For calculating density by thermal expansivity. Units: \\si{\\kelvin}");
          prm.declare_entry ("Minimum strain rate", "1.96e-40", Patterns::List(Patterns::Double(0)),
                             "Stabilizes strain dependent viscosity. Units: \\si{\\per\\second}");
          prm.declare_entry ("Minimum viscosity", "1e10", Patterns::List(Patterns::Double(0)),
                             "Lower cutoff for effective viscosity. Units: \\si{\\pascal\\second}");
          prm.declare_entry ("Maximum viscosity", "1e28", Patterns::List(Patterns::Double(0)),
                             "Upper cutoff for effective viscosity. Units: \\si{\\pascal\\second}");
          prm.declare_entry ("Effective viscosity coefficient", "1.0", Patterns::List(Patterns::Double(0)),
                             "Scaling coefficient for effective viscosity.");
          prm.declare_entry ("Reference viscosity", "1e22", Patterns::List(Patterns::Double(0)),
                             "Reference viscosity for nondimensionalization. Units $Pa s$");
          prm.declare_entry ("Reference compressibility", "4e-12", Patterns::Double (0),
                             "The value of the reference compressibility. Units: \\si{\\per\\pascal}.");

          // averaging parameters
          prm.declare_entry ("Viscosity averaging p", "-1",
                             Patterns::Double(),
                             "This is the p value in the generalized weighed average eqation: "
                             " mean = \\frac{1}{k}(\\sum_{i=1}^k \\big(c_i \\eta_{\\text{eff}_i}^p)\\big)^{\\frac{1}{p}}. "
                             " Units: \\si{\\pascal\\second}");

          // finite difference versus analytical
          prm.declare_entry ("Use analytical derivative", "true",
                             Patterns::Bool(),
                             "This value determines whether to use the analytical or the finite difference derivative "
                             "for the Newton solver. Although the analytical derivative should be much faster to compute, "
                             "the finite difference derivative can be useful when changing the rheology.");

          prm.declare_entry ("Use deviator of strain-rate", "true",
                             Patterns::Bool(),
                             "Whether or not to use the deviator of the strain-rate in the viscosity calculation.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    SpiegelmanMaterial<dim>::parse_parameters (ParameterHandler &prm)
    {
      using namespace Utilities;
      // can't use this->n_compositional_fields(), because some
      // tests never initiate the simulator, but uses the material
      // model directly.
      prm.enter_subsection ("Compositional fields");
      {
        n_fields = prm.get_integer ("Number of fields")+1;

        // temperature and density parameters
        thermal_diffusivity = possibly_extend_from_1_to_N(string_to_double(split_string_list(prm.get("List of conductivities"))),
                                                          n_fields,
                                                          "List of conductivities");
        heat_capacity = possibly_extend_from_1_to_N(string_to_double(split_string_list(prm.get("List of capacities"))),
                                                    n_fields,
                                                    "List of capacities");
        densities = possibly_extend_from_1_to_N(string_to_double(split_string_list(prm.get("List of reference densities"))),
                                                n_fields,
                                                "List of refdens");
        thermal_expansivities = possibly_extend_from_1_to_N(string_to_double(split_string_list(prm.get("Thermal expansivities"))),
                                                            n_fields,
                                                            "Thermal expansivities");

        // Rheological parameters
        cohesion = possibly_extend_from_1_to_N(string_to_double(split_string_list(prm.get("List of cohesions"))),
                                               n_fields,
                                               "List of cohesions");
        phi = possibly_extend_from_1_to_N(string_to_double(split_string_list(prm.get("List of angles of internal friction"))),
                                          n_fields,
                                          "List of angles of internal friction");


        // Todo: make use of the initial viscosity
        ref_visc_list = possibly_extend_from_1_to_N(string_to_double(split_string_list(prm.get("List of initial viscosities"))),
                                                    n_fields,
                                                    "List of initial viscosities");

        sin_phi.resize(n_fields);
        cos_phi.resize(n_fields);
        for (unsigned int c = 0; c < n_fields; ++c)
          {
            sin_phi[c] = std::sin(phi[c] * numbers::PI/180);
            cos_phi[c] = std::cos(phi[c] * numbers::PI/180);
          }

        constant_viscosity = possibly_extend_from_1_to_N(string_to_double(split_string_list(prm.get("List of constant viscosities"))),
                                                         n_fields,
                                                         "List of constant viscosities");
      }
      prm.leave_subsection();

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection ("Spiegelman 2016");
        {

          // Reference and minimum/maximum values
          reference_T = prm.get_double("Reference temperature");
          ref_visc = prm.get_double ("Reference viscosity");
          min_strain_rate = possibly_extend_from_1_to_N(string_to_double(split_string_list(prm.get("Minimum strain rate"))),
                                                        n_fields,
                                                        "Minimum strain rate");
          min_visc = get_vector_double ("Minimum viscosity",n_fields,prm);
          max_visc = get_vector_double ("Maximum viscosity",n_fields,prm);
          veff_coefficient = get_vector_double ("Effective viscosity coefficient",n_fields,prm);
          reference_compressibility  = prm.get_double ("Reference compressibility");
          use_deviator_of_strain_rate = prm.get_bool ("Use deviator of strain-rate");


          // averaging parameters
          viscosity_averaging_p = prm.get_double("Viscosity averaging p");

          use_analytical_derivative = prm.get_bool("Use analytical derivative");

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
    ASPECT_REGISTER_MATERIAL_MODEL(SpiegelmanMaterial,
                                   "spiegelman 2016",
                                   "An implementation of the spiegelman 2016 benchmark paper in gcubed "
                                   "(doi:10.1002/ 2015GC006228). It implements a regularized Drucker Prager "
                                   "material (when the parameter constant viscosity is 0) and a linear material "
                                   "(when the parameter constant viscosity is not 0), which are both needed to "
                                   "implment the Spiegelman benchmark. The derivatives for the Newton solver "
                                   "have been implemented for this material model.")
  }
}
