/*
  Copyright (C) 2015 - 2024 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_melt_visco_plastic_h
#define _aspect_material_model_melt_visco_plastic_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator.h>
#include <aspect/simulator_access.h>
#include <aspect/postprocess/melt_statistics.h>
#include <aspect/melt.h>
#include <aspect/utilities.h>

#include <deal.II/numerics/fe_field_function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/parsed_function.h>

namespace aspect
{
  namespace MaterialModel
  {
    /**
     * Additional output fields for the plastic parameters weakened (or hardened)
     * by strain to be added to the MaterialModel::MaterialModelOutputs structure
     * and filled in the MaterialModel::Interface::evaluate() function.
     */
    template <int dim>
    class PlasticAdditionalOutputs2 : public NamedAdditionalMaterialOutputs<dim>
    {
      public:
        PlasticAdditionalOutputs2(const unsigned int n_points);

        virtual std::vector<double> get_nth_output(const unsigned int idx) const;

        /**
         * Cohesions at the evaluation points passed to
         * the instance of MaterialModel::Interface::evaluate() that fills
         * the current object.
         */
        std::vector<double> cohesions;

        /**
         * Internal angles of friction at the evaluation points passed to
         * the instance of MaterialModel::Interface::evaluate() that fills
         * the current object.
         */
        std::vector<double> friction_angles;

        /**
         * The area where the viscous stress exceeds the plastic yield strength,
         * and viscosity is rescaled back to the yield envelope.
         */
        std::vector<double> yielding;
    };

    /**
     * A material model that implements a simple formulation of the
     * material parameters required for the modeling of melt transport,
     * including a source term for the porosity according to the melting
     * model for dry peridotite of Katz, 2003. This also includes a
     * computation of the latent heat of melting (if the latent heat
     * heating model is active).
     *
     * Most of the material properties are constant, except for the shear,
     * compaction and melt viscosities and the permeability, which depend on
     * the porosity; and the solid and melt densities, which depend on
     * temperature and pressure.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class MeltViscoPlastic : public MaterialModel::MeltFractionModel<dim>,
      public MaterialModel::MeltInterface<dim>,
      public ::aspect::SimulatorAccess<dim>
    {
      public:

        /**
         * Return whether the model is compressible or not.
         * In this case, Keller et al, 2003 assume the fluid and
         * solid phases are each incompressible: "compressibility
         * in the model is accounted for by changes in melt fraction."
         * Thus, this function will always return false.
         */
        bool is_compressible () const override;

        /**
         * @name Reference quantities
         * @{
         */
        double reference_darcy_coefficient () const override;

        void evaluate(const typename Interface<dim>::MaterialModelInputs &in,
                      typename Interface<dim>::MaterialModelOutputs &out) const override;

        virtual void melt_fractions (const MaterialModel::MaterialModelInputs<dim> &in,
                                     std::vector<double> &melt_fractions) const;

        /**
         * @}
         */

        /**
         * @name Functions used in dealing with run-time parameters
         * @{
         */
        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);
        /**
         * @}
         */

        virtual
        void
        create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const override;

      private:

        double ref_temperature;

        std::vector<double> densities;
        std::vector<double> thermal_expansivities;
        std::vector<double> specific_heats;
        std::vector<double> thermal_conductivities;

        double min_strain_rate;
        double ref_strain_rate;
        double min_viscosity;
        double max_viscosity;

        double ref_viscosity;

        /**
         * Enumeration for selecting which viscosity averaging scheme to use.
         */
        MaterialUtilities::CompositionalAveragingOperation viscosity_averaging;

        std::vector<double> linear_viscosities;

        std::vector<double> angles_internal_friction;
        std::vector<double> cohesions;

        double maximum_yield_stress;

        std::vector<double> strength_reductions;

        double melt_density_change;
        double xi_0;
        double eta_f;
        double reference_permeability;
        double alpha_phi;
        double freezing_rate;
        double melting_time_scale;

        /**
         * Parameters for anhydrous melting of peridotite after Katz, 2003
         */

        // for the solidus temperature
        double A1;   // °C
        double A2; // °C/Pa
        double A3; // °C/(Pa^2)

        // for the lherzolite liquidus temperature
        double B1;   // °C
        double B2;   // °C/Pa
        double B3; // °C/(Pa^2)

        // for the liquidus temperature
        double C1;   // °C
        double C2;  // °C/Pa
        double C3; // °C/(Pa^2)

        // for the reaction coefficient of pyroxene
        double r1;     // cpx/melt
        double r2;     // cpx/melt/GPa
        double M_cpx;  // mass fraction of pyroxene

        // melt fraction exponent
        double beta;

        /**
         * Percentage of material that is molten for a given temperature and
         * pressure (assuming equilibrium conditions). Melting model after Katz,
         * 2003, for dry peridotite.
         */
        virtual
        double
        melt_fraction (const double temperature,
                       const double pressure) const;

        /**
         * A function that fills the plastic additional output in the
         * MaterialModelOutputs object that is handed over, if it exists.
         * Does nothing otherwise.
         */
        void fill_plastic_outputs (const unsigned int point_index,
                                   const std::vector<double> &volume_fractions,
                                   const bool plastic_yielding,
                                   const MaterialModel::MaterialModelInputs<dim> &in,
                                   MaterialModel::MaterialModelOutputs<dim> &out) const;

        /**
         * Whether to use a function defined as an input parameter for the melting rate
         */
        bool use_melting_rate_function;

        /**
         * A function object representing the melting rate.
         */
        Functions::ParsedFunction<dim> function;

    };

  }
}

#endif

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
    PlasticAdditionalOutputs2<dim>::PlasticAdditionalOutputs2 (const unsigned int n_points)
      :
      NamedAdditionalMaterialOutputs<dim>(make_plastic_additional_outputs_names()),
      cohesions(n_points, numbers::signaling_nan<double>()),
      friction_angles(n_points, numbers::signaling_nan<double>()),
      yielding(n_points, numbers::signaling_nan<double>())
    {}

    template <int dim>
    std::vector<double>
    PlasticAdditionalOutputs2<dim>::get_nth_output(const unsigned int idx) const
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
    double
    MeltViscoPlastic<dim>::
    reference_darcy_coefficient () const
    {
      // 0.01 = 1% melt
      return reference_permeability * Utilities::fixed_power<3>(0.01) / eta_f;
    }

    template <int dim>
    bool
    MeltViscoPlastic<dim>::
    is_compressible () const
    {
      return false;
    }

    template <int dim>
    double
    MeltViscoPlastic<dim>::
    melt_fraction (const double temperature,
                   const double pressure) const
    {
      // anhydrous melting of peridotite after Katz, 2003
      const double T_solidus  = A1 + 273.15
                                + A2 * pressure
                                + A3 * pressure * pressure;
      const double T_lherz_liquidus = B1 + 273.15
                                      + B2 * pressure
                                      + B3 * pressure * pressure;
      const double T_liquidus = C1 + 273.15
                                + C2 * pressure
                                + C3 * pressure * pressure;

      // melt fraction for peridotite with clinopyroxene
      double peridotite_melt_fraction;
      if (temperature < T_solidus || pressure > 1.3e10)
        peridotite_melt_fraction = 0.0;
      else if (temperature > T_lherz_liquidus)
        peridotite_melt_fraction = 1.0;
      else
        peridotite_melt_fraction = std::pow((temperature - T_solidus) / (T_lherz_liquidus - T_solidus),beta);

      // melt fraction after melting of all clinopyroxene
      const double R_cpx = r1 + r2 * std::max(0.0, pressure);
      const double F_max = M_cpx / R_cpx;

      if (peridotite_melt_fraction > F_max && temperature < T_liquidus)
        {
          const double T_max = std::pow(F_max,1/beta) * (T_lherz_liquidus - T_solidus) + T_solidus;
          peridotite_melt_fraction = F_max + (1 - F_max) * std::pow((temperature - T_max) / (T_liquidus - T_max),beta);
        }
      return peridotite_melt_fraction;
    }

    template <int dim>
    void
    MeltViscoPlastic<dim>::
    melt_fractions (const MaterialModel::MaterialModelInputs<dim> &in,
                    std::vector<double> &melt_fractions) const
    {
      for (unsigned int q=0; q<in.n_evaluation_points(); ++q)
        melt_fractions[q] = melt_fraction(in.temperature[q],
                                          std::max(0.0, in.pressure[q]));
      return;
    }

    template <int dim>
    void
    MeltViscoPlastic<dim>::
    evaluate(const typename Interface<dim>::MaterialModelInputs &in, typename Interface<dim>::MaterialModelOutputs &out) const
    {
      // 1) Initial viscosities and other material properties
      for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
        {
          const std::vector<double> volume_fractions = MaterialUtilities::compute_composition_fractions(in.composition[i]);
          out.viscosities[i] = MaterialUtilities::average_value(volume_fractions, linear_viscosities, viscosity_averaging);

          out.densities[i] = MaterialUtilities::average_value(volume_fractions, densities, MaterialUtilities::CompositionalAveragingOperation::arithmetic);
          out.thermal_expansion_coefficients[i] = MaterialUtilities::average_value(volume_fractions, thermal_expansivities, MaterialUtilities::CompositionalAveragingOperation::arithmetic);
          out.thermal_conductivities[i] = MaterialUtilities::average_value(volume_fractions, thermal_conductivities, MaterialUtilities::CompositionalAveragingOperation::arithmetic);
          out.specific_heat[i] = MaterialUtilities::average_value(volume_fractions, specific_heats, MaterialUtilities::CompositionalAveragingOperation::arithmetic);

          // TODO: compute the actual number
          out.entropy_derivative_pressure[i]    = 0.0;
          out.entropy_derivative_temperature[i] = 0.0;
        }

      // Store the intrinsic viscosities for computing the compaction viscosities later on
      // (Keller et al. eq. 51).
      const std::vector<double> intrinsic_viscosities = out.viscosities;

      // 2) Retrieve fluid pressure and volumetric strain rate
      std::vector<double> fluid_pressures(in.n_evaluation_points());
      std::vector<double> volumetric_strain_rates(in.n_evaluation_points());
      std::vector<double> volumetric_yield_strength(in.n_evaluation_points());

      ReactionRateOutputs<dim> *reaction_rate_out = out.template get_additional_output<ReactionRateOutputs<dim>>();

      if (this->include_melt_transport() )
        {
          if (in.current_cell.state() == IteratorState::valid)
            {
              std::vector<Point<dim>> quadrature_positions(in.n_evaluation_points());

              for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
                quadrature_positions[i] = this->get_mapping().transform_real_to_unit_cell(in.current_cell, in.position[i]);

              // FEValues requires a quadrature and we provide the default quadrature
              // as we only need to evaluate the solution and gradients.
              FEValues<dim> fe_values (this->get_mapping(),
                                       this->get_fe(),
                                       Quadrature<dim>(quadrature_positions),
                                       update_values | update_gradients);

              fe_values.reinit (in.current_cell);

              // get fluid pressure from the current solution
              const FEValuesExtractors::Scalar extractor_pressure = this->introspection().variable("fluid pressure").extractor_scalar();
              fe_values[extractor_pressure].get_function_values (this->get_solution(),
                                                                 fluid_pressures);

              // get volumetric strain rate
              // see Keller et al. eq. 11.
              fe_values[this->introspection().extractors.velocities].get_function_divergences (this->get_solution(),
                  volumetric_strain_rates);
            }

          const double time_scale = this->convert_output_to_years() ? year_in_seconds : 1.0;

          // 3) Get porosity, melt density and update melt reaction terms
          for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
            {
              // get peridotite and porosity field indices
              const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
              const unsigned int peridotite_idx = this->introspection().compositional_index_for_name("peridotite");

              const double old_porosity = in.composition[i][porosity_idx];
              const double maximum_melt_fraction = in.composition[i][peridotite_idx];

              // calculate the melting rate as difference between the equilibrium melt fraction
              // and the solution of the previous time step
              double porosity_change = 0.0;

              // batch melting
              porosity_change = melt_fraction(in.temperature[i], this->get_adiabatic_conditions().pressure(in.position[i]))
                                - std::max(maximum_melt_fraction, 0.0);
              porosity_change = std::max(porosity_change, 0.0);

              // do not allow negative porosity
              porosity_change = std::max(porosity_change, -old_porosity);

              // because depletion is a volume-based, and not a mass-based property that is advected,
              // additional scaling factors on the right hand side apply
              for (unsigned int c=0; c<in.composition[i].size(); ++c)
                {
                  if (use_melting_rate_function == false)
                    {
                      // fill reaction rate outputs
                      if (reaction_rate_out != nullptr)
                        {
                          if (c == peridotite_idx && this->get_timestep_number() > 0)
                            reaction_rate_out->reaction_rates[i][c] = porosity_change / melting_time_scale
                                                                      * (1 - maximum_melt_fraction) / (1 - old_porosity);
                          else if (c == porosity_idx && this->get_timestep_number() > 0)
                            reaction_rate_out->reaction_rates[i][c] = porosity_change / melting_time_scale;
                          else
                            reaction_rate_out->reaction_rates[i][c] = 0.0;
                        }

                      out.reaction_terms[i][c] = 0.0;
                    }
                  else
                    {
                      // Set reaction rates to 0
                      if (reaction_rate_out != nullptr)
                        reaction_rate_out->reaction_rates[i][c] = 0.;
                      //
                      if (c == porosity_idx)
                        out.reaction_terms[i][c] = function.value(in.position[i]) / time_scale * out.densities[i];
                      else if (c == peridotite_idx)
                        out.reaction_terms[i][c] = 0.0;
                    }
                }

              // 4) Reduce shear viscosity due to melt presence
              const double porosity = std::min(1.0, std::max(in.composition[i][porosity_idx],0.0));
              out.viscosities[i] *= std::exp(- alpha_phi * porosity);
            }
        }

      if (in.requests_property(MaterialProperties::viscosity) )
        {
          // 5) Compute plastic weakening of the viscosity
          for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
            {
              // Compute volume fractions
              const std::vector<double> volume_fractions = MaterialUtilities::compute_composition_fractions(in.composition[i]);

              // 4) Compute plastic weakening of visco(elastic) viscosity
              double porosity = 0.0;

              if (this->include_melt_transport() )
                porosity = std::min(1.0, std::max(in.composition[i][this->introspection().compositional_index_for_name("porosity")],0.0));

              // calculate deviatoric strain rate (Keller et al. eq. 13)
              const double edot_ii = ( (this->get_timestep_number() == 0 && in.strain_rate[i].norm() <= std::numeric_limits<double>::min())
                                       ?
                                       ref_strain_rate
                                       :
                                       std::max(std::sqrt(std::fabs(second_invariant(deviator(in.strain_rate[i])))),
                                                min_strain_rate) );

              // compute viscous stress
              const double viscous_stress = 2. * out.viscosities[i] * edot_ii * (1.0 - porosity);

              // In case porosity lies above the melt transport threshold
              // P_effective = P_bulk - P_f = (1-porosity) * P_s + porosity * P_f - P_f = (1-porosity) * (P_s - P_f)
              // otherwise,
              // P_effective = P_bulk, which equals P_solid (which is given by in.pressure[i])
              const double effective_pressure = ((this->include_melt_transport() && this->get_melt_handler().is_melt_cell(in.current_cell))
                                                 ?
                                                 (1. - porosity) * (in.pressure[i] - fluid_pressures[i])
                                                 :
                                                 in.pressure[i]);

              double yield_strength = 0.0;
              double tensile_strength = 0.0;

              for (unsigned int c=0; c< volume_fractions.size(); ++c)
                {
                  const double tensile_strength_c = cohesions[c]/strength_reductions[c];

                  // Convert friction angle from degrees to radians
                  double phi = angles_internal_friction[c] * constants::degree_to_radians;
                  const double transition_pressure = (cohesions[c] * std::cos(phi) - tensile_strength_c) / (1.0 -  std::sin(phi));

                  double yield_strength_c = 0.0;
                  // In case we're not using the Keller et al. formulation,
                  // or the effective pressure is bigger than the transition pressure, use
                  // the normal yield strength formulation
                  if (effective_pressure > transition_pressure || !this->include_melt_transport())
                    yield_strength_c = ( (dim==3)
                                         ?
                                         ( 6.0 * cohesions[c] * std::cos(phi) + 2.0 * effective_pressure * std::sin(phi) )
                                         / ( std::sqrt(3.0) * (3.0 + std::sin(phi) ) )
                                         :
                                         cohesions[c] * std::cos(phi) + effective_pressure * std::sin(phi) );
                  else
                    // Note typo in Keller et al. paper eq. (37) (minus sign)
                    yield_strength_c = tensile_strength_c + effective_pressure;

                  // TODO add different averagings?
                  yield_strength += volume_fractions[c]*yield_strength_c;
                  tensile_strength += volume_fractions[c]*tensile_strength_c;
                }

              bool plastic_yielding = false;

              // If the viscous stress is greater than the yield strength, rescale the viscosity back to yield surface
              // and reaction term for plastic finite strain
              if (viscous_stress >= yield_strength)
                {
                  plastic_yielding = true;
                  out.viscosities[i] = yield_strength / (2.0 * edot_ii);
                }

              // Calculate average internal friction angle and cohesion values
              const double cohesion = MaterialUtilities::average_value(volume_fractions, cohesions, viscosity_averaging);
              const double angle_internal_friction = MaterialUtilities::average_value(volume_fractions, angles_internal_friction, viscosity_averaging);

              PlasticAdditionalOutputs2<dim> *plastic_out = out.template get_additional_output<PlasticAdditionalOutputs2<dim>>();
              if (plastic_out != nullptr)
                {
                  plastic_out->cohesions[i] = cohesion;
                  plastic_out->friction_angles[i] = angle_internal_friction;
                  plastic_out->yielding[i] = plastic_yielding ? 1 : 0;
                }

              // Limit the viscosity with specified minimum and maximum bounds
              out.viscosities[i] = std::min(std::max(out.viscosities[i], min_viscosity), max_viscosity);

              // Compute the volumetric yield strength (Keller et al. eq (38))
              volumetric_yield_strength[i] = viscous_stress - tensile_strength;
            }
        }

      // fill melt outputs if they exist
      MeltOutputs<dim> *melt_out = out.template get_additional_output<MeltOutputs<dim>>();

      if (melt_out != nullptr)
        {
          for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
            {
              const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
              double porosity = std::min(1.0, std::max(in.composition[i][porosity_idx],0.0));
              melt_out->fluid_viscosities[i] = eta_f;
              melt_out->permeabilities[i] = reference_permeability * Utilities::fixed_power<3>(porosity) * Utilities::fixed_power<2>(1.0-porosity);

              melt_out->fluid_densities[i] = out.densities[i] + melt_density_change;
              melt_out->fluid_density_gradients[i] = 0.0;

              const double compaction_pressure = (1.0 - porosity) * (in.pressure[i] - fluid_pressures[i]);

              const double phi_0 = 0.05;
              porosity = std::max(std::min(porosity,0.995),1.e-8);
              // compaction viscosities (Keller et al. eq (51)
              melt_out->compaction_viscosities[i] = intrinsic_viscosities[i] * phi_0 / porosity;

              // visco(elastic) compaction viscosity
              // Keller et al. eq (36)
              // TODO include elastic part
              melt_out->compaction_viscosities[i] *= (1. - porosity);

              // TODO compaction stress evolution parameter
              // Keller et al. eq. (41) and (44)

              // effective compaction viscosity (Keller et al. eq (43) )
              // NB: I've added a minus sign as according to eq 43
              if (in.strain_rate.size() && compaction_pressure < volumetric_yield_strength[i])
                {
                  // the volumetric strain rate might be negative, but will always have the same sign as the volumetric yield strength
                  if (volumetric_strain_rates[i] >= 0)
                    volumetric_strain_rates[i] = std::max(volumetric_strain_rates[i], min_strain_rate);
                  else
                    volumetric_strain_rates[i] = std::min(volumetric_strain_rates[i], -min_strain_rate);
                  melt_out->compaction_viscosities[i] = - volumetric_yield_strength[i] / volumetric_strain_rates[i];
                }

              // Limit the viscosity with specified minimum and maximum bounds
              melt_out->compaction_viscosities[i] = std::min(std::max(melt_out->compaction_viscosities[i], min_viscosity), max_viscosity);
            }
        }


    }


    template <int dim>
    void
    MeltViscoPlastic<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Melt visco plastic");
        {

          prm.declare_entry ("Reference temperature", "293", Patterns::Double(0),
                             "For calculating density by thermal expansivity. Units: $K$");
          prm.declare_entry ("Densities", "3300.",
                             Patterns::List(Patterns::Double(0)),
                             "List of densities for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value.  Units: $kg / m^3$");
          prm.declare_entry ("Thermal expansivities", "4.e-5",
                             Patterns::List(Patterns::Double(0)),
                             "List of thermal expansivities for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value. Units: $1/K$");
          prm.declare_entry ("Specific heats", "1250.",
                             Patterns::List(Patterns::Double(0)),
                             "List of specific heats $C_p$ for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value. Units: $J /kg /K$");
          prm.declare_entry ("Thermal conductivities", "4.7",
                             Patterns::List(Patterns::Double(0)),
                             "List of thermal conductivities for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value. Units: $W/m/K$ ");

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

          prm.declare_entry ("Viscosity averaging scheme", "harmonic",
                             Patterns::Selection("arithmetic|harmonic|geometric|maximum composition "),
                             "When more than one compositional field is present at a point "
                             "with different viscosities, we need to come up with an average "
                             "viscosity at that point.  Select a weighted harmonic, arithmetic, "
                             "geometric, or maximum composition.");

          prm.declare_entry ("Linear viscosities", "1.e22",
                             Patterns::List(Patterns::Double(0)),
                             "List of linear viscosities for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "The values can be used instead of the viscosities derived from the "
                             "base material model. Units: Pa s.");

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

          prm.declare_entry ("Maximum yield stress", "1e12", Patterns::Double(0),
                             "Limits the maximum value of the yield stress determined by the "
                             "drucker-prager plasticity parameters. Default value is chosen so this "
                             "is not automatically used. Values of 100e6--1000e6 $Pa$ have been used "
                             "in previous models. Units: $Pa$");

          prm.declare_entry ("Host rock strength reductions", "4",
                             Patterns::List(Patterns::Double(0)),
                             "List of reduction factors of strength of the host rock under tensile stress, $R$, "
                             "for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "Units: none.");

          prm.declare_entry ("Melt density change", "-500",
                             Patterns::Double (),
                             "Difference between solid density $\\rho_{s}$ and melt/fluid$\\rho_{f}$. "
                             "Units: $kg/m^3$.");
          prm.declare_entry ("Reference bulk viscosity", "1e22",
                             Patterns::Double (0),
                             "The value of the constant bulk viscosity $\\xi_0$ of the solid matrix. "
                             "This viscosity may be modified by both temperature and porosity "
                             "dependencies. Units: $Pa s$.");
          prm.declare_entry ("Reference melt viscosity", "10",
                             Patterns::Double (0),
                             "The value of the constant melt viscosity $\\eta_f$. Units: $Pa s$.");
          prm.declare_entry ("Exponential melt weakening factor", "27",
                             Patterns::Double (0),
                             "The porosity dependence of the viscosity. Units: dimensionless.");
          prm.declare_entry ("Reference permeability", "1e-8",
                             Patterns::Double(),
                             "Reference permeability of the solid host rock."
                             "Units: $m^2$.");
          prm.declare_entry ("Melting time scale for operator splitting", "1e3",
                             Patterns::Double (0),
                             "Because the operator splitting scheme is used, the porosity field can not "
                             "be set to a new equilibrium melt fraction instantly, but the model has to "
                             "provide a melting time scale instead. This time scale defines how fast melting "
                             "happens, or more specifically, the parameter defines the time after which "
                             "the deviation of the porosity from the equilibrium melt fraction will be "
                             "reduced to a fraction of $1/e$. So if the melting time scale is small compared "
                             "to the time step size, the reaction will be so fast that the porosity is very "
                             "close to the equilibrium melt fraction after reactions are computed. Conversely, "
                             "if the melting time scale is large compared to the time step size, almost no "
                             "melting and freezing will occur."
                             "\n\n"
                             "Also note that the melting time scale has to be larger than or equal to the reaction "
                             "time step used in the operator splitting scheme, otherwise reactions can not be "
                             "computed. "
                             "Units: yr or s, depending on the ``Use years in output instead of seconds'' parameter.");

          prm.declare_entry ("A1", "1085.7",
                             Patterns::Double (),
                             "Constant parameter in the quadratic "
                             "function that approximates the solidus "
                             "of peridotite. "
                             "Units: ${}^\\circ C$.");
          prm.declare_entry ("A2", "1.329e-7",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the quadratic function that approximates "
                             "the solidus of peridotite. "
                             "Units: ${}^\\circ C/Pa$.");
          prm.declare_entry ("A3", "-5.1e-18",
                             Patterns::Double (),
                             "Prefactor of the quadratic pressure term "
                             "in the quadratic function that approximates "
                             "the solidus of peridotite. "
                             "Units: ${}^\\circ C/(Pa^2)$.");
          prm.declare_entry ("B1", "1475.0",
                             Patterns::Double (),
                             "Constant parameter in the quadratic "
                             "function that approximates the lherzolite "
                             "liquidus used for calculating the fraction "
                             "of peridotite-derived melt. "
                             "Units: ${}^\\circ C$.");
          prm.declare_entry ("B2", "8.0e-8",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the quadratic function that approximates "
                             "the  lherzolite liquidus used for "
                             "calculating the fraction of peridotite-"
                             "derived melt. "
                             "Units: ${}^\\circ C/Pa$.");
          prm.declare_entry ("B3", "-3.2e-18",
                             Patterns::Double (),
                             "Prefactor of the quadratic pressure term "
                             "in the quadratic function that approximates "
                             "the  lherzolite liquidus used for "
                             "calculating the fraction of peridotite-"
                             "derived melt. "
                             "Units: ${}^\\circ C/(Pa^2)$.");
          prm.declare_entry ("C1", "1780.0",
                             Patterns::Double (),
                             "Constant parameter in the quadratic "
                             "function that approximates the liquidus "
                             "of peridotite. "
                             "Units: ${}^\\circ C$.");
          prm.declare_entry ("C2", "4.50e-8",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the quadratic function that approximates "
                             "the liquidus of peridotite. "
                             "Units: ${}^\\circ C/Pa$.");
          prm.declare_entry ("C3", "-2.0e-18",
                             Patterns::Double (),
                             "Prefactor of the quadratic pressure term "
                             "in the quadratic function that approximates "
                             "the liquidus of peridotite. "
                             "Units: ${}^\\circ C/(Pa^2)$.");
          prm.declare_entry ("r1", "0.5",
                             Patterns::Double (),
                             "Constant in the linear function that "
                             "approximates the clinopyroxene reaction "
                             "coefficient. "
                             "Units: non-dimensional.");
          prm.declare_entry ("r2", "8e-11",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the linear function that approximates "
                             "the clinopyroxene reaction coefficient. "
                             "Units: $1/Pa$.");
          prm.declare_entry ("beta", "1.5",
                             Patterns::Double (),
                             "Exponent of the melting temperature in "
                             "the melt fraction calculation. "
                             "Units: non-dimensional.");
          prm.declare_entry ("Mass fraction cpx", "0.15",
                             Patterns::Double (),
                             "Mass fraction of clinopyroxene in the "
                             "peridotite to be molten. "
                             "Units: non-dimensional.");

          prm.declare_entry ("Use melting rate function", "false",
                             Patterns::Bool (),
                             "Prescribe the melting rate through a function defined as an input parameter. "
                             "Units: None");

          prm.enter_subsection("Melting rate function");
          {
            Functions::ParsedFunction<dim>::declare_parameters (prm, 1);
          }
          prm.leave_subsection();

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    MeltViscoPlastic<dim>::parse_parameters (ParameterHandler &prm)
    {
      //increment by one for background:
      const unsigned int n_fields = this->n_compositional_fields() + 1;

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Melt visco plastic");
        {

          ref_temperature = prm.get_double ("Reference temperature");
          densities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Densities"))),
                                                              n_fields,
                                                              "Densities");
          thermal_expansivities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Thermal expansivities"))),
                                                                          n_fields,
                                                                          "Thermal expansivities");
          specific_heats = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Specific heats"))),
                                                                   n_fields,
                                                                   "Specific heats");
          thermal_conductivities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Thermal conductivities"))),
                                                                           n_fields,
                                                                           "Thermal conductivities");

          min_strain_rate = prm.get_double("Minimum strain rate");
          ref_strain_rate = prm.get_double("Reference strain rate");
          min_viscosity = prm.get_double ("Minimum viscosity");
          max_viscosity = prm.get_double ("Maximum viscosity");

          ref_viscosity = prm.get_double ("Reference viscosity");

          viscosity_averaging = MaterialUtilities::parse_compositional_averaging_operation ("Viscosity averaging scheme",
                                prm);

          linear_viscosities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Linear viscosities"))),
                                                                       n_fields,
                                                                       "Linear viscosities");

          angles_internal_friction = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Angles of internal friction"))),
                                                                             n_fields,
                                                                             "Angles of internal friction");
          cohesions = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Cohesions"))),
                                                              n_fields,
                                                              "Cohesions");

          maximum_yield_stress = prm.get_double("Maximum yield stress");

          strength_reductions = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Host rock strength reductions"))),
                                                                        n_fields,
                                                                        "Host rock strength reductions");

          melt_density_change        = prm.get_double ("Melt density change");
          xi_0                       = prm.get_double ("Reference bulk viscosity");
          eta_f                      = prm.get_double ("Reference melt viscosity");
          reference_permeability     = prm.get_double ("Reference permeability");
          alpha_phi                  = prm.get_double ("Exponential melt weakening factor");
          melting_time_scale         = prm.get_double ("Melting time scale for operator splitting");

          A1              = prm.get_double ("A1");
          A2              = prm.get_double ("A2");
          A3              = prm.get_double ("A3");
          B1              = prm.get_double ("B1");
          B2              = prm.get_double ("B2");
          B3              = prm.get_double ("B3");
          C1              = prm.get_double ("C1");
          C2              = prm.get_double ("C2");
          C3              = prm.get_double ("C3");
          r1              = prm.get_double ("r1");
          r2              = prm.get_double ("r2");
          beta            = prm.get_double ("beta");
          M_cpx           = prm.get_double ("Mass fraction cpx");

          use_melting_rate_function  = prm.get_bool ("Use melting rate function");

          if (use_melting_rate_function)
            {
              prm.enter_subsection("Melting rate function");
              try
                {
                  function.parse_parameters (prm);
                }
              catch (...)
                {
                  std::cerr << "ERROR: FunctionParser failed to parse\n"
                            << "\t'Material model.Melting rate function.Function'\n"
                            << "with expression\n"
                            << "\t'" << prm.get("Function expression") << "'"
                            << "More information about the cause of the parse error \n"
                            << "is shown below.\n";
                  throw;
                }
              prm.leave_subsection();
            }

          AssertThrow(this->introspection().compositional_name_exists("peridotite"),
                      ExcMessage("Material model Melt visco plastic only works if there is a "
                                 "compositional field called peridotite."));

          if (this->convert_output_to_years() == true)
            melting_time_scale *= year_in_seconds;

          AssertThrow(melting_time_scale >= this->get_parameters().reaction_time_step,
                      ExcMessage("The reaction time step " + Utilities::to_string(this->get_parameters().reaction_time_step)
                                 + " in the operator splitting scheme is too large to compute melting rates! "
                                 "You have to choose it in such a way that it is smaller than the 'Melting time scale for "
                                 "operator splitting' chosen in the material model, which is currently "
                                 + Utilities::to_string(melting_time_scale) + "."));
          AssertThrow(melting_time_scale > 0,
                      ExcMessage("The Melting time scale for operator splitting must be larger than 0!"));

          if (this->include_melt_transport())
            {
              AssertThrow(this->introspection().compositional_name_exists("porosity"),
                          ExcMessage("Material model Melt visco plastic with melt transport only "
                                     "works if there is a compositional field called porosity."));
            }

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

    }

    template <int dim>
    void
    MeltViscoPlastic<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      if (out.template get_additional_output<PlasticAdditionalOutputs2<dim>>() == nullptr)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std::make_unique<MaterialModel::PlasticAdditionalOutputs2<dim>> (n_points));
        }
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(MeltViscoPlastic,
                                   "EXPERIMENTAL melt visco plastic",
                                   "WARNING: This material model has not been thoroughly tested, "
                                   "and it is not advised that it be used for any serious models "
                                   "at this point.\n\n"
                                   "This material model implements a simple formulation of the "
                                   "material parameters required for the modeling of melt transport, "
                                   "including a source term for the porosity according to the melting "
                                   "model for dry peridotite of \\cite{KSL2003}. All other material "
                                   "properties are taken from the visco-plastic model.")
  }
}
