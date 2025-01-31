/*
  Copyright (C) 2015 - 2023 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_melt_boukare_h
#define _aspect_material_model_melt_boukare_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/melt.h>

namespace aspect
{
  namespace MaterialModel
  {
    /**
     * Additional output fields for the melt boukare material model.
     */
    template <int dim>
    class BoukareOutputs : public NamedAdditionalMaterialOutputs<dim>
    {
      public:
        BoukareOutputs(const unsigned int n_points);

        std::vector<double> get_nth_output(const unsigned int idx) const override;

        /**
         * Bulk composition of the material.
         */
        std::vector<double> bulk_composition;
        std::vector<double> molar_volatiles_in_melt;
    };

    /**
     * A material model that implements a simplified version of the melting
     * model of Boukare et al. (https://doi.org/10.1002/2015JB011929) for the
     * lowermost mantle. The model parameterizes the composition (which includes
     * the components MgO, FeO and SiO2) as a mixture between two endmembers
     * (one iron-bearing and one magnesium-bearing). The equation of state
     * considers three phases: bridgmanite, ferropericlase, and melt (each with
     * their individual compositions).
     * More details can be found in Dannberg, J., Myhill, R., Gassmöller, R.,
     * & Cottaar, S. (2021). The morphology, evolution and seismic visibility
     * of partial melt at the core–mantle boundary: implications for ULVZs.
     * Geophysical Journal International, 227(2), 1028-1059.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class MeltBoukare : public MaterialModel::MeltInterface<dim>,
      public MaterialModel::MeltFractionModel<dim>,
      public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Initialization function. Computes endmember properties.
         */
        void
        initialize () override;

        /**
         * Return whether the model is compressible or not. In this model,
         * both the melt and solid are compressible.
         */
        bool is_compressible () const override;

        void evaluate(const typename Interface<dim>::MaterialModelInputs &in,
                      typename Interface<dim>::MaterialModelOutputs &out) const override;

        /**
         * Compute the equilibrium melt fractions for the given input conditions.
         * @p in and @p melt_fractions need to have the same size.
         *
         * @param in Object that contains the current conditions.
         * @param melt_fractions Vector of doubles that is filled with the
         * equilibrium melt fraction for each given input conditions.
         */
        void melt_fractions (const MaterialModel::MaterialModelInputs<dim> &in,
                             std::vector<double> &melt_fractions) const override;

        /**
         * @name Reference quantities
         * @{
         */
        double reference_darcy_coefficient () const override;

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
        void
        parse_parameters (ParameterHandler &prm) override;

        /**
         * @}
         */

        void
        create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const override;


      private:
        // Thermodynamic parameters of the endmember components in the equation of state,
        // needed to compute their material properties.
        std::vector<std::string> endmember_names;
        std::vector<double> molar_masses;
        std::vector<double> number_of_atoms;
        std::vector<double> reference_volumes;
        std::vector<double> reference_thermal_expansivities;
        std::vector<double> reference_bulk_moduli;
        std::vector<double> bulk_modulus_pressure_derivatives;
        std::vector<double> bulk_modulus_second_pressure_derivatives;
        std::vector<double> Einstein_temperatures;
        std::vector<double> reference_enthalpies;
        std::vector<double> reference_entropies;
        std::vector<double> reference_specific_heats;
        std::vector<double> specific_heat_linear_coefficients;
        std::vector<double> specific_heat_second_coefficients;
        std::vector<double> specific_heat_third_coefficients;

        std::vector<double> tait_parameters_a;
        std::vector<double> tait_parameters_b;
        std::vector<double> tait_parameters_c;

        // Reference conditions for computing the material properties in the equation of state.
        double reference_temperature;
        double reference_pressure;

        // Rheological properties.
        double eta_0;
        double xi_0;
        double eta_f;
        double thermal_viscosity_exponent;
        double thermal_bulk_viscosity_exponent;
        double alpha_phi;

        // Transport properties that are not computed from the equation of state.
        double thermal_conductivity;
        double reference_permeability;

        // Parameters controlling melting and solidification of melt.
        bool include_melting_and_freezing;
        double melting_time_scale;

        // Parameters defining the molar composition of the two endmember compositions we use
        // in the melting model.
        const double molar_MgO_in_Mg_mantle_endmember = 0.581;
        const double molar_SiO2_in_Mg_mantle_endmember = 0.419;
        const double molar_FeO_in_Fe_mantle_endmember = 0.908;
        const double molar_SiO2_in_Fe_mantle_endmember = 0.092;

        // Parameters describing the melting properties of the two endmembers of the melting model.
        const double melting_reference_pressure = 120.e9;

        // reference melting temperature for Fe and Mg mantle endmember at the reference pressure
        double Fe_mantle_melting_temperature;
        double Mg_mantle_melting_temperature;

        // molar entropy change of melting in J/mol K
        const double Fe_mantle_melting_entropy = 33.77;
        const double Mg_mantle_melting_entropy = 34.33;

        // molar volume change of melting of solid Fe and Mg mantle endmember in m3/mol
        const double Fe_mantle_melting_volume = 1.51e-07;
        const double Mg_mantle_melting_volume = 9.29e-08;

        // Number of moles of atoms mixing on pseudosite in mantle lattice (empirical model fitting the full Boukare model).
        double Fe_number_of_moles;
        double Mg_number_of_moles;

        // Indices for the endmembers in the equation of state.
        unsigned int febdg_idx;
        unsigned int mgbdg_idx;
        unsigned int wus_idx;
        unsigned int per_idx;
        unsigned int femelt_idx;
        unsigned int mgmelt_idx;


        // Parameter describing if an endmember in the equation of state is solid or molten.
        struct EndmemberState
        {
          enum Kind
          {
            solid,
            melt
          };
        };

        std::vector<typename EndmemberState::Kind> endmember_states;


        // Material properties of the different endmembers of the equation of state.
        struct EndmemberProperties
        {
          /**
           * Constructor. Initialize the various arrays of this structure with the
           * given number of endmembers.
           */
          EndmemberProperties(const unsigned int n_endmembers);

          std::vector<double> volumes;
          std::vector<double> gibbs_energies;
          std::vector<double> entropies;
          std::vector<double> thermal_expansivities;
          std::vector<double> bulk_moduli;
          std::vector<double> heat_capacities;
        };


        /**
         * Fill the endmember properties at a single point.
         */
        virtual
        void
        fill_endmember_properties (const typename Interface<dim>::MaterialModelInputs &in,
                                   const unsigned int q,
                                   EndmemberProperties &properties) const;


        /**
         * Calculate the Einstein thermal energy of an endmember component.
         */
        virtual
        double
        endmember_thermal_energy (const double temperature,
                                  const unsigned int endmember_index) const;


        /**
         * Calculate the heat capacity of an endmember component.
         */
        virtual
        double
        endmember_molar_heat_capacity (const double temperature,
                                       const unsigned int endmember_index) const;

        /**
         * Calculate the thermal pressure of an endmember component.
         */
        virtual
        double
        endmember_thermal_pressure (const double temperature,
                                    const unsigned int endmember_index) const;

        /**
         * Calculate the thermal addition to the standard state enthalpy of an endmember component
         * at the reference pressure.
         */
        virtual
        double
        endmember_enthalpy_thermal_addition (const double temperature,
                                             const unsigned int endmember_index) const;

        /**
         * Calculate the thermal addition to the standard state enttropy of an endmember component
         * at the reference pressure.
         */
        virtual
        double
        endmember_entropy_thermal_addition (const double temperature,
                                            const unsigned int endmember_index) const;


        /**
         * Convert from the mole fraction of iron in the solid and melt to the mole fraction of iron
         * and magnesium in each of the two solid phases, bridgmanite and ferropericlase, and the melt
         * phase. In addition, compute the molar fraction of bridgmanite in the solid.
         */
        virtual
        void
        convert_composition_to_fraction_of_endmembers (const double temperature,
                                                       const double molar_Fe_in_solid,
                                                       const double molar_Fe_in_melt,
                                                       const std::vector<double> &endmember_gibbs_energies,
                                                       std::vector<double> &endmember_mole_fractions_per_phase,
                                                       double &molar_bridgmanite_in_solid) const;


        /**
         * Compute the molar fraction of melt from the volume fraction of melt (the porosity).
         * This requires knowing the densities of both solid and melt, which can be computed
         * from the fractions and material properties of the component endmembers.
         */
        virtual
        double
        compute_melt_molar_fraction (const double porosity,
                                     const double bridgmanite_molar_fraction_in_solid,
                                     EndmemberProperties &properties,
                                     const std::vector<double> &endmember_mole_fractions_per_phase) const;


        /**
         * Make sure that when the property @p value is updated by adding @p change_of_value,
         * the new value is between zero and one. Otherwise, throw an error message.
         */
        virtual
        double
        assert_update_is_within_0_and_1 (const double value,
                                         const double change_of_value) const;


        /**
         * Compute the equilibrium melt fraction for a given @p temperature, @p pressure,
         * and @p bulk_composition and the corresponding composition of melt and solid. The
         * melting model is based on Phipps Morgan, Jason. "Thermodynamics of pressure release
         * melting of a veined plum pudding mantle." Geochemistry, Geophysics, Geosystems 2.4
         * (2001), and additionally includes volatiles.
         */
        virtual
        double
        melt_fraction (const double temperature,
                       const double pressure,
                       const double bulk_composition,
                       double &molar_volatiles_in_melt,
                       double &solid_composition,
                       double &melt_composition) const;
    };

  }
}

#endif
