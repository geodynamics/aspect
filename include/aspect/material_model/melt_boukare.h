/*
  Copyright (C) 2015 - 2019 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_melt_global_h
#define _aspect_material_model_melt_global_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/postprocess/melt_statistics.h>
#include <aspect/melt.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * Additional output fields for the melt boukare material model.
     */
    template <int dim>
    class BoukareOutputs : public NamedAdditionalMaterialOutputs<dim>
    {
      public:
        BoukareOutputs(const unsigned int n_points);

        virtual std::vector<double> get_nth_output(const unsigned int idx) const;

        /**
         * Bulk composition of the material.
         */
        std::vector<double> bulk_composition;
        std::vector<double> molar_volatiles_in_melt;
    };

    /**
     * A material model that implements a simple formulation of the
     * material parameters required for the modelling of melt transport
     * in a global model, including a source term for the porosity according
     * a simplified linear melting model.
     *
     * The model is considered incompressible, following the definition
     * described in Interface::is_compressible.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class MeltBoukare : public MaterialModel::MeltInterface<dim>, public ::aspect::SimulatorAccess<dim>, public MaterialModel::MeltFractionModel<dim>
    {
      public:
        /**
          * Initialization function. Computes endmember properties.
          */
        virtual
        void
        initialize ();

        /**
         * Return whether the model is compressible or not.  Incompressibility
         * does not necessarily imply that the density is constant; rather, it
         * may still depend on temperature or pressure. In the current
         * context, compressibility means whether we should solve the continuity
         * equation as $\nabla \cdot (\rho \mathbf u)=0$ (compressible Stokes)
         * or as $\nabla \cdot \mathbf{u}=0$ (incompressible Stokes).
         */
        virtual bool is_compressible () const;

        virtual void evaluate(const typename Interface<dim>::MaterialModelInputs &in,
                              typename Interface<dim>::MaterialModelOutputs &out) const;

        /**
         * Compute the equilibrium melt fractions for the given input conditions.
         * @p in and @p melt_fractions need to have the same size.
         *
         * @param in Object that contains the current conditions.
         * @param melt_fractions Vector of doubles that is filled with the
         * equilibrium melt fraction for each given input conditions.
         */
        virtual void melt_fractions (const MaterialModel::MaterialModelInputs<dim> &in,
                                     std::vector<double> &melt_fractions) const;

        /**
         * @name Reference quantities
         * @{
         */
        virtual double reference_viscosity () const;

        virtual double reference_darcy_coefficient () const;


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
        create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const;


      private:
        // properties of the endmember components
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

        double reference_temperature;
        double reference_pressure;

        double eta_0;
        double xi_0;
        double eta_f;
        double thermal_viscosity_exponent;
        double thermal_bulk_viscosity_exponent;

        double thermal_conductivity;
        double reference_permeability;
        double alpha_phi;

        bool include_melting_and_freezing;
        double melting_time_scale;

        // melting model parameters
        const double molar_MgO_in_Mg_mantle_endmember = 0.581;
        const double molar_SiO2_in_Mg_mantle_endmember = 0.419;
        const double molar_FeO_in_Fe_mantle_endmember = 0.908;
        const double molar_SiO2_in_Fe_mantle_endmember = 0.092;

        // number of moles of atoms mixing on pseudosite in mantle lattice (empirical model fitting the full boukare model)
        double Fe_number_of_moles;
        double Mg_number_of_moles;


        // names of the endmembers
        unsigned int febdg_idx;
        unsigned int mgbdg_idx;
        unsigned int wus_idx;
        unsigned int per_idx;
        unsigned int femelt_idx;
        unsigned int mgmelt_idx;
        unsigned int simelt_idx;

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
         * Fill the endmember properties at a single quadrature point.
         */
        virtual
        void
        fill_endmember_properties (const typename Interface<dim>::MaterialModelInputs &in,
                                   const unsigned int q,
                                   EndmemberProperties &properties) const;


        struct EndmemberState
        {
          enum Kind
          {
            solid,
            melt
          };
        };

        typename EndmemberState::Kind density_formulation;
        std::vector<typename EndmemberState::Kind> endmember_states;

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
         * Convert from the mole fraction of iron in the solid to the mole fraction of iron in the
         * two solid phases, bridgmanite and ferropericlase, and the molar fraction of bridgmanite
         * in the solid.
         */
        virtual
        void
        convert_to_fraction_of_endmembers_in_solid (const double temperature,
                                                    const double pressure,
                                                    const double molar_Fe_in_solid,
                                                    const std::vector<double> &endmember_gibbs_energies,
                                                    double &molar_FeSiO3_in_bridgmanite,
                                                    double &molar_FeO_in_ferropericlase,
                                                    double &molar_bridgmanite_in_solid) const;

        /**
         * Convert from the mole fraction of iron in the solid to the mole fraction of iron in the
         * two solid phases, bridgmanite and ferropericlase, and the mass fraction of bridgmanite
         * in the solid.
         */
        virtual
        double
        compute_melt_molar_fraction (const double porosity,
                                     const double bridgmanite_molar_fraction_in_solid,
                                     EndmemberProperties &properties,
                                     const std::vector<double> &endmember_mole_fractions_per_phase) const;

        virtual
        double
        melt_fraction (const double temperature,
                       const double pressure,
                       const double bulk_composition,
                       double &molar_volatiles_in_melt,
                       double &solid_composition,
                       double &melt_composition) const;

        virtual
        double
        limit_update_to_0_and_1 (const double value,
                                 const double change_of_value) const;


        const double melting_reference_pressure = 120.e9;       // Pa

        double Fe_mantle_melting_temperature = 3424.5;          // Kelvin at the reference pressure - reference melting temperature for Fe mantle endmember
        double Mg_mantle_melting_temperature = 4821.2;          // Kelvin at reference pressure - reference melting temperature for Mg mantle endmember

        const double Fe_mantle_melting_entropy = 33.77;         // molar entropy change of melting in J/mol K
        const double Mg_mantle_melting_entropy = 34.33;         // molar entropy change of melting in J/mol K

        const double Fe_mantle_melting_volume = 1.51e-07;       // molar volume change of melting of solid Fe mantle endmember in m3/mol
        const double Mg_mantle_melting_volume = 9.29e-08;       // molar volume change of melting volume of solid Mg mantle endmember in m3/mol
    };

  }
}

#endif
