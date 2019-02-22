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
          * Initialization function. Computes endmember propoerties.
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

        // entropy change upon melting
        double peridotite_melting_entropy_change;

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

        virtual
        double
        melt_fraction (const double temperature,
                       const double pressure,
                       const double bulk_composition,
                       double &solid_composition,
                       double &melt_composition) const;

        virtual
        double
        limit_update_to_0_and_1 (const double value,
                                 const double change_of_value) const;
    };

  }
}

#endif
