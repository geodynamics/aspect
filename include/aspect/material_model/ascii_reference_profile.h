/*
  Copyright (C) 2016 - 2022 by the authors of the ASPECT code.

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


#ifndef _aspect_model_ascii_reference_profile_h
#define _aspect_model_ascii_reference_profile_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>
namespace aspect
{
  namespace MaterialModel
  {
    /**
     * A material model that reads in a reference profile from an ascii file
     * and computes properties based on this profile.
     *
     * The viscosity is computed as
     * \f[
     * \eta(z,T) = \eta_r(z) \eta_0 \exp\left(-A \frac{T - T_\text{adi}}{T_\text{adi}}\right)."
     * \f]
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class AsciiReferenceProfile : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor. Initialize variables.
         */
        AsciiReferenceProfile ();

        /**
         * Initialization function. This function is called once at the
         * beginning of the program after parse_parameters is run and after
         * the SimulatorAccess (if applicable) is initialized.
         */
        void
        initialize () override;

        void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                      MaterialModel::MaterialModelOutputs<dim> &out) const override;

        /**
         * @name Qualitative properties one can ask a material model
         * @{
         */

        /**
         * Return whether the model is compressible or not.
         */
        bool is_compressible () const override;
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

        /**
         * Add the named outputs for seismic velocities.
         */
        void
        create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const override;


      private:
        /**
         * Use truncated anelastic approximation?
         */
        bool tala;

        /**
         * The reference viscosity
         */
        double viscosity;

        /**
         * The constant $A$ in the temperature dependence of viscosity
         * $\exp\left(-A \frac{T - T_\text{adi}}{T_\text{adi}}\right).$
         */
        double thermal_viscosity_exponent;

        /**
         * A list of constants that make up the piece-wise constant function
         * $\eta_r(z)$, which determines the depth dependence of viscosity,
         * and is multiplied with the reference viscosity and the
         * temperature dependence to compute the viscosity $\eta(z,T)$.
         */
        std::vector<double> viscosity_prefactors;

        /**
         * A list of depths that determine the locations of the jumps in
         * the piece-wise constant function $\eta_r(z)$, which describes the
         * depth dependence of viscosity.
         */
        std::vector<double> transition_depths;

        double thermal_conductivity;

        /**
         * Object containing the data profile.
         */
        aspect::Utilities::AsciiDataProfile<dim> profile;

        /**
         * The column indices of the temperature, pressure, and density column
         * in the data file.
         */
        unsigned int density_index;
        unsigned int thermal_expansivity_index;
        unsigned int specific_heat_index;
        unsigned int compressibility_index;

        /**
         * The column indices of the seismic velocities and their temperature
         * derivatives columns in the data file.
         */
        unsigned int seismic_vp_index;
        unsigned int seismic_vs_index;
        unsigned int seismic_dvp_dT_index;
        unsigned int seismic_dvs_dT_index;
    };

  }
}

#endif
