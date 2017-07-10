/*
  Copyright (C) 2016 - 2017 by the authors of the ASPECT code.

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


#ifndef _aspect_material_model_nondimensional_h
#define _aspect_material_model_nondimensional_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * A material model for incompressible (using the Boussinesq approximation)
     * and compressible computations (with ALA or TALA) for a nondimensionalized
     * problem. The viscosity is (optionally) depth and temperature dependent.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class Nondimensional : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:

        /**
         * Initialize.
         */
        virtual
        void initialize ();

        /**
         * Evaluate the material model.
         */
        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const;

        /**
         * @name Qualitative properties one can ask a material model
         * @{
         */

        /**
         * Return whether the model is compressible or not.
         */
        virtual bool is_compressible () const;
        /**
         * @}
         */

        /**
         * @name Reference quantities
         * @{
         */
        virtual double reference_viscosity () const;
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

      private:
        /**
         * Use truncated anelastic approximation?
         */
        bool tala;

        /**
         * Are we using ALA/TALA or BA? (determined by the choice of Di)
         */
        bool compressible;

        /**
         * Parameter describing the temperature prefactor of viscosity.
         */
        double exponential_viscosity_temperature_prefactor;

        /**
         * Parameter describing the depth prefactor of viscosity.
         */
        double exponential_viscosity_depth_prefactor;

        /**
         * The surface density.
         */
        double reference_rho;

        /**
         * The nondimensional numbers (Dissipation number,
         * Rayleigh number, Grueneisen parameter).
         */
        double Di, Ra, gamma;

        /**
         * The constant specific heat
         */
        double reference_specific_heat;

    };

  }
}

#endif
