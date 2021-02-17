/*
  Copyright (C) 2014 - 2019 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_multicomponent_compressible_h
#define _aspect_material_model_multicomponent_compressible_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/material_model/equation_of_state/multicomponent_compressible.h>


namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * A material model which is intended for use with multiple compositional
     * fields. Each compositional field is meant to be a single rock type,
     * where the value of the field at a point is interpreted to be the
     * mass fraction of that rock type. If the sum of the compositional field
     * mass fractions is less than one, then the remainder of the mass is
     * assumed to be ``background mantle''.  If the sum of the compositional
     * field mass fractions is greater than one, then they are renormalized
     * to sum to one and there is no background mantle.
     *
     * For each material parameter the user supplies a comma delimited list of
     * length N+1, where N is the number of compositional fields.  The
     * additional field corresponds to the value for background mantle.  They
     * should be ordered ``background, composition1, composition2...''
     *
     * If a single value is given, then all the compositional fields are given
     * that value. Other lengths of lists are not allowed.  The material
     * parameters for each compositional field are derived from the
     * multicomponent compressible equation of state, and are pressure and
     * temperature dependent.
     *
     * When more than one field is present at a point, they are either averaged
     * self-consistently (for equation of state properties) or arithmetically
     * (for thermal conductivity). For the special case of viscosity, the user
     * may select from a variety of averaging schemes:
     * arithmetic, harmonic, geometric, or by selecting the
     * viscosity of the composition with the greatest volume fraction.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class MulticomponentCompressible : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:

        /**
         * Function to compute the material properties in @p out given the
         * inputs in @p in. If MaterialModelInputs.strain_rate has the length
         * 0, then the viscosity does not need to be computed.
         */
        void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                      MaterialModel::MaterialModelOutputs<dim> &out) const override;

        /**
         * @name Qualitative properties one can ask a material model
         * @{
         */

        /**
         * This model is compressible, so this returns true.
         */
        bool is_compressible () const override;
        /**
         * @}
         */

        /**
         * @name Reference quantities
         * @{
         */
        double reference_viscosity () const override;
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

      private:
        /**
         * Enumeration for selecting which viscosity averaging scheme to use.
         */
        MaterialUtilities::CompositionalAveragingOperation viscosity_averaging;

        /**
         * Vector for field viscosities, read from parameter file.
         */
        std::vector<double> viscosities;

        /**
         * Vector for field thermal conductivities, read from parameter file.
         */
        std::vector<double> thermal_conductivities;

        EquationOfState::MulticomponentCompressible<dim> equation_of_state;
    };

  }
}

#endif
