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

#ifndef _aspect_material_model_dynamic_friction_h
#define _aspect_material_model_dynamic_friction_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/material_model/rheology/drucker_prager.h>
#include <aspect/material_model/equation_of_state/multicomponent_incompressible.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * This model is for use with an arbitrary number of compositional fields, where each field
     * represents a rock type which can have completely different properties from the others.
     * Each rock type itself has constant material properties, with the exception of viscosity
     * which is modified according to a Drucker-Prager yield criterion. Unlike the drucker prager
     * or visco plastic material models, the angle of internal friction is a function of velocity.
     * This relationship is similar to rate-and-state friction constitutive relationships, which
     * are applicable to the strength of rocks during earthquakes. The formulation used here is
     * derived from van Dinther et al. 2013, JGR.

     * For each material parameter the user supplies a comma delimited list of
     * length N+1, static friction of coefficient, dynamic friction of coefficient,
     * cohesions and background viscosity (to calculate viscous stresses )where N
     * is the number of compositional fields.  The
     * additional field corresponds to the value for background mantle.  They
     * should be ordered ``background, composition1, composition2...''
     *
     * If a single value is given, then all the compositional fields are given
     * that value. Other lengths of lists are not allowed.  For a given
     * compositional field the material parameters are treated as constant,
     * except density, which varies linearly with temperature according to the
     * thermal expansivity.
     *
     * When more than one field is present at a point, they are averaged
     * arithmetically. An exception is viscosity, which may be averaged
     * arithmetically, harmonically, geometrically, or by selecting the
     * viscosity of the composition with the greatest volume fraction.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class DynamicFriction : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
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
         * This model is not compressible, so this returns false.
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
         * From a list of static friction of coefficient, dynamic friction of
         * coefficient, cohesions and background viscosity for N + 1 fields
         * (background mantle and N compositions) , we compute viscosities for
         * drucker prager model with coefficient of friction dependent on the
         * strain rate.
         */
        const std::vector<double> compute_viscosities(const double pressure,
                                                      const SymmetricTensor<2,dim> &strain_rate) const;

        /**
         * Reference temperature for thermal expansion.  All components use
         * the same reference_T.
         */
        double reference_T;

        MaterialUtilities::CompositionalAveragingOperation viscosity_averaging;

        EquationOfState::MulticomponentIncompressible<dim> equation_of_state;

        /*
          * Objects for computing plastic stresses, viscosities, and additional outputs
          */
        Rheology::DruckerPrager<dim> drucker_prager_plasticity;

        /**
         * The dynamic coefficient of friction
         */
        std::vector<double> mu_d;

        /**
         * The applied viscosity bounds
         */
        double minimum_viscosity;
        double maximum_viscosity;

        /**
         * The static coefficient of friction
         */
        std::vector<double> mu_s;

        /**
         * Vector for field viscosities, read from parameter file.
         */
        std::vector<double> cohesions;

        /**
         * Vector for field viscosities, read from parameter file.
         */
        std::vector<double> background_viscosities;

        /**
         * The reference strain rate used as a first estimate
         */
        double reference_strain_rate;
        double minimum_strain_rate;

        /**
         * Vector for field thermal conductivities, read from parameter file.
         */
        std::vector<double> thermal_conductivities;
    };

  }
}

#endif
