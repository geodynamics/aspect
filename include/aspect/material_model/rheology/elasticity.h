/*
  Copyright (C) 2019 - 2020 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_rheology_elasticity_h
#define _aspect_material_model_rheology_elasticity_h

#include <aspect/global.h>
#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * Additional output fields for the elastic shear modulus to be added to
     * the MaterialModel::MaterialModelOutputs structure and filled in the
     * MaterialModel::Interface::evaluate() function.
     */
    template <int dim>
    class ElasticAdditionalOutputs : public NamedAdditionalMaterialOutputs<dim>
    {
      public:
        explicit ElasticAdditionalOutputs(const unsigned int n_points);

        std::vector<double> get_nth_output(const unsigned int idx) const override;

        /**
         * Elastic shear moduli at the evaluation points passed to
         * the instance of MaterialModel::Interface::evaluate() that fills
         * the current object.
         */
        std::vector<double> elastic_shear_moduli;
    };



    namespace Rheology
    {
      template <int dim>
      class Elasticity : public ::aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * Declare the parameters this function takes through input files.
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Read the parameters from the parameter file.
           */
          void
          parse_parameters (ParameterHandler &prm);

          /**
           * Create the additional material model outputs object that contains the
           * elastic shear moduli.
           */
          void
          create_elastic_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const;

          /**
           * Given the stress of the previous time step in the material model inputs @p in,
           * the elastic shear moduli @p average_elastic_shear_moduli a each point,
           * and the (viscous) viscosities given in the material model outputs object @p out,
           * fill an additional material model outputs objects with the elastic force terms.
           */
          void
          fill_elastic_force_outputs (const MaterialModel::MaterialModelInputs<dim> &in,
                                      const std::vector<double> &average_elastic_shear_moduli,
                                      MaterialModel::MaterialModelOutputs<dim> &out) const;

          /**
           * Given the stress of the previous time step in the material model inputs @p in,
           * the elastic shear moduli @p average_elastic_shear_moduli a each point,
           * and the (viscous) viscosities given in the material model outputs object @p out,
           * compute an update to the elastic stresses and use it to fill the reaction terms
           * material model output property.
           */
          void
          fill_reaction_outputs (const MaterialModel::MaterialModelInputs<dim> &in,
                                 const std::vector<double> &average_elastic_shear_moduli,
                                 MaterialModel::MaterialModelOutputs<dim> &out) const;

          /**
           * Return the values of the elastic shear moduli for each composition used in the
           * rheology model.
           */
          const std::vector<double> &
          get_elastic_shear_moduli () const;

          /**
           * Given the (viscous or visco-plastic) viscosity and the shear modulus, compute the viscoelastic
           * viscosity (eqn 28 in Moresi et al., 2003, J. Comp. Phys.).
           */
          double
          calculate_viscoelastic_viscosity (const double viscosity,
                                            const double shear_modulus) const;

          /**
           * Calculate the square root of the second moment invariant for the deviatoric
           * strain rate tensor, including viscoelastic stresses.
           */
          double
          calculate_viscoelastic_strain_rate (const SymmetricTensor<2,dim> &strain_rate,
                                              const SymmetricTensor<2,dim> &stress,
                                              const double shear_modulus) const;

          /**
           * Compute the elastic time step.
           */
          double
          elastic_timestep () const;

        private:
          /**
           * Vector for field elastic shear moduli, read from parameter file.
           */
          std::vector<double> elastic_shear_moduli;

          /**
           * Bool indicating whether to use a fixed material time scale in the
           * viscoelastic rheology for all time steps (if true) or to use the
           * actual (variable) advection time step of the model (if false). Read
           * from parameter file.
           */
          bool use_fixed_elastic_time_step;

          /**
           * Bool indicating whether to use a stress averaging scheme to account
           * for differences between the numerical and fixed elastic time step
           * (if true). When set to false, the viscoelastic stresses are not
           * modified to account for differences between the viscoelastic time
           * step and the numerical time step. Read from parameter file.
           */
          bool use_stress_averaging;

          /**
           * Double for fixed elastic time step value, read from parameter file
           */
          double fixed_elastic_time_step;

          /**
           * A stabilization factor for the elastic stresses that influences how 
           * fast elastic stresses adjust to deformation. 1.0 is equivalent to no 
           * stabilization, and infinity is equivalent to not applying elastic 
           * stresses at all. The factor is multiplied with the computational
           * time step to create a time scale. 
           */
          double stabilization_time_scale_factor;
      };
    }
  }
}
#endif
