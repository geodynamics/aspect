/*
  Copyright (C) 2019 - 2024 by the authors of the ASPECT code.

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

#include <deal.II/matrix_free/fe_point_evaluation.h>

namespace aspect
{
  namespace MaterialModel
  {
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
          fill_elastic_outputs (const MaterialModel::MaterialModelInputs<dim> &in,
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
           * Calculates the effective elastic viscosity (this is the equivalent viscosity of
           * a material which was unstressed at the end of the previous timestep).
           */
          double
          calculate_elastic_viscosity (const double shear_modulus) const;

          /**
           * Given the (viscous or visco-plastic) viscosity and the shear modulus, compute the viscoelastic
           * viscosity (eqn 28 in Moresi et al., 2003, J. Comp. Phys.).
           */
          double
          calculate_viscoelastic_viscosity (const double viscosity,
                                            const double shear_modulus) const;

          /**
           * Calculate the effective deviatoric strain rate tensor,
           * which equals the true deviatoric strain rate plus
           * a fictional strain rate which would arise from stored elastic stresses.
           * In ASPECT, this additional strain rate is
           * supported by a fictional body force.
           * This formulation allows the use of an isotropic effective viscosity
           * by ensuring that the resulting strain rate tensor is equal to the
           * total current stress tensor multiplied by a scalar.
           */
          SymmetricTensor<2,dim>
          calculate_viscoelastic_strain_rate (const SymmetricTensor<2,dim> &strain_rate,
                                              const SymmetricTensor<2,dim> &stored_stress,
                                              const double shear_modulus) const;

          /**
           * Compute the elastic time step.
           */
          double
          elastic_timestep () const;

        private:
          /**
           * Viscosity of a damper used to stabilize elasticity.
           * A value of 0 Pas is equivalent to not using a damper.
           */
          double elastic_damper_viscosity;

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

          /**
           * We cache the evaluator that is necessary to evaluate the old velocity
           * gradients. They are required to compute the elastic stresses, but
           * are not provided by the material model.
           * By caching the evaluator, we can avoid recreating it every time we
           * need it.
           */
          mutable std::unique_ptr<FEPointEvaluation<dim, dim>> evaluator;
      };
    }
  }
}
#endif
