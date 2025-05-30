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
     * Additional output fields for the elastic shear modulus and other
     * elastic outputs to be added to the MaterialModel::MaterialModelOutputs
     * structure and filled in the MaterialModel::Interface::evaluate() function.
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

        /**
        * Elastic viscosity at the evaluation points passed to
        * the instance of MaterialModel::Interface::evaluate() that fills
        * the current object.
        */
        std::vector<double> elastic_viscosity;

        /**
        * The deviatoric stress of the current timestep, so including
        * the rotation, advection and stress update, at the evaluation points
        * passed to the instance of MaterialModel::Interface::evaluate()
        * that fills the current object.
        */
        std::vector<SymmetricTensor<2,dim>> deviatoric_stress;
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
           * Create the two additional material model output objects that contain the
           * elastic shear moduli, elastic viscosity, ratio of computational to elastic timestep,
           * and deviatoric stress of the current timestep and the reaction rates.
           */
          void
          create_elastic_additional_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const;

          /**
           * Given the stress of the previous time step in the material model inputs @p in,
           * the elastic shear moduli @p average_elastic_shear_moduli at each point,
           * and the (viscous) viscosities given in the material model outputs object @p out,
           * fill a material model outputs objects with the elastic force terms, viscoelastic
           * strain rate and viscous dissipation.
           *
           * Two sets of stresses are available from @p in:
           * 1) the stress tensor components stress_0_advected, which represent the stress from the previous
           * timestep $t$ rotated and advected into the current timestep $t+\Delta t_c$; and
           * 2) the stress tensor components stress_old, which represent the stress from the previous
           * timestep $t$ advected into the current timestep $t+\Delta t_c$.
           * Rotation of the stresses has been applied through the reaction_terms filled in the function
           * fill_reaction_outputs that are used in the right hand side of the respective compositional
           * field advection equations or to update the particles after they have been advected.
           * Advection of the stresses occurs through solving the respective field advection equations,
           * or by advecting the particles carrrying the stresses.
           * By the time the elastic force terms and the viscoelastic strain rate are required to assemble
           * the Stokes system, the stresses in @p in have thus been rotated and/or advected.
           */
          void
          fill_elastic_outputs (const MaterialModel::MaterialModelInputs<dim> &in,
                                const std::vector<double> &average_elastic_shear_moduli,
                                MaterialModel::MaterialModelOutputs<dim> &out) const;

          /**
           * Given the stress of the previous time step in the material model inputs @p in,
           * the elastic shear moduli @p average_elastic_shear_moduli at each point,
           * and the (viscous) viscosities given in the material model outputs object @p out,
           * fill the additional material model outputs (ElasticAdditionalOutputs) in @p out with the
           * average shear modulus, elastic viscosity, and the deviatoric stress of the current timestep.
           *
           * Two sets of stresses are available from @p in:
           * 1) the stress tensor components stress_0_advected, which represent the stress from the previous
           * timestep $t$ rotated and advected into the current timestep $t+\Delta t_c$; and
           * 2) the stress tensor components stress_old, which represent the stress from the previous
           * timestep $t$ advected into the current timestep $t+\Delta t_c$.
           * The stress update to the total deviatoric stress of timestep $t+\Delta t_c$ only occurs
           * at the beginning of the next timestep through an operator splitting procedure if the
           * stresses are tracked on compositional fields or through a direct update of the stresses
           * stored on particles. In both cases, this is done through the reaction_rates computed in
           * the function fill_reaction_rates. In case the full deviatoric stress is needed earlier, e.g.
           * during postprocessing of timestep $t+\Delta t_c$, this function computes it.
           */
          void
          fill_elastic_additional_outputs (const MaterialModel::MaterialModelInputs<dim> &in,
                                           const std::vector<double> &average_elastic_shear_moduli,
                                           MaterialModel::MaterialModelOutputs<dim> &out) const;

          /**
           * Given the stress of the previous time step in the material model inputs @p in,
           * the elastic shear moduli @p average_elastic_shear_moduli at each point,
           * and the (viscous) viscosities given in the material model outputs object @p out,
           * compute an update to the elastic stresses and use it to fill the reaction terms
           * material model output property.
           *
           * The reaction terms are used to applied the rotation of the stresses from the previous
           * timestep $t$ to the current timestep $t+\Delta t_c$. As such, the reaction_terms are
           * only non-zero for the first set of stresses that are tracked on compositional fields
           * or particles.
           * The reaction terms are an update to the stresses, requiring subtracting the old stresses.
           * In case of compositional fields, this requires evaluating the solution of the previous
           * timestep, which has been updated to the full deviatoric stress at the beginning of the
           * current timestep through operator splitting, instead of the first set of stresses in @p in
           * (which are the current linearization point).
           * When using particles, however, the first set of stresses in @p in represents the
           * full deviatoric stress of the last timestep.
           */
          void
          fill_reaction_outputs (const MaterialModel::MaterialModelInputs<dim> &in,
                                 const std::vector<double> &average_elastic_shear_moduli,
                                 MaterialModel::MaterialModelOutputs<dim> &out) const;

          /**
           * Given the stress of the previous time step in the material model inputs @p in,
           * the elastic shear moduli @p average_elastic_shear_moduli at each point,
           * and the (viscous) viscosities given in the material model outputs object @p out,
           * compute the update to the elastic stresses of the previous timestep and use it
           * to fill the reaction rates material model output property in @p out.
           *
           * Both sets of stresses tracked on compositional fields or particles are updated
           * through the reaction rates. They are updated to the same total deviatoric stress.
           * Later in the timestep, however, the first set will be rotated and advected,
           * the second only advected. The stresses in @p in are based on 'solution' or
           * 'old_solution', which at the time of operator splitting are the same.
           */
          void
          fill_reaction_rates (const MaterialModel::MaterialModelInputs<dim> &in,
                               const std::vector<double> &average_elastic_shear_moduli,
                               MaterialModel::MaterialModelOutputs<dim> &out) const;

          /**
           * Return the values of the elastic shear moduli for each composition used in the
           * rheology model.
           */
          const std::vector<double> &
          get_elastic_shear_moduli () const;

          /**
           * Calculate the effective elastic viscosity (this is the equivalent viscosity of
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
           *
           * Stress tensor components @p stress_0_advected represent the stress from the previous
           * timestep $t$ rotated and advected into the current timestep $t+\Delta t_c$.
           * Stress tensor components @p stress_old represent the stress from the previous
           * timestep $t$ advected into the current timestep $t+\Delta t_c$.
           * By the time the viscoelastic strain rate is required to assemble
           * the Stokes system, the stresses have already been rotated and/or advected.
           */
          SymmetricTensor<2,dim>
          calculate_viscoelastic_strain_rate (const SymmetricTensor<2,dim> &strain_rate,
                                              const SymmetricTensor<2, dim> &stress_0_advected,
                                              const SymmetricTensor<2, dim> &stress_old,
                                              const double viscosity_pre_yield,
                                              const double shear_modulus) const;

          /**
           * Compute the elastic time step.
           */
          double
          elastic_timestep () const;

          /**
           * Calculate the ratio between the computational timestep and
           * the elastic timestep.
           */
          double
          calculate_timestep_ratio() const;

        private:
          /**
           * Get the stored stress of the previous timestep. For fields, use a
           * composition evaluator of the old solution. For particles, get the
           * stress directly from the particles, which is available from in.composition
           * by default. However, it can be specified from the input file that the
           * particle property plugin is to use the stress field solution at the
           * particle location instead.
           */
          std::vector<SymmetricTensor<2, dim>>
          retrieve_stress_previous_timestep (const MaterialModel::MaterialModelInputs<dim> &in,
                                             const std::vector<Point<dim>> &quadrature_positions) const;

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
           * Double for fixed elastic time step value, read from parameter file.
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
           * We cache the evaluators that are necessary to evaluate the velocity
           * gradients and the old compositions. They are required to compute the elastic stresses,
           * but are not provided by the material model.
           * By caching the evaluators, we can avoid recreating them every time we need them.
           */
          mutable std::unique_ptr<FEPointEvaluation<dim, dim>> evaluator;
          static constexpr unsigned int n_independent_components = SymmetricTensor<2, dim>::n_independent_components;
          mutable std::unique_ptr<FEPointEvaluation<n_independent_components, dim>> evaluator_composition;

      };
    }
  }
}
#endif
