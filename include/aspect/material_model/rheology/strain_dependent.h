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

#ifndef _aspect_material_model_rheology_strain_dependent_h
#define _aspect_material_model_rheology_strain_dependent_h

#include <aspect/global.h>
#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>

#include<deal.II/fe/component_mask.h>
#include <deal.II/matrix_free/fe_point_evaluation.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace Rheology
    {
      /**
       * Enumeration for selecting which type of weakening mechanism to use.
       * For none, no strain weakening occurs.
       * Otherwise, the material can be weakened based on the second
       * invariant of the full finite strain tensor, the total accumulated
       * strain, or the plastic strain and viscous strain can be tracked
       * separately and used only for the corresponding (plastic or viscous)
       * part of the viscosity computation.
       */
      enum WeakeningMechanism
      {
        none,
        finite_strain_tensor,
        total_strain,
        plastic_weakening_with_plastic_strain_only,
        plastic_weakening_with_total_strain_only,
        plastic_weakening_with_plastic_strain_and_viscous_weakening_with_viscous_strain,
        viscous_weakening_with_viscous_strain_only
      };

      /**
       * Enumeration for selecting which type of healing mechanism to use.
       * For the case no healing, no strain healing occurs.
       * Otherwise, the strain is healed (reduced) as a function of temperature,
       * a user defined healing time scale and additional parameters.
       * Future models could consider strain healing formulations that are a function
       * of time, deformation rate, or other parameters.
       */
      enum HealingMechanism
      {
        no_healing,
        temperature_dependent
      };

      template <int dim>
      class StrainDependent : public ::aspect::SimulatorAccess<dim>
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
           * A function that computes by how much the rheologic parameters change
           * if strain weakening is applied. Given a vector @p composition of all
           * fields, it returns reduction factors for the cohesion, friction angle
           * and the prefactor of the viscous flow law(s) for the compositional
           * field with index @p j. The reason all fields are passed is that
           * the weakening factors can depend on the values of fields that track
           * different measures of previously applied strain.
           */
          std::array<double, 3>
          compute_strain_weakening_factors(const std::vector<double> &composition,
                                           const unsigned int j) const;

          /**
           * @deprecated: Deprecated version of the function of the same
           * name described above.
           */
          DEAL_II_DEPRECATED
          std::array<double, 3>
          compute_strain_weakening_factors(const unsigned int j,
                                           const std::vector<double> &composition) const;

          /**
           * A function that alters the viscous weakening factor based on the
           * temperature field.
           */
          std::array<double, 3>
          apply_temperature_dependence_to_strain_weakening_factors(const std::array<double, 3> &weakening_factors,
                                                                   const double temperature,
                                                                   const unsigned int j) const;

          /**
           * A function that computes the strain healing (reduction in accumulated strain)
           */
          double
          calculate_strain_healing (const MaterialModel::MaterialModelInputs<dim> &in,
                                    const unsigned int j) const;

          /**
           * A function that computes by how much the cohesion and internal friction
           * angle for a given compositional field are weakened under the influence
           * of a given strain.
           */
          std::pair<double, double>
          calculate_plastic_weakening (const double strain_ii,
                                       const unsigned int j) const;

          /**
           * A function that computes by how much the diffusion and dislocation
           * prefactors for a given compositional field are weakened under the
           * influence of a given strain.
           */
          double
          calculate_viscous_weakening (const double strain_ii,
                                       const unsigned int j) const;

          /**
           * Whether to use the temperature-activated viscous strain weakening.
           */
          bool use_temperature_activated_strain_softening;

          /**
           * A function that fills the reaction terms for the finite strain tensor in
           * MaterialModelOutputs object that is handed over. It assumes the first
           * component of the finite strain tensor is named 's11' and all other
           * components follow this compositional field.
           */
          void compute_finite_strain_reaction_terms (const MaterialModel::MaterialModelInputs<dim> &in,
                                                     MaterialModel::MaterialModelOutputs<dim> &out) const;
          /**
           * A function that fills the reaction terms for the finite strain invariant(s) in
           * MaterialModelOutputs object that is handed over.
           */
          void fill_reaction_outputs (const MaterialModel::MaterialModelInputs<dim> &in,
                                      const int i,
                                      const double min_strain_rate,
                                      const bool plastic_yielding,
                                      MaterialModel::MaterialModelOutputs<dim> &out) const;

          /**
           * A function that returns a ComponentMask, which indicates that components
           * associated with strain should be excluded during the volume fraction computation.
           */
          ComponentMask get_strain_composition_mask() const;

          /**
           * A function that returns the selected type of strain weakening mechanism.
           */
          WeakeningMechanism
          get_weakening_mechanism () const;

          /**
           * A function that returns the selected type of strain healing mechanism.
           */
          HealingMechanism
          get_healing_mechanism () const;

        private:

          WeakeningMechanism weakening_mechanism;

          HealingMechanism healing_mechanism;

          /**
           * The start of the strain interval (plastic or total strain)
           * within which cohesion and angle of friction should be weakened.
           */
          std::vector<double> start_plastic_strain_weakening_intervals;

          /**
           * The end of the strain interval (plastic or total strain)
           * within which cohesion and angle of friction should be weakened.
           */
          std::vector<double> end_plastic_strain_weakening_intervals;

          /**
           * The factor specifying the amount of weakening of the
           * cohesion over the prescribed strain interval (plastic or total strain).
           */
          std::vector<double> cohesion_strain_weakening_factors;

          /**
           * The factor specifying the amount of weakening of the
           * internal friction angles over the prescribed strain interval
           * (plastic or total strain).
           */
          std::vector<double> friction_strain_weakening_factors;

          /**
           * The start of the strain interval (viscous or total strain)
           * within which cohesion and angle of friction should be weakened.
           */
          std::vector<double> start_viscous_strain_weakening_intervals;

          /**
           * The end of the strain interval (viscous or total strain)
           * within which cohesion and angle of friction should be weakened.
           */
          std::vector<double> end_viscous_strain_weakening_intervals;

          /**
           * The factor specifying the amount of weakening over
           * the prescribed strain interval (viscous or total strain).
           */
          std::vector<double> viscous_strain_weakening_factors;

          /**
           * The four temperatures that parameterize the temperature-activated strain softening.
           * These can be different for each compositional field.
           * ------            -------- 1
           *       \          /
           *        \        /
           *         \______/ _ _ _ _ _ viscous_weakening_factor[2]
           *
           * ----------------------------> T
           *     T0  T1   T2  T3
           */
          std::vector<double> viscous_strain_weakening_T0;
          std::vector<double> viscous_strain_weakening_T1;
          std::vector<double> viscous_strain_weakening_T2;
          std::vector<double> viscous_strain_weakening_T3;

          /**
           * The healing rate used in the temperature dependent strain healing model.
           */
          double strain_healing_temperature_dependent_recovery_rate;

          /**
           * A prefactor of viscosity used in the strain healing calculation.
           */
          double strain_healing_temperature_dependent_prefactor;

          /**
           * We cache the evaluators that are necessary to evaluate the velocity
           * gradients and compositions.
           * By caching the evaluator, we can avoid recreating them
           * every time we need it.
           */
          mutable std::unique_ptr<FEPointEvaluation<dim, dim>> evaluator;
          mutable std::vector<std::unique_ptr<FEPointEvaluation<1, dim>>> composition_evaluators;
      };
    }
  }
}
#endif
