/*
  Copyright (C) 2019 by the authors of the ASPECT code.

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

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

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
           * if strain weakening is applied. Given a compositional field with
           * the index j and a vector of all compositional fields, it returns
           * reduction factors for the cohesion, friction angle and the prefactor
           * of the viscous flow law(s) used in the computation for that composition.
           */
          std::array<double, 3>
          compute_strain_weakening_factors(const unsigned int j,
                                           const std::vector<double> &composition) const;

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
           * The healing rate used in the temperature dependent strain healing model.
           */
          double strain_healing_temperature_dependent_recovery_rate;

          /**
           * A prefactor of viscosity used in the strain healing calculation.
           */
          double strain_healing_temperature_dependent_prefactor;
      };
    }
  }
}
#endif


