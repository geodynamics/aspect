/*
  Copyright (C) 2020 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_rheology_visco_plastic_h
#define _aspect_material_model_rheology_visco_plastic_h

#include <aspect/global.h>
#include <aspect/material_model/interface.h>
#include <aspect/material_model/utilities.h>
#include <aspect/material_model/rheology/strain_dependent.h>
#include <aspect/material_model/rheology/diffusion_creep.h>
#include <aspect/material_model/rheology/dislocation_creep.h>
#include <aspect/material_model/rheology/frank_kamenetskii.h>
#include <aspect/material_model/rheology/peierls_creep.h>
#include <aspect/material_model/rheology/constant_viscosity_prefactors.h>
#include <aspect/material_model/rheology/drucker_prager.h>
#include <aspect/material_model/rheology/elasticity.h>
#include <aspect/simulator_access.h>

#include<deal.II/fe/component_mask.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    namespace Rheology
    {

      template <int dim>
      class ViscoPlastic : public ::aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * Constructor.
           */
          ViscoPlastic();

          /**
           * This function calculates viscosities assuming that all the compositional fields
           * experience the same strain rate (isostrain).
           */
          std::pair<std::vector<double>, std::vector<bool> >
          calculate_isostrain_viscosities ( const MaterialModel::MaterialModelInputs<dim> &in,
                                            const unsigned int i,
                                            const std::vector<double> &volume_fractions,
                                            const std::vector<double> &phase_function_values = std::vector<double>(),
                                            const std::vector<unsigned int> &n_phases_per_composition =
                                              std::vector<unsigned int>()) const;

          /**
           * A function that fills the viscosity derivatives in the
           * MaterialModelOutputs object that is handed over, if they exist.
           * Does nothing otherwise.
           */
          void compute_viscosity_derivatives(const unsigned int point_index,
                                             const std::vector<double> &volume_fractions,
                                             const std::vector<double> &composition_viscosities,
                                             const MaterialModel::MaterialModelInputs<dim> &in,
                                             MaterialModel::MaterialModelOutputs<dim> &out,
                                             const std::vector<double> &phase_function_values = std::vector<double>(),
                                             const std::vector<unsigned int> &n_phases_per_composition =
                                               std::vector<unsigned int>()) const;


          /**
           * A function that returns a ComponentMask that represents all compositional
           * fields that should be considered 'volumetric', that is representing a
           * physical proportion of the material, e.g. volume fraction of peridotite
           * (as opposed to non-volumetric quantities like the amount of finite-strain).
           */
          ComponentMask get_volumetric_composition_mask() const;

          /**
           * Declare the parameters this function takes through input files.
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Read the parameters this class declares from the parameter file.
           * If @p expected_n_phases_per_composition points to a vector of
           * unsigned integers this is considered the number of phase transitions
           * for each compositional field and will be checked against the parsed
           * parameters.
           */
          void
          parse_parameters (ParameterHandler &prm,
                            const std::shared_ptr<std::vector<unsigned int>> &expected_n_phases_per_composition =
                              std::shared_ptr<std::vector<unsigned int>>());



          double ref_visc;
          double min_strain_rate;

          /**
           * Enumeration for selecting which viscosity averaging scheme to use.
           */
          MaterialUtilities::CompositionalAveragingOperation viscosity_averaging;


          Rheology::StrainDependent<dim> strain_rheology;
          /*
           * Input parameters for the drucker prager plasticity.
           */
          Rheology::DruckerPragerParameters drucker_prager_parameters;
          /**
           * Object for computing viscoelastic viscosities and stresses.
           */
          Rheology::Elasticity<dim> elastic_rheology;

          /**
           * Whether to include viscoelasticity in the constitutive formulation.
           */
          bool use_elasticity;


        private:

          double ref_strain_rate;
          double min_visc;
          double max_visc;

          /**
           * Enumeration for selecting which type of viscous flow law to use.
           * Select between diffusion, dislocation or composite.
           */
          enum ViscosityScheme
          {
            diffusion,
            dislocation,
            frank_kamenetskii,
            composite
          } viscous_flow_law;

          /**
           * Enumeration for selecting which type of yield mechanism to use.
           * Select between Drucker Prager and stress limiter.
           */
          enum YieldScheme
          {
            stress_limiter,
            drucker_prager
          } yield_mechanism;

          /**
           * Whether to allow negative pressures to be used in the computation
           * of plastic yield stresses and viscosities. If false, the minimum
           * pressure in the plasticity formulation will be set to zero.
           */
          bool allow_negative_pressures_in_plasticity;
          std::vector<double> exponents_stress_limiter;

          /**
           * temperature gradient added to temperature used in the flow law.
           */
          double adiabatic_temperature_gradient_for_viscosity;

          /**
           * Objects for computing viscous creep viscosities.
           */
          Rheology::DiffusionCreep<dim> diffusion_creep;
          Rheology::DislocationCreep<dim> dislocation_creep;
          std::unique_ptr<Rheology::FrankKamenetskii<dim> > frank_kamenetskii_rheology;

          /**
           * Whether to include peierls creep in the constitutive formulation.
           */
          bool use_peierls_creep;

          /**
           * Objects for computing peierls creep viscosities.
           */
          std::unique_ptr<Rheology::PeierlsCreep<dim> > peierls_creep;

          /**
           * Object for computing the viscosity multiplied by a constant prefactor.
           * This multiplication step is done just prior to calculating the effective
           * viscoelastic viscosity or plastic viscosity.
           */
          Rheology::ConstantViscosityPrefactors<dim> constant_viscosity_prefactors;

          /*
           * Objects for computing plastic stresses, viscosities, and additional outputs
           */
          Rheology::DruckerPrager<dim> drucker_prager_plasticity;

      };
    }
  }
}
#endif
