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

#ifndef _aspect_material_model_rheology_composite_visco_plastic_h
#define _aspect_material_model_rheology_composite_visco_plastic_h

#include <aspect/global.h>
#include <aspect/material_model/interface.h>
#include <aspect/material_model/utilities.h>
#include <aspect/material_model/rheology/diffusion_creep.h>
#include <aspect/material_model/rheology/dislocation_creep.h>
#include <aspect/material_model/rheology/peierls_creep.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    namespace Rheology
    {

      template <int dim>
      class CompositeViscoPlastic : public ::aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * Constructor.
           */
          CompositeViscoPlastic();

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

          /**
           * Compute the viscosity based on the composite viscous creep law.
           * If @p expected_n_phases_per_composition points to a vector of
           * unsigned integers this is considered the number of phase transitions
           * for each compositional field and viscosity will be first computed on
           * each phase and then averaged for each compositional field.
           */
          double
          compute_viscosity (const double pressure,
                             const double temperature,
                             const unsigned int composition,
                             const SymmetricTensor<2,dim> &strain_rate,
                             std::vector<double> &partial_strain_rates,
                             const std::vector<double> &phase_function_values = std::vector<double>(),
                             const std::vector<unsigned int> &n_phases_per_composition = std::vector<unsigned int>()) const;

          /**
            * Compute the strain rate and first stress derivative
            * as a function of stress based on the composite viscous creep law.
            * If @p expected_n_phases_per_composition points to a vector of
            * unsigned integers this is considered the number of phase transitions
            * for each compositional field and viscosity will be first computed on
            * each phase and then averaged for each compositional field.
            */
          std::pair<double, double>
          compute_strain_rate_and_derivative (const double creep_stress,
                                              const double pressure,
                                              const double temperature,
                                              const DiffusionCreepParameters diffusion_creep_parameters,
                                              const DislocationCreepParameters dislocation_creep_parameters,
                                              const PeierlsCreepParameters peierls_creep_parameters,
                                              const double min_viscosity,
                                              const double max_viscosity) const;

        private:

          /**
           * Whether to use different creep mechanisms
           */
          bool use_diffusion_creep;
          bool use_dislocation_creep;
          bool use_peierls_creep;

          /**
           * Pointers to objects for computing viscous creep viscosities.
           */
          std::unique_ptr<Rheology::DiffusionCreep<dim>> diffusion_creep;
          std::unique_ptr<Rheology::DislocationCreep<dim>> dislocation_creep;
          std::unique_ptr<Rheology::PeierlsCreep<dim>> peierls_creep;

          std::vector<double> min_viscosities;
          std::vector<double> max_viscosities;

          double min_strain_rate;
          double strain_rate_residual_threshold;
          unsigned int stress_max_iteration_number;
      };
    }
  }
}
#endif
