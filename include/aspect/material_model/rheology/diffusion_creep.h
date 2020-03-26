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

#ifndef _aspect_material_model_rheology_diffusion_creep_h
#define _aspect_material_model_rheology_diffusion_creep_h

#include <aspect/global.h>
#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    namespace Rheology
    {
      template <int dim>
      class DiffusionCreep : public ::aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * Constructor.
           */
          DiffusionCreep();

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
          parse_parameters (ParameterHandler &prm, const std::shared_ptr<std::vector<unsigned int>> expected_n_phases_per_composition =
                              std::shared_ptr<std::vector<unsigned int>>());

          /**
           * Compute the viscosity based on the diffusion creep law in the
           * presence of different material phases as determined by the phase function
           * parameters @p gammas. TODO: extend, how many gammas? for which composition?
           */
          double
          compute_viscosity (const double pressure,
                             const double temperature,
                             const unsigned int composition,
                             const std::pair<std::vector<double>, const std::vector<unsigned int>> &gamma_inputs= {}) const;

        private:

          /**
           * List of diffusion creep prefactors A.
           */
          std::vector<double> prefactors_diffusion;

          /**
           * List of diffusion creep grain size exponenents m.
           */
          std::vector<double> grain_size_exponents_diffusion;

          /**
           * List of diffusion creep activation energies E.
           */
          std::vector<double> activation_energies_diffusion;

          /**
           * List of diffusion creep activation volumes V.
           */
          std::vector<double> activation_volumes_diffusion;

          /**
           * Diffusion creep grain size d.
           */
          std::vector<double> grain_size;
      };
    }
  }
}
#endif
