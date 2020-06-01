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

#ifndef _aspect_material_model_rheology_peierls_creep_h
#define _aspect_material_model_rheology_peierls_creep_h

#include <aspect/global.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    namespace Rheology
    {
      template <int dim>
      class PeierlsCreep : public ::aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * Constructor.
           */
          PeierlsCreep();

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
           * Compute the viscosity based on the selected Peierls creep flow law.
           */
          double
          compute_viscosity (const double strain_rate,
                             const double temperature,
                             const unsigned int composition) const;

        private:
          /**
           * Enumeration for selecting which type of Peierls creep flow law to use.
           * Currently, the only available option directly calculates the viscosity
           * using an approximation to the model of Mei et al., 2010. The full
           * derivation for this approximation (derived by Magali Billen) can be
           * found at https://ucdavis.app.box.com/s/scaorcblr9u294836pgk4hyf60gv7psr.
           */
          enum PeierlsCreepScheme
          {
            mei_viscosity_approx
          } peierls_creep_flow_law;

          /**
           * List of Peirls creep prefactors A.
           */
          std::vector<double> prefactors_peierls;

          /**
           * List of Peierls creep activation energies E.
           */
          std::vector<double> activation_energies_peierls;

          /**
           * Peierls stress
           */
          double peierls_stress;

          /**
           * Peierls fitting parameter
           */
          double peierls_fitting_parameter;

          /**
           * Peierls fitting exponent
           */
          double peierls_fitting_exponent;

      };
    }
  }
}
#endif
