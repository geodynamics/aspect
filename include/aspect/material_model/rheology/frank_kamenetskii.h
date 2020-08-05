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

#ifndef _aspect_material_model_rheology_frank_kamenetskii_h
#define _aspect_material_model_rheology_frank_kamenetskii_h

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
      /**
       * A class that computes a Frank-Kamenetskii viscosity approximation
       * of the form:
       * viscosity = A * exp(E * 0.5 * (1.0-(T/ref_T)))
       * A: prefactor of viscosity, E: adjusted viscosity ratio,
       * ref_T: reference temperature, T: temperature.
       */

      template <int dim>
      class FrankKamenetskii : public ::aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * Constructor.
           */
          FrankKamenetskii();

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
           * Compute the viscosity based on the Frank-Kamenetskii approximation.
           */
          double
          compute_viscosity (const double temperature,
                             const unsigned int composition) const;

        private:
          /**
           * List of Frank-Kamenetskii viscosity ratios (E).
           */
          std::vector<double> viscosity_ratios_frank_kamenetskii;

          /**
           * List of Frank-Kamenetskii prefactors (A).
           */
          std::vector<double> prefactors_frank_kamenetskii;
      };
    }
  }
}
#endif
