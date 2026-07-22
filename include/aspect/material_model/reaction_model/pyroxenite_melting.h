/*
  Copyright (C) 2024 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_reaction_model_pyroxenite_melting_h
#define _aspect_material_model_reaction_model_pyroxenite_melting_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace ReactionModel
    {

      /**
       * A melt model that calculates melt fraction and entropy change
       * according to the melting model for dry peridotite of Katz, 2003.
       * This also includes a computation of the latent heat of melting (if the latent heat
       * heating model is active).
       *
       * These functions can be used in the calculation of melting and melt transport
       * in the melt_simple material model and can be extended to other material models
       *
       * @ingroup ReactionModel
       */
      template <int dim>
      class PyroxeniteMelting : public ::aspect::SimulatorAccess<dim>
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
           * Percentage of material that is molten for a given @p temperature and
           * @p pressure (assuming equilibrium conditions). Melting model after Katz,
           * 2003, for dry peridotite.
           */
          double
          melt_fraction (const double temperature,
                         const double pressure) const;

        private:
          /**
           * Parameters for melting of pyroxenite after Sobolev et al., 2011
           */

          // for the melting temperature
          double D1;    // °C
          double D2;  // °C/Pa
          double D3; // °C/(Pa^2)

          // for the melt-fraction dependence of productivity
          double E1;
          double E2;
      };
    }

  }
}

#endif
