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

#ifndef _aspect_material_model_reaction_model_crust_and_lithosphere_formation_h
#define _aspect_material_model_reaction_model_crust_and_lithosphere_formation_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    namespace ReactionModel
    {

      /**
      * A simplified model to calculate the change in composition upon melting
      * of average mantle as it approaches the surface to produce a basaltic crust
      * and a harzburgitic lithosphere. The model assumes that the crust is
      * generated at a constant depth, and that the lithosphere is generated
      * below the crust at a constant depth. The reaction producing crust and
      * lithosphere only occurs in material that is upwelling, but does not take
      * into account the temperature of the upwelling material.
      *
      * @ingroup ReactionModel
      */
      template <int dim>
      class CrustLithosphereFormation : public ::aspect::SimulatorAccess<dim>
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
           * Compute the change in composition for the basalt and harzburgite chemical
           * fields upon melting as mantle material reaches the surface. We assume that
           * all upwelling material is converted to basalt or harzburgite as it reaches
           * the crustal and lithospheric depths, respectively. The reaction terms are
           * computed for as many points as are provided in @p in and they are stored
           * in the material model outputs object @p out.
           */
          void
          calculate_reaction_terms (const typename Interface<dim>::MaterialModelInputs  &in,
                                    typename Interface<dim>::MaterialModelOutputs       &out) const;

        private:
          /**
           * Parameters controlling where the generation of crust and lithosphere
           * occurs. Crust is generated above the crustal_thickness, and lithosphere
           * is generated below the crustal_thickness and down to a depth that is the
           * sum of crustal_thickness and lithosphere_thickness.
           */
          double crustal_thickness;
          double lithosphere_thickness;

          /**
           * The indices of the compositional fields that store the basalt and
           * harzburgite chemical compositions.
           */
          unsigned int basalt_index;
          unsigned int harzburgite_index;
      };
    }

  }
}

#endif
