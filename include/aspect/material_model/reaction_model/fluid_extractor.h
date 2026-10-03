/*
  Copyright (C) 2026 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_reaction_model_fluid_extractor_h
#define _aspect_material_model_reaction_model_fluid_extractor_h

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
      * A simple reaction model which extracts fluid from the system above a
      * user-specified depth. This is intended to be used in conjunction with the
      * reactive fluid transport material model, and is designed to be used with the
      * Tian 2019 solubility model or the Katz 2003 mantle melting model.
      *
      * Above the extraction depth, the extraction rates can be specified to be equal
      * to the current porosity, or to increase with decreasing depth.
      *
      * @ingroup ReactionModel
      */
      template <int dim>
      class FluidExtractor : public ::aspect::SimulatorAccess<dim>
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
           * Compute the rate of extraction of fluid above the user-specified extraction depth,
           * and based on the desired extraction method.
           */
          void
          calculate_reaction_rate_outputs (const typename Interface<dim>::MaterialModelInputs  &in,
                                           typename Interface<dim>::MaterialModelOutputs       &out) const;

          /**
           * Set the parameters for the fluid extractor model. The parameters are parsed from the
           * parent reaction model, and this function is used to pass them to the fluid extractor model.
           */
          void
          set_parameters (const std::string &scheme_name,
                          const double depth,
                          const std::string &method);

        private:
          /**
           * Depth controlling where fluid extraction occurs. Fluid is extracted
           * above the extraction_depth.
           */
          double extraction_depth;

          /**
           * Parameters controlling how the rate of extraction above the extraction depth occurs.
           */
          std::string extraction_method;

          /**
           * Name of the reaction scheme selected by the parent reactive-fluid-transport model.
           */
          std::string reaction_scheme_name;
      };
    }

  }
}

#endif
