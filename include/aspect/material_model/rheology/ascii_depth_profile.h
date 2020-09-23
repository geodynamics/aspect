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


#ifndef _aspect_material_model_rheology_ascii_depth_profile_h
#define _aspect_material_model_rheology_ascii_depth_profile_h

#include <aspect/material_model/interface.h>
#include <aspect/utilities.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/point.h>


namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    namespace Rheology
    {

      /**
       * A class that implements a prescribed viscosity with depth determined from
       * an AsciiDataProfile input file.
       *
       * @ingroup Rheology
       */
      template <int dim>
      class AsciiDepthProfile : public Utilities::AsciiDataProfile<dim> , public SimulatorAccess<dim>
      {
        public:
          /**
           * Constructor. Initialize variables.
           */
          AsciiDepthProfile ();

          /**
           * Initialization function.
           */
          void initialize ();

          // avoid -Woverloaded-virtual:
          using Utilities::AsciiDataProfile<dim>::initialize;

          /**
          * Return the viscosity at a given point of the domain.
          */
          double compute_viscosity (const double depth) const;

          /**
           * Declare the parameters for the input files.
           */
          static
          void
          declare_parameters (ParameterHandler  &prm,
                              const std::string &subsection_name = "Ascii data model");

          /**
           * Read the parameters from the parameter file.
           */
          void
          parse_parameters (ParameterHandler &prm,
                            const std::string &subsection_name = "Ascii data model");

        private:

          /**
           * The column index of the viscosity in the data file.
           */
          unsigned int viscosity_index;
      };
    }
  }
}


#endif