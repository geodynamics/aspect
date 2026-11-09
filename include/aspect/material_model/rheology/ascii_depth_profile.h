// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2020-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


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
