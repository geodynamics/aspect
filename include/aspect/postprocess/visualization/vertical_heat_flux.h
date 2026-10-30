// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2021 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_postprocess_visualization_vertical_heat_flux_h
#define _aspect_postprocess_visualization_vertical_heat_flux_h

#include <aspect/postprocess/visualization.h>
#include <aspect/simulator_access.h>

#include <deal.II/numerics/data_postprocessor.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      /**
       * A class derived from DataPostprocessor that takes an output vector
       * and computes a variable that represents the density at every point.
       *
       * The member functions are all implementations of those declared in the
       * base class. See there for their meaning.
       */
      template <int dim>
      class VerticalHeatFlux
        : public DataPostprocessorScalar<dim>,
          public SimulatorAccess<dim>,
          public Interface<dim>
      {
        public:
          /**
           * Constructor
           */
          VerticalHeatFlux ();

          /**
           * Evaluate the vertical heat flux for the current cell.
           */
          void
          evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                                std::vector<Vector<double>> &computed_quantities) const override;
      };
    }
  }
}

#endif
