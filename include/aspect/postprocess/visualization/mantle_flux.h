/*
  Copyright (C) 2011 - 2026 by the authors of the ASPECT code.

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

#ifndef _aspect_postprocess_visualization_mantle_flux_h
#define _aspect_postprocess_visualization_mantle_flux_h

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
       * A visualization postprocessor that writes the local values used by
       * the mantle flux statistics postprocessor. The structure field is 1
       * in plumes, -1 in slabs, and 0 in the background mantle.
       *
       * @ingroup Postprocessing
       * @ingroup Visualization
       */
      template <int dim>
      class MantleFlux : public DataPostprocessor<dim>,
        public SimulatorAccess<dim>,
        public Interface<dim>
      {
        public:
          /**
           * Constructor.
           */
          MantleFlux ();

          /**
           * Return the names of the fields written to graphical output.
           */
          std::vector<std::string>
          get_names () const override;

          /**
           * Mark every output field as a scalar.
           */
          std::vector<DataComponentInterpretation::DataComponentInterpretation>
          get_data_component_interpretation () const override;

          /**
           * Return units that follow the global choice to use seconds or
           * years in output.
           */
          std::string
          get_physical_units () const override;

          /**
           * Request the solution information needed by the material model.
           */
          UpdateFlags
          get_needed_update_flags () const override;

          /**
           * Compute plume, slab, and local flux values for graphical output.
           */
          void
          evaluate_vector_field (const DataPostprocessorInputs::Vector<dim> &input_data,
                                 std::vector<Vector<double>> &computed_quantities) const override;

          /**
           * Ensure that mantle flux statistics is active and available.
           */
          std::list<std::string>
          required_other_postprocessors () const override;
      };
    }
  }
}

#endif
