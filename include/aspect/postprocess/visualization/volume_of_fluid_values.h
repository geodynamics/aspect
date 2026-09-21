// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2016-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_postprocess_visualization_volume_of_fluid_values_h
#define _aspect_postprocess_visualization_volume_of_fluid_values_h

#include <aspect/postprocess/visualization.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      /**
       * Plugin to output volume fractions and reconstructed interface data
       * from the Volume of Fluid Interface Tracking data.
       */
      template <int dim>
      class VolumeOfFluidValues
        : public DataPostprocessor<dim>,
          public SimulatorAccess<dim>,
          public Interface<dim>
      {
        public:
          /**
           * Standard constructor
           */
          VolumeOfFluidValues ();

          /**
           * Get the list of names for the components that will be produced by
           * this postprocessor
           */
          std::vector<std::string>
          get_names () const override;

          /**
           * Get the list of component interpretations for the components that
           * will be produced by this postprocessor
           */
          std::vector<DataComponentInterpretation::DataComponentInterpretation>
          get_data_component_interpretation () const override;

          /**
           * Get required update flags
           */
          UpdateFlags
          get_needed_update_flags () const override;

          /**
           * Produce that data based on provided solution data
           */
          void
          evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                                std::vector<Vector<double>> &computed_quantities) const override;

          /**
           * Declare the parameters this class takes through input files.
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Read the parameters this class declares from the parameter file.
           */
          void
          parse_parameters (ParameterHandler &prm) override;

        private:
          /**
           * Stored list of names for the produced components
           */
          std::vector<std::string> volume_of_fluid_names;

          /**
           * Stored list of interpretations for produced components
           */
          std::vector<DataComponentInterpretation::DataComponentInterpretation> interpretations;

          /**
           * If true, the data output will include a field that has the
           * reconstructed fluid interface as the zero contour
           *
           * This variable is read from the parameter file through a parameter called 'Output interface reconstruction contour'.
           */
          bool include_contour;

          /**
           * If true, the data output will include the normal vector for the reconstructed interface
           *
           * This variable is read from the parameter file through a parameter called 'Output interface normals'.
           */
          bool include_normal;
      };
    }
  }
}

#endif
