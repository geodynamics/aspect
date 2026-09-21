// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_postprocess_visualization_compositional_vector_h
#define _aspect_postprocess_visualization_compositional_vector_h

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
       * A class derived from DataPostprocessor that can output sets of dim
       * compositional fields as vectors.
       */
      template <int dim>
      class CompositionalVector
        : public DataPostprocessor<dim>,
          public SimulatorAccess<dim>,
          public Interface<dim>
      {
        public:
          CompositionalVector ();

          std::vector<std::string>
          get_names () const override;

          std::vector<DataComponentInterpretation::DataComponentInterpretation>
          get_data_component_interpretation () const override;

          UpdateFlags
          get_needed_update_flags () const override;

          void
          evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                                std::vector<Vector<double>> &computed_quantities) const override;

          static
          void
          declare_parameters (ParameterHandler &prm);

          void
          parse_parameters (ParameterHandler &prm) override;

        private:
          /**
           * Names of vector fields.
           *
           * This variable is read from the parameter file through a parameter called 'Names of vectors'.
           */
          std::vector<std::string> vector_names;

          /**
           * Stores sets of compositional field indices to
           * be visualized together as vector fields.
           */
          std::vector<std::vector<unsigned int>> sets;
      };
    }
  }
}

#endif
