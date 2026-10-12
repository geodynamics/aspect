// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2016-2022 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_postprocess_visualization_maximum_horizontal_compressive_stress_h
#define _aspect_postprocess_visualization_maximum_horizontal_compressive_stress_h

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
       * A class that computes a field of horizontal vectors that
       * represent the direction of maximal horizontal compressive
       * stress. For an exact definition, see the documentation of
       * this plugin in the manual.
       *
       * The member functions are all implementations of those declared in the
       * base class. See there for their meaning.
       */
      template <int dim>
      class MaximumHorizontalCompressiveStress
        : public DataPostprocessor<dim>,
          public SimulatorAccess<dim>,
          public Interface<dim>
      {
        public:
          /**
           * Constructor.
           */
          MaximumHorizontalCompressiveStress();

          void
          evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                                std::vector<Vector<double>> &computed_quantities) const override;

          /**
           * Return the vector of strings describing the names of the computed
           * quantities. Given the purpose of this class, this is a vector
           * with entries all equal to the name of the plugin.
           */
          std::vector<std::string> get_names () const override;

          /**
           * This functions returns information about how the individual
           * components of output files that consist of more than one data set
           * are to be interpreted. The returned value is
           * DataComponentInterpretation::component_is_scalar repeated
           * SymmetricTensor::n_independent_components times. (These
           * components should really be part of a symmetric tensor, but
           * deal.II does not allow marking components as such.)
           */
          std::vector<DataComponentInterpretation::DataComponentInterpretation>
          get_data_component_interpretation () const override;

          /**
           * Return which data has to be provided to compute the derived
           * quantities. The flags returned here are the ones passed to the
           * constructor of this class.
           */
          UpdateFlags get_needed_update_flags () const override;
      };
    }
  }
}

#endif
