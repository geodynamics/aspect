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


#include <aspect/postprocess/visualization/prescribed_dilation.h>
#include <aspect/prescribed_dilation/interface.h>



namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      PrescribedDilation<dim>::
      PrescribedDilation ()
        :
        DataPostprocessorScalar<dim> ("prescribed_dilation",
                                      update_quadrature_points),
        Interface<dim>("1/s")
      {}



      template <int dim>
      void
      PrescribedDilation<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double>> &computed_quantities) const
      {
        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert (computed_quantities[0].size() == 1,                   ExcInternalError());
        Assert (input_data.solution_values[0].size() == this->introspection().n_components,           ExcInternalError());

        const aspect::PrescribedDilation::Manager<dim> &prescribed_dilation = this->get_prescribed_dilation_manager();
        MaterialModel::MaterialModelInputs<dim> in(input_data,
                                                   this->introspection());
        aspect::PrescribedDilation::PrescribedDilationOutputs<dim> out(in.n_evaluation_points());

        prescribed_dilation.evaluate (in, out);

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          computed_quantities[q](0) = out.dilation[q];
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(PrescribedDilation,
                                                  "prescribed dilation",
                                                  "A visualization output postprocessor that outputs prescribed dilation. "
                                                  "Physical units: $\\si{\\per\\seconds}$.")
    }
  }
}
