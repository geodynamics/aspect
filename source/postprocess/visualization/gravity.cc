// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#include <aspect/postprocess/visualization/gravity.h>
#include <aspect/gravity_model/interface.h>



namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      Gravity<dim>::
      Gravity ()
        :
        DataPostprocessorVector<dim> ("gravity",
                                      update_quadrature_points),
        Interface<dim>("m/s/s")
      {}



      template <int dim>
      void
      Gravity<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double>> &computed_quantities) const
      {
        const unsigned int n_quadrature_points = input_data.evaluation_points.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert (computed_quantities[0].size() == dim,    ExcInternalError());

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            const Tensor<1,dim> g = this->get_gravity_model().gravity_vector (input_data.evaluation_points[q]);
            for (unsigned int k=0; k<dim; ++k)
              computed_quantities[q](k) = g[k];
          }
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(Gravity,
                                                  "gravity",
                                                  "A visualization output object that outputs the gravity vector."
                                                  "\n\n"
                                                  "Physical units: $\\frac{\\text{m}}{\\text{s}^2}$.")
    }
  }
}
