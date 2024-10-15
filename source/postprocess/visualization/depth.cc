/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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


#include <aspect/postprocess/visualization/depth.h>
#include <aspect/geometry_model/interface.h>



namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      Depth<dim>::
      Depth ()
        :
        DataPostprocessorScalar<dim> ("depth",
                                      update_quadrature_points),
        Interface<dim>("m")
      {}



      template <int dim>
      void
      Depth<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double>> &computed_quantities) const
      {
        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert (computed_quantities[0].size() == 1,                   ExcInternalError());
        Assert (input_data.solution_values[0].size() == this->introspection().n_components,           ExcInternalError());

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            computed_quantities[q](0) = this->get_geometry_model().depth (input_data.evaluation_points[q]);
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(Depth,
                                                  "depth",
                                                  "A visualization output postprocessor that outputs "
                                                  "the depth for all points inside the domain, as "
                                                  "determined by the geometry model."
                                                  "\n\n"
                                                  "It is worth comparing this visualization postprocessor with the "
                                                  "one called ``surface elevation''. The current one is used to visualize a volume "
                                                  "variable, whereas the latter only outputs information on "
                                                  "the surface. Moreover ``depth'' is based on a member function of "
                                                  "the geometry models that is documented as never returning a "
                                                  "number less than zero -- in other words, it returns the depth "
                                                  "of an evaluation point with regard to a reference surface that "
                                                  "defines a zero depth, but for points that lie above this "
                                                  "reference surface, it returns zero. As a consequence, it cannot "
                                                  "be used to visualize positive elevations, whereas the the one "
                                                  "called ``surface elevation'' can."
                                                  "\n\n"
                                                  "Physical units: \\si{\\meter}.")
    }
  }
}
