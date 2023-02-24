/*
  Copyright (C) 2011 - 2022 by the authors of the ASPECT code.

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


#include <aspect/postprocess/visualization/depth_including_mesh_deformation.h>
#include <aspect/geometry_model/interface.h>



namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      DepthIncMeshDef<dim>::
      DepthIncMeshDef ()
        :
        DataPostprocessorScalar<dim> ("depth",
                                      update_quadrature_points),
        Interface<dim>("m")
      {}



      template <int dim>
      void
      DepthIncMeshDef<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double>> &computed_quantities) const
      {
        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert (computed_quantities[0].size() == 1,                   ExcInternalError());
        Assert (input_data.solution_values[0].size() == this->introspection().n_components,           ExcInternalError());

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            computed_quantities[q](0) = this->get_geometry_model().depth_including_mesh_deformation (input_data.evaluation_points[q]);
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(DepthIncMeshDef,
                                                  "depth including mesh deformation",
                                                  "A visualization output postprocessor that outputs "
                                                  "the depth for all points inside the domain, as "
                                                  "determined by the current model surface. This plugin "
                                                  "will include changes to the surface from mesh deformation."
                                                  "\n\n"
                                                  "Physical units: \\si{\\meter}.")
    }
  }
}
