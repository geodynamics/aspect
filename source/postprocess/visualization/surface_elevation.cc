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
#include <aspect/simulator.h>
#include <aspect/postprocess/visualization/surface_elevation.h>

namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      SurfaceElevation<dim>::
      SurfaceElevation ()
        :
        DataPostprocessorScalar<dim> ("surface_elevation",
                                      update_quadrature_points),
        Interface<dim>("m")
      {}



      template <int dim>
      void
      SurfaceElevation<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double>> &computed_quantities) const
      {
        // Initialize everything to zero, so that we can ignore faces we are
        // not interested in (namely, those not labeled as 'top' or 'bottom')
        for (auto &quantity : computed_quantities)
          quantity(0) = 0;

        // We only want to output dynamic topography at the top
        // boundary. We know that this class will only be called on faces
        // (because it is derived from SurfaceOnlyVisualization), so it is safe
        // to query both cell and face associated with the DataPostprocessorInputs
        // object.
        const types::boundary_id boundary_id
          = input_data.template get_cell<dim>()->face(input_data.get_face_number())->boundary_id();
        if (this->get_geometry_model().translate_id_to_symbol_name (boundary_id) == "top")
          for (unsigned int q=0; q<input_data.evaluation_points.size(); ++q)
            computed_quantities[q](0)
              = this->get_geometry_model().height_above_reference_surface(input_data.evaluation_points[q]);
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(SurfaceElevation,
                                                  "surface elevation",
                                                  "This postprocessor is used to visualize the elevation of points "
                                                  "on the surface of the geometry relative to a ``reference "
                                                  "elevation'' defined by an undeformed geometry. It can be used, "
                                                  "for example, to visualize an initial topography field, as well "
                                                  "as the result of dynamic surface deformation due to a free "
                                                  "surface."
                                                  "\n\n"
                                                  "The surface elevation is computed only for those parts of the "
                                                  "boundary that the geometry description marks as ``top''. "
                                                  "On all other parts of the boundary, the class outputs a zero "
                                                  "elevation."
                                                  "\n\n"
                                                  "It is worth comparing this visualization postprocessor with the "
                                                  "one called ``depth''. The latter is used to visualize a volume "
                                                  "variable, whereas the current one only outputs information on "
                                                  "the surface. Moreover ``depth'' is based on a member function of "
                                                  "the geometry models that is documented as never returning a "
                                                  "number less than zero -- in other words, it returns the depth "
                                                  "of an evaluation point with regard to a reference surface that "
                                                  "defines a zero depth, but for points that lie above this "
                                                  "reference surface, it returns zero. As a consequence, it cannot "
                                                  "be used to visualize positive elevations, whereas the current "
                                                  "visualization postprocessor can."
                                                  "\n\n"
                                                  "Finally, it is worth pointing out the ``topography'' "
                                                  "postprocessor (not a visualization postprocessor) that returns "
                                                  "the surface elevation as a point cloud into a text file. The "
                                                  "information is comparable to what the current object creates, "
                                                  "but it is not as easily used to visualize information."
                                                  "\n\n"
                                                  "Physical units: \\si{\\meter}.")
    }
  }
}
