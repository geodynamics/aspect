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
                                                  "Physical units: \\si{\\meter}.")
    }
  }
}
