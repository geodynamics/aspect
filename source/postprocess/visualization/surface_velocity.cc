/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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
#include <aspect/utilities.h>
#include <aspect/postprocess/visualization/surface_velocity.h>

namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      SurfaceVelocity<dim>::
      SurfaceVelocity ()
        :
        DataPostprocessorVector<dim> ("surface_velocity",
        			update_values | update_quadrature_points | update_gradients)
      {}

      template <int dim>
      void
      SurfaceVelocity<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double> > &computed_quantities) const
      {

    	Assert ((computed_quantities[0].size() == dim), ExcInternalError());
        auto cell = input_data.template get_cell<DoFHandler<dim> >();

        const double velocity_scaling_factor =
                      this->convert_output_to_years() ? year_in_seconds : 1.0;

        // We only want to output dynamic topography at the top and bottom
        // boundary, so only compute it if the current cell has
        // a face at the top or bottom boundary.
        bool cell_at_top_boundary = false;
        for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
          if (cell->at_boundary(f) &&
              (this->get_geometry_model().translate_id_to_symbol_name (cell->face(f)->boundary_id()) == "top"))
            cell_at_top_boundary = true;

        if (cell_at_top_boundary)
          {
            // add strain rate in future
            for (unsigned int q=0; q<computed_quantities.size(); ++q)
            	for (unsigned int d = 0; d < dim; ++d)
            	{
            		computed_quantities[q](d) = input_data.solution_values[q][d] * velocity_scaling_factor;
            	}
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(SurfaceVelocity,
                                                  "surface velocity",
                                                  "A visualization for caculating the surface velocity to compare"
												  "with the GPS velocities.")
    }
  }
}
