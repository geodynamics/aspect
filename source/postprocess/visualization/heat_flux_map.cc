/*
  Copyright (C) 2011 - 2018 by the authors of the ASPECT code.

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


#include <aspect/postprocess/visualization/heat_flux_map.h>
#include <aspect/postprocess/heat_flux_map.h>
#include <aspect/geometry_model/interface.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {

      template <int dim>
      std::pair<std::string, Vector<float> *>
      HeatFluxMap<dim>::execute () const
      {
        std::pair<std::string, Vector<float> *>
        return_value ("heat_flux_map",
                      new Vector<float>(this->get_triangulation().n_active_cells()));

        std::vector<std::vector<std::pair<double, double> > > heat_flux_and_area =
          internal::compute_heat_flux_through_boundary_faces (*this);

        // loop over all of the surface cells and evaluate the heatflux
        typename DoFHandler<dim>::active_cell_iterator
        cell = this->get_dof_handler().begin_active(),
        endc = this->get_dof_handler().end();

        for (; cell!=endc; ++cell)
          if (cell->is_locally_owned())
            {
              (*return_value.second)(cell->active_cell_index()) = 0;

              for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
                if (cell->at_boundary(f) &&
                    (this->get_geometry_model().translate_id_to_symbol_name (cell->face(f)->boundary_id()) == "top" ||
                     this->get_geometry_model().translate_id_to_symbol_name (cell->face(f)->boundary_id()) == "bottom"))
                  {
                    // add heatflow for this face
                    (*return_value.second)(cell->active_cell_index()) += heat_flux_and_area[cell->active_cell_index()][f].first /
                        heat_flux_and_area[cell->active_cell_index()][f].second;
                  }
            }

        return return_value;
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

      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(HeatFluxMap,
                                                  "heat flux map",
                                                  "A visualization output object that generates output for "
                                                  "the heat flux density across each boundary. The heat flux density "
                                                  "is computed in outward direction, i.e., from the domain to the "
                                                  "outside, using the formula $-k \\nabla T \\cdot \\mathbf n$, where "
                                                  "$k$ is the thermal conductivity as reported by the material "
                                                  "model, $T$ is the temperature, and $\\mathbf n$ is the outward "
                                                  "normal. Note that the quantity so computed does not include "
                                                  "any energy transported across the boundary by material "
                                                  "transport in cases where $\\mathbf u \\cdot \\mathbf n \\neq 0$."
                                                  "At the edge of the domain (e.g. top left corner) the heat flux "
                                                  "density is calculated as the sum of the heat flux across both "
                                                  "boundaries of the cell (e.g. top and left boundary) divided "
                                                  "by the sum of both face areas. The integrated heatflux for each "
                                                  "boundary can be obtained from the heat flux statistics postprocessor.")
    }
  }
}
