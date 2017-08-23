/*
  Copyright (C) 2011 - 2016 by the authors of the ASPECT code.

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
#include <aspect/simulator_access.h>

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

        // create a quadrature formula based on the temperature element alone.
        const QGauss<dim-1> quadrature_formula (this->get_fe().base_element(this->introspection().base_elements.temperature).degree+1);

        FEFaceValues<dim> fe_face_values (this->get_mapping(),
                                          this->get_fe(),
                                          quadrature_formula,
                                          update_gradients      | update_values |
                                          update_normal_vectors |
                                          update_q_points       | update_JxW_values);

        typename MaterialModel::Interface<dim>::MaterialModelInputs in(fe_face_values.n_quadrature_points, this->n_compositional_fields());
        typename MaterialModel::Interface<dim>::MaterialModelOutputs out(fe_face_values.n_quadrature_points, this->n_compositional_fields());

        std::vector<Tensor<1,dim> > temperature_gradients (quadrature_formula.size());
        std::vector<std::vector<double> > composition_values (this->n_compositional_fields(),std::vector<double> (quadrature_formula.size()));

        // loop over all of the surface cells and evaluate the heatflux
        typename DoFHandler<dim>::active_cell_iterator
        cell = this->get_dof_handler().begin_active(),
        endc = this->get_dof_handler().end();

        unsigned int cell_index = 0;
        for (; cell!=endc; ++cell,++cell_index)
          if (cell->is_locally_owned())
            if (cell->at_boundary())
              {
                double normal_flux = 0;
                double face_area = 0;

                for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
                  if (cell->at_boundary(f))
                    {
                      fe_face_values.reinit (cell, f);

                      // Temperature gradients needed for heat flux.
                      fe_face_values[this->introspection().extractors.temperature].get_function_gradients (this->get_solution(),
                          temperature_gradients);

                      // Set use_strain_rate to false since we don't need viscosity.
                      in.reinit(fe_face_values, cell, this->introspection(), this->get_solution(), false);
                      this->get_material_model().evaluate(in, out);


                      // Calculate the normal conductive heat flux given by the formula
                      //   j = - k * n . grad T

                      for (unsigned int q=0; q<fe_face_values.n_quadrature_points; ++q)
                        {
                          const double thermal_conductivity
                            = out.thermal_conductivities[q];

                          normal_flux += -thermal_conductivity *
                                         (temperature_gradients[q] * fe_face_values.normal_vector(q)) *
                                         fe_face_values.JxW(q);
                          face_area += fe_face_values.JxW(q);
                        }
                    }

                // store final position and heatflow
                (*return_value.second)(cell_index) = normal_flux / face_area;

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
