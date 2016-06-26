/*
  Copyright (C) 2011 - 2014 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
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

        // loop over all of the surface cells and if one less than h/3 away from
        // the top surface, evaluate the stress at its center
        typename DoFHandler<dim>::active_cell_iterator
        cell = this->get_dof_handler().begin_active(),
        endc = this->get_dof_handler().end();

        unsigned int cell_index = 0;
        for (; cell!=endc; ++cell)
          if (cell->is_locally_owned())
            if (cell->at_boundary())
              {
                // see if the cell is at the *top* boundary, not just any boundary
                unsigned int top_face_idx = numbers::invalid_unsigned_int;
                {
                  for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
                    if (cell->at_boundary(f) && this->get_geometry_model().depth (cell->face(f)->center()) < cell->face(f)->minimum_vertex_distance()/3)
                      {
                        top_face_idx = f;
                        break;
                      }
                }
                if (top_face_idx == numbers::invalid_unsigned_int)
                  {
                    (*return_value.second)(cell_index) = 0;
                    continue;
                  }
                fe_face_values.reinit (cell, top_face_idx);


                // get the various components of the solution, then
                // evaluate the material properties there
                fe_face_values[this->introspection().extractors.temperature].get_function_gradients (this->get_solution(),
                    temperature_gradients);

                for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                  fe_face_values[this->introspection().extractors.compositional_fields[c]].get_function_values(this->get_solution(),
                      composition_values[c]);


                for (unsigned int i=0; i<fe_face_values.n_quadrature_points; ++i)
                  {
                    for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                      in.composition[i][c] = composition_values[c][i];
                  }

                this->get_material_model().evaluate(in, out);

                // Calculate the normal conductive heat flux given by the formula
                //   j = - k * n . grad T

                double normal_flux = 0;
                for (unsigned int q=0; q<fe_face_values.n_quadrature_points; ++q)
                  {
                    const double thermal_conductivity
                      = out.thermal_conductivities[q];

                    normal_flux += -thermal_conductivity *
                                   (temperature_gradients[q] * fe_face_values.normal_vector(q)) *
                                   fe_face_values.JxW(q);
                  }

                // get the location for each gridpoint at the top boundary
                fe_face_values.reinit(cell, top_face_idx);
                const Point<dim> midpoint_at_surface = cell->face(top_face_idx)->center();

                (*return_value.second)(cell_index) = normal_flux;

              }

        return return_value;
      }

      template <int dim>
      void
      HeatFluxMap<dim>::
      declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
        }
        prm.leave_subsection();
      }

      template <int dim>
      void
      HeatFluxMap<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
        }
        prm.leave_subsection();
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
                                                  "the heatflux across the top boundary. The heat flux is computed in "
                                                  "outward direction, i.e., from the domain to the outside, "
                                                  "using the formula $k \\nabla T \\cdot \\mathbf n$, where "
                                                  "$k$ is the thermal conductivity as reported by the material "
                                                  "model, $T$ is the temperature, and $\\mathbf n$ is the outward "
                                                  "normal. Note that the quantity so computed does not include "
                                                  "any energy transported across the boundary by material "
                                                  "transport in cases where $\\mathbf u \\cdot \\mathbf n \\neq 0$.")
    }
  }
}
