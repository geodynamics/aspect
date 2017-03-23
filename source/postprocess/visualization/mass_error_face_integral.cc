/*
  Copyright (C) 2017 by the authors of the ASPECT code.

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


#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <aspect/postprocess/visualization/mass_error_face_integral.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      std::pair<std::string, Vector<float> *>
      MassErrorFaceIntegral<dim>::execute() const
      {
        std::pair<std::string, Vector<float> *>
        return_value ("mass_error_face_integral",
                      new Vector<float>(this->get_triangulation().n_active_cells()));

        const unsigned int quadrature_degree = this->get_fe().base_element(this->introspection().base_elements.velocities).degree+1;

        const QGauss<dim-1> quadrature_formula_face (quadrature_degree);
        const QGauss<dim> quadrature_formula (quadrature_degree);

        const unsigned int n_q_points = quadrature_formula.size();
        const unsigned int n_q_points_face = quadrature_formula_face.size();

        /* We need three different FEValues objects for the operations below.
         * fe_values_support evaluates the material model at the support points
         * of the hypothetical density element declared above.
         * fe_values directly evaluates the velocity at the quadrature points
         * fe_values_density is used to interpolate the density from the support
         * points to the quadrature points.
         */

        FEFaceValues<dim> fe_face_values (this->get_mapping(),
                                          this->get_fe(),
                                          quadrature_formula_face,
                                          update_values |
                                          update_gradients |
                                          update_JxW_values |
                                          update_normal_vectors |
                                          update_q_points);

        FEValues<dim> fe_values (this->get_mapping(),
                                 this->get_fe(),
                                 quadrature_formula,
                                 update_JxW_values);

        std::vector<std::vector<double> > composition_values (this->n_compositional_fields(),std::vector<double> (quadrature_formula.size()));

        MaterialModel::MaterialModelInputs<dim> in(n_q_points_face, this->n_compositional_fields());
        MaterialModel::MaterialModelOutputs<dim> out(n_q_points_face, this->n_compositional_fields());

        typename DoFHandler<dim>::active_cell_iterator
        cell = this->get_dof_handler().begin_active(),
        endc = this->get_dof_handler().end();
        unsigned int cell_index = 0;
        for (; cell!=endc; ++cell,++cell_index)
          if (cell->is_locally_owned())
            {
              fe_values.reinit(cell);

              double cell_volume = 0.0;

              for (unsigned int q = 0; q < n_q_points; ++q)
                cell_volume += fe_values.JxW(q);

              double div_rho_u = 0.0;

              for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
                {
                  fe_face_values.reinit(cell,f);

                  fe_face_values[this->introspection().extractors.velocities].get_function_values (this->get_solution(),
                      in.velocity);
                  fe_face_values[this->introspection().extractors.temperature].get_function_values (this->get_solution(),
                      in.temperature);
                  fe_face_values[this->introspection().extractors.pressure].get_function_values (this->get_solution(),
                                                                                                 in.pressure);
                  fe_face_values[this->introspection().extractors.pressure].get_function_gradients (this->get_solution(),
                      in.pressure_gradient);
                  for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                    fe_face_values[this->introspection().extractors.compositional_fields[c]].get_function_values(this->get_solution(),
                        composition_values[c]);

                  in.position = fe_face_values.get_quadrature_points();

                  in.strain_rate.resize(0);

                  for (unsigned int i=0; i<n_q_points_face; ++i)
                    {
                      for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                        in.composition[i][c] = composition_values[c][i];
                    }
                  in.cell = &cell;

                  this->get_material_model().evaluate(in, out);

                  for (unsigned int q = 0; q < n_q_points_face; ++q)
                    div_rho_u  += out.densities[q] * (in.velocity[q] * fe_face_values.normal_vector(q))
                                  * fe_face_values.JxW(q);
                }

              (*return_value.second)(cell_index) = div_rho_u / cell_volume;
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(MassErrorFaceIntegral,
                                                  "mass error face integral",
                                                  "A postprocessor that computes some statistics about the "
                                                  "error of the mass conservation equation. More precisely, it "
                                                  "computes the norm and the integral of the equation "
                                                  "$\\rho \\Nabla \\mathbf u + \\Nabla \\rho \\cdot \\mathbf u$. "
                                                  "Therefore, its results indicate how accurate the chosen approximation "
                                                  "for compressibility is for the given model setup.")
    }
  }
}
