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


#include <aspect/postprocess/visualization/averaged_strain_rate.h>

namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {


      template <int dim>
      AveragedStrainRate<dim>::AveragedStrainRate()
      {}

      template <int dim>
      std::pair<std::string, Vector<float> *>
      AveragedStrainRate<dim>::execute() const
      {
        std::pair<std::string, Vector<float> *>
        return_value ("averaged_strain_rate",
                      new Vector<float>(this->get_triangulation().n_active_cells()));

        const unsigned int quadrature_degree = 5;
        const QMidpoint<1> quadrature_formula_mp;
        const QIterated<dim> quadrature_formula(quadrature_formula_mp, quadrature_degree);
        const unsigned int n_q_points = quadrature_formula.size();

        FEValues<dim> fe_values(this->get_mapping(), this->get_fe(),
                                quadrature_formula,
                                update_gradients | update_quadrature_points | update_JxW_values);

        // Set up cell iterator for looping
        typename DoFHandler<dim>::active_cell_iterator
        cell = this->get_dof_handler().begin_active(),
        endc = this->get_dof_handler().end();


        std::vector<SymmetricTensor<2,dim> > strain_rates(n_q_points);

        for (; cell != endc; ++cell)
          {
            if (cell->is_locally_owned())
              {
                fe_values.reinit(cell);
                fe_values[this->introspection().extractors.velocities].get_function_symmetric_gradients (this->get_solution(), strain_rates);

                double area = 0.0;
                double strain_rate_value = 0.0;

                for (unsigned int q=0; q<n_q_points; ++q)
                  {
                    const SymmetricTensor<2,dim> strain_rate = strain_rates[q];
                    const SymmetricTensor<2,dim> compressible_strain_rate
                      = (this->get_material_model().is_compressible()
                         ?
                         strain_rate - 1./3 * trace(strain_rate) * unit_symmetric_tensor<dim>()
                         :
                         strain_rate);

                    strain_rate_value += std::sqrt(compressible_strain_rate*compressible_strain_rate)
                                         * fe_values.JxW(q);

                    area += fe_values.JxW(q);

                  }

                (*return_value.second)(cell->active_cell_index()) = strain_rate_value/area;
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(AveragedStrainRate,
                                                  "averaged strain rate",
                                                  "A visualization output object that generates the "
                                                  "the average norm of the strain rate , i.e., for the quantity "
                                                  "$\\sqrt{\\varepsilon(\\mathbf u):\\varepsilon(\\mathbf u)}$ "
                                                  "in the incompressible case and "
                                                  "$\\sqrt{[\\varepsilon(\\mathbf u)-\\tfrac 13(\\textrm{tr}\\;\\varepsilon(\\mathbf u))\\mathbf I]:"
                                                  "[\\varepsilon(\\mathbf u)-\\tfrac 13(\\textrm{tr}\\;\\varepsilon(\\mathbf u))\\mathbf I]}$ "
                                                  "in the compressible case.")
    }
  }
}
