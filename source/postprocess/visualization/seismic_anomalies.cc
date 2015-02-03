/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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


#include <aspect/postprocess/visualization/seismic_anomalies.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      namespace
      {
        /**
         * Computes a running average of a vector.
         * http://en.wikipedia.org/wiki/Moving_average
         */
        void
        compute_running_average(std::vector<double> &values,
                                const int ncells)
        {
          std::vector<double> temp(values.size());

          for (unsigned int idx=0; idx<values.size(); idx++)
            {
              double sum = 0;
              for (int isum=-ncells; isum<=ncells; isum++)
                {
                  sum += values[std::max(0,std::min((int) values.size()-1,(int) idx+isum))];
                }
              temp[idx] = sum/((double) (ncells*2+1));
            }
          values = temp;
        }



        /**
         * Copy the values of the compositional fields at the quadrature point
         * q given as input parameter to the output vector
         * composition_values_at_q_point.
         */
        void
        extract_composition_values_at_q_point (const std::vector<std::vector<double> > &composition_values,
                                               const unsigned int                      q,
                                               std::vector<double>                    &composition_values_at_q_point)
        {
          for (unsigned int k=0; k < composition_values_at_q_point.size(); ++k)
            composition_values_at_q_point[k] = composition_values[k][q];
        }
      }



      template <int dim>
      std::pair<std::string, Vector<float> *>
      SeismicVsAnomaly<dim>::execute() const
      {
        std::pair<std::string, Vector<float> *>
        return_value ("Vs_anomaly",
                      new Vector<float>(this->get_triangulation().n_active_cells()));


        const unsigned int npoints = 2; // window in running average half-width of window
        std::vector<double> Vs_depth_average(50);

        this->get_depth_average_Vs(Vs_depth_average);
        compute_running_average(Vs_depth_average, npoints);

        const unsigned int num_slices = Vs_depth_average.size();
        const double max_depth = this->get_geometry_model().maximal_depth();

        // evaluate a single point per cell
        const QMidpoint<dim> quadrature_formula;
        const unsigned int n_q_points = quadrature_formula.size();

        FEValues<dim> fe_values (this->get_mapping(),
                                 this->get_fe(),
                                 quadrature_formula,
                                 update_values   |
                                 update_quadrature_points );

        std::vector<double> pressure_values(n_q_points);
        std::vector<double> temperature_values(n_q_points);
        std::vector<std::vector<double> > composition_values (this->n_compositional_fields(),
                                                              std::vector<double> (n_q_points));
        std::vector<double> composition_values_at_q_point (this->n_compositional_fields());

        typename DoFHandler<dim>::active_cell_iterator
        cell = this->get_dof_handler().begin_active(),
        endc = this->get_dof_handler().end();

        unsigned int cell_index = 0;
        for (; cell!=endc; ++cell,++cell_index)
          if (cell->is_locally_owned())
            {
              fe_values.reinit (cell);
              fe_values[this->introspection().extractors.pressure].get_function_values (this->get_solution(),
                                                                                        pressure_values);
              fe_values[this->introspection().extractors.temperature].get_function_values (this->get_solution(),
                                                                                           temperature_values);
              for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                fe_values[this->introspection().extractors.compositional_fields[c]].get_function_values(this->get_solution(),
                    composition_values[c]);


              extract_composition_values_at_q_point (composition_values,
                                                     0,
                                                     composition_values_at_q_point);

              const double Vs = this->get_material_model().seismic_Vs(temperature_values[0],
                                                                      pressure_values[0],
                                                                      composition_values_at_q_point,
                                                                      fe_values.quadrature_point(0));
              const double depth = this->get_geometry_model().depth(fe_values.quadrature_point(0));
              const unsigned int idx = static_cast<unsigned int>((depth*num_slices)/max_depth);
              Assert(idx<num_slices, ExcInternalError());

              // compute the deviation from the average in per cent
              (*return_value.second)(cell_index) = (Vs - Vs_depth_average[idx])/Vs_depth_average[idx]*1e2;
            }

        return return_value;
      }


      template <int dim>
      std::pair<std::string, Vector<float> *>
      SeismicVpAnomaly<dim>::execute() const
      {
        std::pair<std::string, Vector<float> *>
        return_value ("Vp_anomaly",
                      new Vector<float>(this->get_triangulation().n_active_cells()));


        const unsigned int npoints = 2; // window in running average half-width of window
        std::vector<double> Vp_depth_average(50);

        this->get_depth_average_Vp(Vp_depth_average);
        compute_running_average(Vp_depth_average, npoints);

        const unsigned int num_slices = Vp_depth_average.size();
        const double max_depth = this->get_geometry_model().maximal_depth();

        // evaluate a single point per cell
        const QMidpoint<dim> quadrature_formula;
        const unsigned int n_q_points = quadrature_formula.size();

        FEValues<dim> fe_values (this->get_mapping(),
                                 this->get_fe(),
                                 quadrature_formula,
                                 update_values   |
                                 update_quadrature_points );

        std::vector<double> pressure_values(n_q_points);
        std::vector<double> temperature_values(n_q_points);
        std::vector<std::vector<double> > composition_values (this->n_compositional_fields(),
                                                              std::vector<double> (n_q_points));
        std::vector<double> composition_values_at_q_point (this->n_compositional_fields());

        typename DoFHandler<dim>::active_cell_iterator
        cell = this->get_dof_handler().begin_active(),
        endc = this->get_dof_handler().end();

        unsigned int cell_index = 0;
        for (; cell!=endc; ++cell,++cell_index)
          if (cell->is_locally_owned())
            {
              fe_values.reinit (cell);
              fe_values[this->introspection().extractors.pressure].get_function_values (this->get_solution(),
                                                                                        pressure_values);
              fe_values[this->introspection().extractors.temperature].get_function_values (this->get_solution(),
                                                                                           temperature_values);
              for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                fe_values[this->introspection().extractors.compositional_fields[c]].get_function_values(this->get_solution(),
                    composition_values[c]);


              extract_composition_values_at_q_point (composition_values,
                                                     0,
                                                     composition_values_at_q_point);

              const double Vp = this->get_material_model().seismic_Vp(temperature_values[0],
                                                                      pressure_values[0],
                                                                      composition_values_at_q_point,
                                                                      fe_values.quadrature_point(0));
              const double depth = this->get_geometry_model().depth(fe_values.quadrature_point(0));
              const unsigned int idx = static_cast<unsigned int>((depth*num_slices)/max_depth);
              Assert(idx<num_slices, ExcInternalError());

              // compute the deviation from the average in per cent
              (*return_value.second)(cell_index) = (Vp - Vp_depth_average[idx])/Vp_depth_average[idx]*1e2;
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(SeismicVsAnomaly,
                                                  "Vs anomaly",
                                                  "A visualization output object that generates output "
                                                  "showing the anomaly in the seismic shear wave "
                                                  "speed $V_s$ as a spatially variable function with one "
                                                  "value per cell. This anomaly is shown as a percentage "
                                                  "change relative to the average value of $V_s$ at "
                                                  "the depth of this cell.")

      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(SeismicVpAnomaly,
                                                  "Vp anomaly",
                                                  "A visualization output object that generates output "
                                                  "showing the anomaly in the seismic compression wave "
                                                  "speed $V_p$ as a spatially variable function with one "
                                                  "value per cell. This anomaly is shown as a percentage "
                                                  "change relative to the average value of $V_p$ at "
                                                  "the depth of this cell.")
    }
  }
}
