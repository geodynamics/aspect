/*
  Copyright (C) 2011 - 2020 by the authors of the ASPECT code.

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


#include <aspect/postprocess/visualization/temperature_anomaly.h>
#include <aspect/boundary_temperature/constant.h>
#include <aspect/lateral_averaging.h>
#include <aspect/geometry_model/interface.h>



namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      TemperatureAnomaly<dim>::
      TemperatureAnomaly ()
        :
        DataPostprocessorScalar<dim> ("temperature_anomaly",
                                      update_values | update_quadrature_points )
      {
      }

      template <int dim>
      void
      TemperatureAnomaly<dim>::
      update ()
      {
        std::vector<double> temperature_depth_average(n_slices);
        this->get_lateral_averaging().get_temperature_averages(temperature_depth_average);

        // Estimates of the lateral temperature average at each depth are required
        // for all cell depths, including those
        // shallower than the midpoint of the first slice or
        // deeper than the midpoint of the last slice.

        padded_temperature_depth_average.resize(n_slices+2);

        if ( extrapolate_surface )
          {
            padded_temperature_depth_average[0] = 2.*temperature_depth_average[0] - temperature_depth_average[1];
          }
        else
          {
            Assert( this->has_boundary_temperature(),ExcInternalError());
            // retrieve the minimum temperature from the fixed temperature boundary indicators
            const double surface_temperature = this->get_boundary_temperature_manager().minimal_temperature(this->get_fixed_temperature_boundary_indicators());
            padded_temperature_depth_average[0] = 2.*surface_temperature - temperature_depth_average[0];
          }
        if ( extrapolate_bottom )
          {
            padded_temperature_depth_average[n_slices+1] = 2.*temperature_depth_average[n_slices-1] - temperature_depth_average[n_slices-2];
          }
        else
          {
            const double bottom_temperature = this->get_boundary_temperature_manager().maximal_temperature(this->get_fixed_temperature_boundary_indicators());
            padded_temperature_depth_average[n_slices+1] = 2.*bottom_temperature - temperature_depth_average[n_slices-1];
          }
        std::copy ( temperature_depth_average.begin(), temperature_depth_average.end(), padded_temperature_depth_average.begin() + 1 );

      }


      template <int dim>
      void
      TemperatureAnomaly<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double> > &computed_quantities) const
      {
        const double max_depth = this->get_geometry_model().maximal_depth();
        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert (computed_quantities[0].size() == 1,                   ExcInternalError());
        Assert (input_data.solution_values[0].size() == this->introspection().n_components,           ExcInternalError());

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            const double temperature = input_data.solution_values[q][this->introspection().component_indices.temperature];
            const double depth = this->get_geometry_model().depth (input_data.evaluation_points[q]);
            // calculate the depth-average cell containing this point. Note that cell centers are offset +0.5 cells in depth.
            const double slice_depth = (depth*n_slices)/max_depth + 0.5;
            const unsigned int idx = static_cast<unsigned int>(slice_depth);
            const double fractional_slice = slice_depth - static_cast<double>(idx);
            Assert(idx<n_slices+1, ExcInternalError());
            const double depth_average_temperature= (1. - fractional_slice)*padded_temperature_depth_average[idx] + fractional_slice*padded_temperature_depth_average[idx+1];
            computed_quantities[q](0) = temperature - depth_average_temperature;
          }
      }

      template <int dim>
      void
      TemperatureAnomaly<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Visualization");
          {
            prm.enter_subsection("Temperature anomaly");
            {
              prm.declare_entry ("Number of depth slices","20",
                                 Patterns::Integer (1),
                                 "Number of depth slices used to define "
                                 "average temperature.");
              prm.declare_entry ("Use maximal temperature for bottom","true",
                                 Patterns::Bool(),
                                 "If true, use the specified boundary temperatures as average temperatures at the surface. "
                                 "If false, extrapolate the temperature gradient between the first and second cells to the surface. "
                                 "This option will only work for models with a fixed surface temperature. ");
              prm.declare_entry ("Use minimal temperature for surface","true",
                                 Patterns::Bool(),
                                 "Whether to use the minimal specified boundary temperature as the bottom boundary temperature. "
                                 "This option will only work for models with a fixed bottom boundary temperature. ");

            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
      template <int dim>
      void
      TemperatureAnomaly<dim>::parse_parameters(ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Visualization");
          {
            prm.enter_subsection("Temperature anomaly");
            {
              n_slices = prm.get_integer("Number of depth slices");
              extrapolate_surface = !prm.get_bool("Use minimal temperature for surface");
              extrapolate_bottom = !prm.get_bool("Use maximal temperature for bottom");
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(TemperatureAnomaly,
                                                  "temperature anomaly",
                                                  "A visualization output postprocessor that outputs the temperature minus the depth-average of the temperature."
                                                  "The average temperature is calculated using the lateral averaging function from the ``depth average'' "
                                                  "postprocessor and interpolated linearly between the layers specified through ``Number of depth slices''")
    }
  }
}
