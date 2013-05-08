/*
  Copyright (C) 2011, 2012 by the authors of the ASPECT code.

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
/*  $Id: temperature_statistics.cc 1252 2012-10-05 19:52:00Z dannberg $  */


#include <aspect/postprocess/plume_statistics.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    PlumeStatistics<dim>::execute (TableHandler &statistics)
    {
      // create a quadrature formula based on the temperature element alone.
      // be defensive about determining that what we think is the temperature
      // element is it in fact
      Assert (this->get_dof_handler().get_fe().n_base_elements() == 3+(this->n_compositional_fields()>0 ? 1 : 0),
              ExcNotImplemented());


      // use a quadrature formula that has one point at
      // the location of each degree of freedom in the
      // temperature element
      const QIterated<dim> quadrature_formula (QTrapez<1>(),
    		                                   this->get_dof_handler().get_fe().base_element(2).degree);
      const unsigned int n_q_points = quadrature_formula.size();

      FEValues<dim> fe_values (this->get_mapping(),
    		                   this->get_fe(),
                               quadrature_formula,
                               update_values   |
                               update_quadrature_points |
                               update_JxW_values);

      const FEValuesExtractors::Scalar temperature (dim+1);
      std::vector<double> temperature_values(n_q_points);
      std::vector<double> compositional_values(n_q_points);

      // also need position of the quadrature points
      Point<dim,double> position;
      Point<dim,double> boundary_point (true);
      double depth;
      const double analysis_depth = 300000;
      const double analysis_temperature = 100;
      double maximal_depth = this->get_geometry_model().maximal_depth();

      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();

      double local_temperature_integral = 0;
      double local_compositional_integral = 0;
      double local_volume_integral = 0;
      double local_max_radius = 0;
      double local_max_temperature = 0;

      // compute the integral quantities by quadrature
      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);
            fe_values[this->introspection().extractors.temperature].get_function_values (this->get_solution(),
                                                                                         temperature_values);
            if (this->n_compositional_fields()>0)
              fe_values[this->introspection().extractors.compositional_fields[0]].get_function_values (this->get_solution(),
                  compositional_values);

            for (unsigned int q=0; q<n_q_points; ++q)
            {
              // check if depth is less than 300 km
              position = fe_values.quadrature_point(q);
              depth = this->get_geometry_model().depth(position);
              const double nonadiabatic_temperature = temperature_values[q] - this->get_adiabatic_conditions().temperature(position);

              if(depth < analysis_depth && nonadiabatic_temperature > analysis_temperature)
              {
                // integrate temperature and volume
                // record maximal temperature and radius (local maxima)
                local_temperature_integral += nonadiabatic_temperature*fe_values.JxW(q);
                if (this->n_compositional_fields()>0)
                  local_compositional_integral += compositional_values[q]*fe_values.JxW(q);
                local_volume_integral += fe_values.JxW(q);
                boundary_point[dim-1] = maximal_depth-depth;
                local_max_radius = std::max(local_max_radius,position.distance(boundary_point));
                local_max_temperature = std::max(local_max_temperature,nonadiabatic_temperature);
              }
            }
          }

      // compute the sum over all processors
      const double global_temperature_integral
        = Utilities::MPI::sum (local_temperature_integral, this->get_mpi_communicator());
      const double global_compositional_integral
        = Utilities::MPI::sum (local_compositional_integral, this->get_mpi_communicator());
      const double global_volume_integral
        = Utilities::MPI::sum (local_volume_integral, this->get_mpi_communicator());

      double global_max_temperature = 0;
      double global_max_radius = 0;

      // now do the reductions that are
      // min/max operations. do them in
      // one communication by multiplying
      // one value by -1
      {
        double local_values[2] = {local_max_temperature, local_max_radius};
        double global_values[2];

        Utilities::MPI::max (local_values, this->get_mpi_communicator(), global_values);

        global_max_temperature = global_values[0];
        global_max_radius = global_values[1];
      }

      // maximal nonadiabatic plume temperature
      statistics.add_value ("Maximal plume temperature (K)",
                            global_max_temperature);
      // average nonadiabatic plume temperature
      statistics.add_value ("Average plume temperature (K)",
                            global_temperature_integral / global_volume_integral);
      // average plume composition
      statistics.add_value ("Average plume composition",
        		            global_compositional_integral / global_volume_integral);
      statistics.add_value ("Plume volume (m³)",
    		                global_volume_integral);
      statistics.add_value ("Maximal radius",
    		                global_max_radius);

      // also make sure that the other columns filled by the this object
      // all show up with sufficient accuracy and in scientific notation
      {
        const char *columns[] = { "Maximal temperature (K)",
                                  "Average temperature (K)",
                                  "Average plume composition",
                                  "Plume volume (m³)",
                                  "Maximal radius"
                                };
        for (unsigned int i=0; i<sizeof(columns)/sizeof(columns[0]); ++i)
          {
            statistics.set_precision (columns[i], 8);
            statistics.set_scientific (columns[i], true);
          }
      }

      std::ostringstream output;
      output.precision(4);
      output << global_max_temperature << " K, "
             << global_temperature_integral / global_volume_integral << " K, "
		     << global_compositional_integral << ", "
		     << global_volume_integral << " m³";

      return std::pair<std::string, std::string> ("Plume temperature max/avg/Composition/Volume:",
                                                  output.str());
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(PlumeStatistics,
                                  "plume statistics",
                                  "A postprocessor that computes some statistics about "
                                  "a rising mantle plume, such as volume, average"
                                  "temperature etc.")
  }
}
