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



#include <aspect/postprocess/velocity_boundary_statistics.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace
    {
      /**
       * Given a string #s, return it in the form ' ("s")' if nonempty.
       * Otherwise just return the empty string itself.
       */
      std::string parenthesize_if_nonempty (const std::string &s)
      {
        if (s.size() > 0)
          return " (\"" + s + "\")";
        else
          return "";
      }
    }


    template <int dim>
    std::pair<std::string,std::string>
    VelocityBoundaryStatistics<dim>::execute (TableHandler &statistics)
    {
      // create a quadrature formula for the velocity.
      const QGauss<dim-1> quadrature_formula (this->get_fe().base_element(this->introspection().base_elements.velocities).degree+1);

      FEFaceValues<dim> fe_face_values (this->get_mapping(),
                                        this->get_fe(),
                                        quadrature_formula,
                                        update_values |
                                        update_q_points);

      std::vector<Tensor<1,dim> > velocities (fe_face_values.n_quadrature_points);

      std::map<types::boundary_id, double> local_max_vel;
      std::map<types::boundary_id, double> local_min_vel;

      const std::set<types::boundary_id>
      boundary_indicators
        = this->get_geometry_model().get_used_boundary_indicators ();
      for (std::set<types::boundary_id>::const_iterator
           p = boundary_indicators.begin();
           p != boundary_indicators.end(); ++p)
        {
          local_max_vel[*p] = -std::numeric_limits<double>::max();
          local_min_vel[*p] = std::numeric_limits<double>::max();
        }

      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();

      // for every surface face that is part of a geometry boundary
      // and that is owned by this processor,
      // compute the maximum and minimum velocity magnitude.
      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
            if (cell->face(f)->at_boundary())
              {
                fe_face_values.reinit (cell, f);
                // extract velocities
                fe_face_values[this->introspection().extractors.velocities].get_function_values (this->get_solution(),
                    velocities);
                // determine the max and min velocity on the face
                double local_max = -std::numeric_limits<double>::max();
                double local_min = std::numeric_limits<double>::max();
                for (unsigned int q=0; q<fe_face_values.n_quadrature_points; ++q)
                  {
                    const double vel_mag = velocities[q].norm();
                    local_max = std::max(vel_mag,
                                         local_max);
                    local_min = std::min(vel_mag,
                                         local_min);
                  }
                // then merge them with the min/max velocities we found for other faces with the same boundary indicator
                local_max_vel[cell->face(f)->boundary_indicator()] = std::max(local_max,
                                                                              local_max_vel[cell->face(f)->boundary_indicator()]);
                local_min_vel[cell->face(f)->boundary_indicator()] = std::min(local_min,
                                                                              local_min_vel[cell->face(f)->boundary_indicator()]);
              }

      // now communicate to get the global values
      std::map<types::boundary_id, double> global_max_vel;
      std::map<types::boundary_id, double> global_min_vel;
      {
        // first collect local values in the same order in which they are listed
        // in the set of boundary indicators
        std::vector<double> local_max_values;
        std::vector<double> local_min_values;
        for (std::set<types::boundary_id>::const_iterator
             p = boundary_indicators.begin();
             p != boundary_indicators.end(); ++p)
          {
            local_max_values.push_back (local_max_vel[*p]);
            local_min_values.push_back (-local_min_vel[*p]);
          }
        // then collect contributions from all processors
        // now do the reductions over all processors. we can use Utilities::MPI::max
        // for the maximal values. unfortunately, there is currently no matching
        // Utilities::MPI::min function, so we already negated the argument, now take the maximum
        // as well, then negate it all again
        std::vector<double> global_max_values;
        Utilities::MPI::max (local_max_values, this->get_mpi_communicator(), global_max_values);
        std::vector<double> global_min_values;
        Utilities::MPI::max (local_min_values, this->get_mpi_communicator(), global_min_values);

        // and now take them apart into the global map again
        unsigned int index = 0;
        for (std::set<types::boundary_id>::const_iterator
             p = boundary_indicators.begin();
             p != boundary_indicators.end(); ++p, ++index)
          {
            global_max_vel[*p] = global_max_values[index];
            global_min_vel[*p] = -global_min_values[index];
          }
      }

      // now add the computed max and min velocities to the statistics object
      // and create a single string that can be output to the screen
      std::ostringstream screen_text;
      unsigned int index = 0;
      for (std::map<types::boundary_id, double>::const_iterator
           p = global_max_vel.begin(), a = global_min_vel.begin();
           p != global_max_vel.end() && a != global_min_vel.end(); ++p, ++a, ++index)
        {
          if (this->convert_output_to_years() == true)
            {
              const std::string name_max = "Maximum velocity magnitude on boundary with indicator "
                                           + Utilities::int_to_string(p->first)
                                           + parenthesize_if_nonempty(this->get_geometry_model()
                                                                      .translate_id_to_symbol_name (p->first))
                                           + " (m/yr)";
              statistics.add_value (name_max, p->second*year_in_seconds);
              const std::string name_min = "Minimum velocity magnitude on boundary with indicator "
                                           + Utilities::int_to_string(a->first)
                                           + parenthesize_if_nonempty(this->get_geometry_model()
                                                                      .translate_id_to_symbol_name (p->first))
                                           + " (m/yr)";
              statistics.add_value (name_min, a->second*year_in_seconds);
              // also make sure that the other columns filled by the this object
              // all show up with sufficient accuracy and in scientific notation
              statistics.set_precision (name_max, 8);
              statistics.set_scientific (name_max, true);
              statistics.set_precision (name_min, 8);
              statistics.set_scientific (name_min, true);
            }
          else
            {
              const std::string name_max = "Maximum velocity magnitude on boundary with indicator "
                                           + Utilities::int_to_string(p->first)
                                           + parenthesize_if_nonempty(this->get_geometry_model()
                                                                      .translate_id_to_symbol_name (p->first))
                                           + " (m/s)";
              statistics.add_value (name_max, p->second);
              const std::string name_min = "Minimum velocity magnitude on boundary with indicator "
                                           + Utilities::int_to_string(a->first)
                                           + parenthesize_if_nonempty(this->get_geometry_model()
                                                                      .translate_id_to_symbol_name (p->first))
                                           + " (m/s)";
              statistics.add_value (name_min, a->second);
              // also make sure that the other columns filled by the this object
              // all show up with sufficient accuracy and in scientific notation
              statistics.set_precision (name_max, 8);
              statistics.set_scientific (name_max, true);
              statistics.set_precision (name_min, 8);
              statistics.set_scientific (name_min, true);
            }

          // finally have something for the screen
          screen_text.precision(4);
          if (this->convert_output_to_years() == true)
            {
              screen_text << p->second *year_in_seconds << " m/yr, "
                          << a->second *year_in_seconds << " m/yr"
                          << (index == global_max_vel.size()-1 ? "" : ", ");
            }
          else
            {
              screen_text << p->second << " m/s, "
                          << a->second << " m/s"
                          << (index == global_max_vel.size()-1 ? "" : ", ");
            }

        }

      return std::pair<std::string, std::string> ("Max and min velocity along boundary parts:",
                                                  screen_text.str());
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(VelocityBoundaryStatistics,
                                  "velocity boundary statistics",
                                  "A postprocessor that computes some statistics about "
                                  "the velocity along the boundaries. For each boundary "
                                  "indicator (see your geometry description for which boundary "
                                  "indicators are used), the min and max velocity magnitude is computed.")
  }
}
