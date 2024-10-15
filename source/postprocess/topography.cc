/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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


#include <aspect/postprocess/topography.h>
#include <aspect/geometry_model/box.h>
#include <aspect/geometry_model/sphere.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/geometry_model/chunk.h>
#include <aspect/simulator.h>
#include <aspect/global.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

#include <cmath>
#include <limits>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    Topography<dim>::execute (TableHandler &statistics)
    {
      //Disallow use of the plugin with sphere geometry model
      AssertThrow(!Plugins::plugin_type_matches<const GeometryModel::Sphere<dim>>(this->get_geometry_model()),
                  ExcMessage("Topography postprocessor is not yet implemented "
                             "for the sphere geometry model. "
                             "Consider using a box, spherical shell, or chunk.") );

      const types::boundary_id relevant_boundary = this->get_geometry_model().translate_symbolic_boundary_name_to_id ("top");

      // Define a quadrature rule that has points only on the vertices of
      // a face, so that we can evaluate the solution at the surface
      // nodes.
      const QTrapezoid<dim-1> face_corners;
      FEFaceValues<dim> face_vals (this->get_mapping(), this->get_fe(), face_corners, update_quadrature_points);

      // have a stream into which we write the data. the text stream is then
      // later sent to processor 0
      std::ostringstream output_stats;
      std::ostringstream output_file;

      // On processor 0, write the file header
      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {
          output_file << "# "
                      << ((dim==2)? "x y" : "x y z")
                      << " topography" << std::endl;
        }

      // Choose stupidly large values for initialization
      double local_max_height = std::numeric_limits<double>::lowest();
      double local_min_height = std::numeric_limits<double>::max();

      // loop over all of the surface cells and save the elevation to stored_value
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned() && cell->at_boundary())
          for (const unsigned int face_no : cell->face_indices())
            if (cell->face(face_no)->at_boundary())
              {
                if ( cell->face(face_no)->boundary_id() != relevant_boundary)
                  continue;

                face_vals.reinit( cell, face_no);

                for (unsigned int corner = 0; corner < face_corners.size(); ++corner)
                  {
                    const Point<dim> vertex = face_vals.quadrature_point(corner);
                    const double elevation = this->get_geometry_model().height_above_reference_surface(vertex);
                    if (write_to_file)
                      output_file << vertex << ' '<< elevation << std::endl;
                    if ( elevation > local_max_height)
                      local_max_height = elevation;
                    if ( elevation < local_min_height)
                      local_min_height = elevation;
                  }
              }

      //Calculate min/max topography across all processes
      const double max_topography = Utilities::MPI::max(local_max_height, this->get_mpi_communicator());
      const double min_topography = Utilities::MPI::min(local_min_height, this->get_mpi_communicator());

      //Write results to statistics file
      statistics.add_value ("Minimum topography (m)",
                            min_topography);
      statistics.add_value ("Maximum topography (m)",
                            max_topography);
      const char *columns[] = { "Minimum topography (m)",
                                "Maximum topography (m)"
                              };
      for (auto &column : columns)
        {
          statistics.set_precision (column, 8);
          statistics.set_scientific (column, true);
        }

      output_stats.precision(4);
      output_stats << min_topography << " m, "
                   << max_topography << " m";

      // Write the solution to file

      // if this is the first time we get here, set the last output time
      // to the current time - output_interval. this makes sure we
      // always produce data during the first time step
      if (std::isnan(last_output_time))
        {
          last_output_time = this->get_time() - output_interval;
        }

      // Just return stats if text output is not required at all or not needed at this time
      if (!write_to_file || ((this->get_time() < last_output_time + output_interval)
                             && (this->get_timestep_number() != 0)))
        return std::pair<std::string, std::string> ("Topography min/max:",
                                                    output_stats.str());

      std::string filename = this->get_output_directory() +
                             "topography." +
                             Utilities::int_to_string(this->get_timestep_number(), 5);
      if (this->get_parameters().run_postprocessors_on_nonlinear_iterations)
        filename.append("." + Utilities::int_to_string (this->get_nonlinear_iteration(), 4));

      Utilities::collect_and_write_file_content(filename,
                                                output_file.str(),
                                                this->get_mpi_communicator());

      // if output_interval is positive, then update the last supposed output
      // time
      if (output_interval > 0)
        {
          // We need to find the last time output was supposed to be written.
          // this is the last_output_time plus the largest positive multiple
          // of output_intervals that passed since then. We need to handle the
          // edge case where last_output_time+output_interval==current_time,
          // we did an output and std::floor sadly rounds to zero. This is done
          // by forcing std::floor to round 1.0-eps to 1.0.
          const double magic = 1.0+2.0*std::numeric_limits<double>::epsilon();
          last_output_time = last_output_time + std::floor((this->get_time()-last_output_time)/output_interval*magic) * output_interval/magic;
        }

      return std::pair<std::string, std::string> ("Topography min/max:",
                                                  output_stats.str());
    }

    template <int dim>
    void
    Topography<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Topography");
        {
          prm.declare_entry ("Output to file", "false",
                             Patterns::Bool(),
                             "Whether or not to write topography to a text file named named "
                             "'topography.NNNNN' in the output directory");

          prm.declare_entry ("Time between text output", "0.",
                             Patterns::Double (0.),
                             "The time interval between each generation of "
                             "text output files. A value of zero indicates "
                             "that output should be generated in each time step. "
                             "Units: years if the "
                             "'Use years in output instead of seconds' parameter is set; "
                             "seconds otherwise.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    Topography<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Topography");
        {
          write_to_file        = prm.get_bool ("Output to file");
          output_interval = prm.get_double ("Time between text output");
          if (this->convert_output_to_years())
            output_interval *= year_in_seconds;
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }




  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(Topography,
                                  "topography",
                                  "A postprocessor intended for use with a deforming top surface.  After every step "
                                  "it loops over all the vertices on the top surface and determines the "
                                  "maximum and minimum topography relative to a reference datum (initial "
                                  "box height for a box geometry model or initial radius for a "
                                  "sphere/spherical shell geometry model). "
                                  "If 'Topography.Output to file' is set to true, also outputs topography "
                                  "into text files named `topography.NNNNN' in the output directory, "
                                  "where NNNNN is the number of the time step.\n"
                                  "The file format then consists of lines with Euclidean coordinates "
                                  "followed by the corresponding topography value."
                                  "Topography is printed/written in meters."
                                  "\n\n"
                                  "It is worth comparing this postprocessor with the "
                                  "visualization postprocessor called ``surface elevation''. The latter is "
                                  "used to *visualize* the surface elevation in graphical form, by "
                                  "outputting it into the same files that the solution and other "
                                  "postprocessed variables are written, to then be used for "
                                  "visualization using programs such as VisIt or Paraview. In "
                                  "contrast, the current postprocessor generates the same *kind* "
                                  "of information, but instead writes it as a point cloud that "
                                  "can then more easily be processed using tools other "
                                  "than visualization programs.")
  }
}
