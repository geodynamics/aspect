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


#include <aspect/simulator.h>
#include <aspect/global.h>
#include <aspect/mesh_deformation/free_surface.h>
#include <aspect/utilities.h>
#include <aspect/structured_data.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/postprocess/sea_level.h>
#include <aspect/postprocess/geoid.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

#include <cmath>
#include <limits>

// TODO:
// Include geoid anomaly of surface load.
// Compute new surface gravity per location, dependent on surface topography and Earth model instead of outer_radius of geometry model.
// Include option for time-changing ocean geometry, now it is a fixed ocean geometry.
// The ice height data is now only loaded once from the first time step. It is not constant over time. Implement an option to load data for current timestep (as in boundary traction).
// Implement iterative structure with solvers.

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    void
    SeaLevel<dim>::initialize()
    {
      Assert(false, ExcNotImplemented(""));
    }

    template <>
    void
    SeaLevel<3>::initialize()
    {
      const int dim = 3;

      topography_lookup = std::make_unique<Utilities::StructuredDataLookup<dim-1>>(/* n_components = */1,
                          /* scale = */1.0);
      topography_lookup->load_file(data_directory_topography+data_file_name_topography,this->get_mpi_communicator());

      ice_height_lookup = std::make_unique<Utilities::StructuredDataLookup<dim-1>>(/* n_components = */1,
                          /* scale = */1.0);
      ice_height_lookup->load_file(data_directory_ice_height+data_file_name_ice_height,this->get_mpi_communicator());
    }



    template <int dim>
    double
    SeaLevel<dim>::compute_nonuniform_sea_level_change(const Point<dim> &/*position*/) const
    {
      Assert(false, ExcNotImplemented());
      return 0;
    }

    template <>
    double
    SeaLevel<3>::compute_nonuniform_sea_level_change(const Point<3> &position) const
    {
      const int dim = 3;

      const Postprocess::Geoid<dim> &geoid =
        this->get_postprocess_manager().template get_matching_active_plugin<Postprocess::Geoid<dim>>();

      const double geoid_displacement = geoid.evaluate(position); // TODO: check sign of geoid_displacement
      const double topography = this->get_geometry_model().height_above_reference_surface(position);

      Point<dim-1> long_colat;
      const std::array<double,dim> spherical_position = aspect::Utilities::Coordinates::cartesian_to_spherical_coordinates(position);
      for (unsigned int d = 1; d < dim; ++d)
        long_colat[d-1] = spherical_position[d];

      const double topography_init = topography_lookup->get_data(long_colat,0);
      int ocean_mask = 0;
      if (topography_init < 0.0)
        ocean_mask = 1;
      else
        ocean_mask = 0;

      const double nonuniform_sea_level_change = (geoid_displacement-topography)*ocean_mask;

      return nonuniform_sea_level_change;
    }



    template <int dim>
    double
    SeaLevel<dim>::compute_sea_level_offset()
    {
      Assert(false, ExcNotImplemented());
      return 0;
    }

    template <>
    double
    SeaLevel<3>::compute_sea_level_offset()
    {
      const int dim = 3;

      const types::boundary_id top_boundary_id = this->get_geometry_model().translate_symbolic_boundary_name_to_id("top");

      const Postprocess::Geoid<dim> &geoid =
        this->get_postprocess_manager().template get_matching_active_plugin<Postprocess::Geoid<dim>>();

      const unsigned int quadrature_degree = this->introspection().polynomial_degree.temperature;
      const QGauss<dim-1> quadrature_formula_face(quadrature_degree);

      FEFaceValues<dim> fe_face_values (this->get_mapping(),
                                        this->get_fe(),
                                        quadrature_formula_face,
                                        update_values |
                                        update_quadrature_points |
                                        update_JxW_values);

      double integral_ocean_mask = 0;
      double integral_ice_height = 0;
      double integral_topo_geoid = 0;

      // Loop over all of the boundary cells and if one is at the
      // surface, evaluate the topography/geoid displacement/ocean mask/ice height there.
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned() && cell->at_boundary())
          {
            unsigned int face_idx = numbers::invalid_unsigned_int;
            bool at_upper_surface = false;
            {
              for (const unsigned int f : cell->face_indices())
                {
                  if (cell->at_boundary(f) && cell->face(f)->boundary_id() == top_boundary_id)
                    {
                      face_idx = f;
                      at_upper_surface = true;
                      break;
                    }
                }

              // If the cell is not at the top boundary, jump to the next cell.
              if (at_upper_surface == false)
                continue;
            }

            // Focus on the boundary cell's upper face if on the top boundary.
            fe_face_values.reinit(cell,face_idx);

            // If the cell is at the top boundary, add its contributions to the topography/geoid displacement/ocean mask/ice height storage vectors.
            for (unsigned int q=0; q<fe_face_values.n_quadrature_points; ++q)
              {
                const Point<dim> current_position = fe_face_values.quadrature_point(q);
                const double topography = this->get_geometry_model().height_above_reference_surface(current_position);

                const double geoid_displacement = geoid.evaluate(current_position);

                Point<dim-1> long_colat;
                const std::array<double,dim> spherical_position = aspect::Utilities::Coordinates::cartesian_to_spherical_coordinates(current_position);
                for (unsigned int d = 1; d < dim; ++d)
                  long_colat[d-1] = spherical_position[d];

                const double topography_init = topography_lookup->get_data(long_colat,0);
                int ocean_mask = 0;
                if (topography_init < 0.0)
                  ocean_mask = 1;
                else
                  ocean_mask = 0;

                const double ice_height = ice_height_lookup->get_data(long_colat,0);

                // Compute required integrals for the sea level offset.
                integral_ocean_mask += ocean_mask*fe_face_values.JxW(q);
                integral_ice_height += ice_height*density_ice*(1-ocean_mask)*fe_face_values.JxW(q);
                integral_topo_geoid += (geoid_displacement-topography)*ocean_mask*fe_face_values.JxW(q);
              }
          }

      integral_ocean_mask = Utilities::MPI::sum (integral_ocean_mask, this->get_mpi_communicator());
      integral_ice_height = Utilities::MPI::sum (integral_ice_height, this->get_mpi_communicator());
      integral_topo_geoid = Utilities::MPI::sum (integral_topo_geoid, this->get_mpi_communicator());

      const double sea_level_offset = -1./integral_ocean_mask*(1./density_water*integral_ice_height+integral_topo_geoid);

      return sea_level_offset;
    }



    template <int dim>
    double
    SeaLevel<dim>::compute_total_surface_pressure(const Point<dim> &/*position*/) const
    {
      Assert(false, ExcNotImplemented());
      return 0;
    }

    template <>
    double
    SeaLevel<3>::compute_total_surface_pressure(const Point<3> &position) const
    {
      // This function should be called from the boundary traction.
      const int dim = 3;

      const double nonuniform_sea_level_change = compute_nonuniform_sea_level_change(position);

      Point<dim-1> long_colat;
      const std::array<double,dim> spherical_position = aspect::Utilities::Coordinates::cartesian_to_spherical_coordinates(position);
      for (unsigned int d = 1; d < dim; ++d)
        long_colat[d-1] = spherical_position[d];
      const double ice_height = ice_height_lookup->get_data(long_colat,0);

      const GeometryModel::SphericalShell<dim> &geometry_model = Plugins::get_plugin_as_type<const GeometryModel::SphericalShell<dim>> (this->get_geometry_model());

      // TODO: instead of outer_radius use free surface topography to get different gravity per surface point
      const double outer_radius = geometry_model.outer_radius();
      Point<dim> surface_point;
      surface_point[0] = outer_radius;
      const double surface_gravity = this->get_gravity_model().gravity_vector(surface_point).norm();

      const double total_surface_pressure = surface_gravity*((nonuniform_sea_level_change+sea_level_offset)*density_water+ice_height*density_ice);

      return total_surface_pressure;
    }



    template <int dim>
    std::pair<std::string,std::string>
    SeaLevel<dim>::execute (TableHandler &/*statistics*/)
    {
      Assert(false, ExcNotImplemented());
      return {{""},{""}};
    }

    template <>
    std::pair<std::string,std::string>
    SeaLevel<3>::execute (TableHandler &statistics)
    {
      const int dim = 3;

      // Disallow use of the plugin for any other than a 3D spherical shell geometry.
      AssertThrow (Plugins::plugin_type_matches<const GeometryModel::SphericalShell<dim>>(this->get_geometry_model()) && dim == 3,
                   ExcMessage("The geoid postprocessor is currently only implemented for the 3D spherical shell geometry model."));

      const types::boundary_id top_boundary_id = this->get_geometry_model().translate_symbolic_boundary_name_to_id("top");

      const Postprocess::Geoid<dim> &geoid =
        this->get_postprocess_manager().template get_matching_active_plugin<Postprocess::Geoid<dim>>();

      // Get the sea level offset (constant for every location).
      sea_level_offset = compute_sea_level_offset();

      // Define a quadrature rule that has points only on the vertices of
      // a face, so that we can evaluate the solution at the surface
      // nodes.
      const QTrapezoid<dim-1> face_corners;

      FEFaceValues<dim> face_vals (this->get_mapping(), this->get_fe(), face_corners, update_quadrature_points);

      // Have a stream into which we write the data. The text stream is then later sent to processor 0.
      std::ostringstream output_stats;
      std::ostringstream output_file;

      // Choose stupidly large values for initialization.
      double local_max_height = -std::numeric_limits<double>::max();
      double local_min_height = std::numeric_limits<double>::max();

      // Loop over all of the surface cells and save the desired data to the output file.
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned() && cell->at_boundary())
          for (const unsigned int face_no : cell->face_indices())
            if (cell->face(face_no)->at_boundary())
              {
                if (cell->face(face_no)->boundary_id() != top_boundary_id)
                  continue;

                face_vals.reinit(cell, face_no);

                for (unsigned int corner = 0; corner < face_corners.size(); ++corner)
                  {
                    const Point<dim> vertex = face_vals.quadrature_point(corner);
                    const double nonuniform_sea_level_change = compute_nonuniform_sea_level_change(vertex);
                    const double geoid_displacement = geoid.evaluate(vertex); // TODO: check sign of geoid_displacement
                    const double topography = this->get_geometry_model().height_above_reference_surface(vertex);

                    Point<dim-1> long_colat;
                    const std::array<double,dim> spherical_position = aspect::Utilities::Coordinates::cartesian_to_spherical_coordinates(vertex);
                    for (unsigned int d = 1; d < dim; ++d)
                      long_colat[d-1] = spherical_position[d];

                    const double ice_height = ice_height_lookup->get_data(long_colat,0);
                    const double topography_init = topography_lookup->get_data(long_colat,0);
                    int ocean_mask = 0;
                    if (topography_init < 0.0)
                      ocean_mask = 1;
                    else
                      ocean_mask = 0;

                    if (write_to_file)
                      {
                        const std::array<double,dim> scoord = aspect::Utilities::Coordinates::cartesian_to_spherical_coordinates(vertex);
                        output_file << scoord[0] << ' ' << scoord[1] << ' ' << scoord[2] << ' ' << nonuniform_sea_level_change << ' ' << sea_level_offset << ' ' << geoid_displacement << ' ' << topography << ' ' << ice_height << ' ' << ocean_mask << std::endl;
                      }
                    if (nonuniform_sea_level_change > local_max_height)
                      local_max_height = nonuniform_sea_level_change;
                    if (nonuniform_sea_level_change < local_min_height)
                      local_min_height = nonuniform_sea_level_change;
                  }
              }

      // Calculate min/max non uniform sea_level change across all processes.
      const double max_nonuniform_sea_level_change = Utilities::MPI::max(local_max_height, this->get_mpi_communicator());
      const double min_nonuniform_sea_level_change = Utilities::MPI::min(local_min_height, this->get_mpi_communicator());

      // Write results to statistics file.
      statistics.add_value ("Minimum non-uniform sea_level change (m)",
                            min_nonuniform_sea_level_change);
      statistics.add_value ("Maximum non-uniform sea_level change (m)",
                            max_nonuniform_sea_level_change);
      const char *columns[] = { "Minimum non-uniform sea_level change (m)",
                                "Maximum non-uniform sea_level change (m)"
                              };
      for (const auto &column : columns)
        {
          statistics.set_precision (column, 8);
          statistics.set_scientific (column, true);
        }

      output_stats.precision(4);
      output_stats << min_nonuniform_sea_level_change << " m, "
                   << max_nonuniform_sea_level_change << " m";

      // Write the solution to file.

      // If this is the first time we get here, set the last output time
      // to the current time - output_interval. This makes sure we
      // always produce data during the first time step.
      if (std::isnan(last_output_time))
        {
          last_output_time = this->get_time() - output_interval;
        }

      // Just return stats if text output is not required at all or not needed at this time.
      if (!write_to_file || ((this->get_time() < last_output_time + output_interval)
                             && (this->get_timestep_number() != 0)))
        return {"Non-uniform sea level change min/max:", output_stats.str()};

      const unsigned int max_data_length = Utilities::MPI::max (output_file.str().size()+1,
                                                                this->get_mpi_communicator());

      const unsigned int mpi_tag = 777;

      // On processor 0, collect all of the data the individual processors sent
      // and concatenate them into one file.
      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {
          std::string filename = this->get_output_directory() +
                                 "nonuniform_sea_level_change." +
                                 Utilities::int_to_string(this->get_timestep_number(), 5);
          if (this->get_parameters().run_postprocessors_on_nonlinear_iterations)
            filename.append("." + Utilities::int_to_string (this->get_nonlinear_iteration(), 4));

          std::ofstream file (filename);

          file << "# "
               << "R phi theta"
               << " nonuniform_sea_level_change"
               << " sea_level_offset"
               << " geoid_displacement"
               << " topography"
               << " ice_height"
               << " ocean_mask" << std::endl;

          // First write out the data we have created locally.
          file << output_file.str();

          std::string tmp;
          tmp.resize (max_data_length, '\0');

          // Then loop through all of the other processors and collect
          // data, then write it to the file.
          for (unsigned int p=1; p<Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()); ++p)
            {
              MPI_Status status;
              // Get the data. Note that MPI says that an MPI_Recv may receive
              // less data than the length specified here. Since we have already
              // determined the maximal message length, we use this feature here
              // rather than trying to find out the exact message length with
              // a call to MPI_Probe.
              const int ierr = MPI_Recv (&tmp[0], max_data_length, MPI_CHAR, p, mpi_tag,
                                         this->get_mpi_communicator(), &status);
              AssertThrowMPI(ierr);

              // Output the string. Note that 'tmp' has length max_data_length,
              // but we only wrote a certain piece of it in the MPI_Recv, ended
              // by a \0 character. Write only this part by outputting it as a
              // C string object, rather than as a std::string.
              file << tmp.c_str();
            }
        }
      else
        // On other processors, send the data to processor zero. include the \0
        // character at the end of the string.
        {
          output_file << "\0";
          const int ierr = MPI_Send (&output_file.str()[0], output_file.str().size()+1, MPI_CHAR, 0, mpi_tag,
                                     this->get_mpi_communicator());
          AssertThrowMPI(ierr);
        }

      // if output_interval is positive, then update the last supposed output time
      if (output_interval > 0)
        {
          // We need to find the last time output was supposed to be written.
          // This is the last_output_time plus the largest positive multiple
          // of output_intervals that passed since then. We need to handle the
          // edge case where last_output_time+output_interval==current_time,
          // we did an output and std::floor sadly rounds to zero. This is done
          // by forcing std::floor to round 1.0-eps to 1.0.
          const double magic = 1.0+2.0*std::numeric_limits<double>::epsilon();
          last_output_time = last_output_time + std::floor((this->get_time()-last_output_time)/output_interval*magic) * output_interval/magic;
        }

      return {"Non-uniform sea level change min/max:", output_stats.str()};
    }



    template <int dim>
    std::list<std::string>
    SeaLevel<dim>::required_other_postprocessors() const
    {
      return {"geoid"};
    }



    template <int dim>
    void
    SeaLevel<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Sea level");
        {
          prm.declare_entry ("Ice density", "931",
                             Patterns::Double(0.),
                             "The density of ice [kg/m3]");
          prm.declare_entry ("Water density", "1000",
                             Patterns::Double(0.),
                             "The density of water [kg/m3]");

          prm.declare_entry ("Data directory topography",
                             "$ASPECT_SOURCE_DIR/data/postprocess/sea-level/",
                             Patterns::DirectoryName (),
                             "The name of a directory that contains the topography "
                             "ascii data. This path may either be absolute (if starting with a "
                             "`/') or relative to the current directory. The path may also "
                             "include the special text `$ASPECT_SOURCE_DIR' which will be "
                             "interpreted as the path in which the ASPECT source files were "
                             "located when ASPECT was compiled. This interpretation allows, "
                             "for example, to reference files located in the `data/' subdirectory "
                             "of ASPECT.");
          prm.declare_entry ("Data file name topography", "shell_3d_topo_top.0.txt",
                             Patterns::Anything (),
                             "The file name of the topography ascii data. For the ascii data, "
                             "provide file in the same format as described in "
                             "'ascii data' initial composition plugin." );

          prm.declare_entry ("Data directory ice height",
                             "$ASPECT_SOURCE_DIR/data/postprocess/sea-level/",
                             Patterns::DirectoryName (),
                             "The name of a directory that contains the ice height [m] "
                             "ascii data. This path may either be absolute (if starting with a "
                             "`/') or relative to the current directory. The path may also "
                             "include the special text `$ASPECT_SOURCE_DIR' which will be "
                             "interpreted as the path in which the ASPECT source files were "
                             "located when ASPECT was compiled. This interpretation allows, "
                             "for example, to reference files located in the `data/' subdirectory "
                             "of ASPECT.");
          prm.declare_entry ("Data file name ice height", "shell_3d_ice_top.0.txt",
                             Patterns::Anything (),
                             "The file name of the ice height ascii data. For the ascii data, "
                             "provide file in the same format as described in "
                             "'ascii data' initial composition plugin." );

          prm.declare_entry ("Output to file", "false",
                             Patterns::List(Patterns::Bool()),
                             "Whether or not to write sea level to a text file named named "
                             "'sea_level.NNNNN' in the output directory");
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
    SeaLevel<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Sea level");
        {
          density_ice = prm.get_double ("Ice density");
          density_water = prm.get_double ("Water density");
          data_directory_topography = Utilities::expand_ASPECT_SOURCE_DIR(prm.get ("Data directory topography"));
          data_file_name_topography = prm.get ("Data file name topography");
          data_directory_ice_height = Utilities::expand_ASPECT_SOURCE_DIR(prm.get ("Data directory ice height"));
          data_file_name_ice_height = prm.get ("Data file name ice height");
          write_to_file = prm.get_bool ("Output to file");
          output_interval = prm.get_double ("Time between text output");
          if (this->convert_output_to_years())
            output_interval *= year_in_seconds;
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    template <class Archive>
    void SeaLevel<dim>::serialize (Archive &ar, const unsigned int)
    {
      ar &last_output_time
      & sea_level_offset;
    }


    template <int dim>
    void
    SeaLevel<dim>::save (std::map<std::string, std::string> &status_strings) const
    {
      std::ostringstream os;

      // Serialize into a stringstream. Put the following into a code
      // block of its own to ensure the destruction of the 'oa'
      // archive triggers a flush() on the stringstream so we can
      // query the completed string below.
      {
        aspect::oarchive oa (os);
        oa << (*this);
      }

      status_strings["SeaLevel"] = os.str();
    }


    template <int dim>
    void
    SeaLevel<dim>::load (const std::map<std::string, std::string> &status_strings)
    {
      // see if something was saved
      if (status_strings.find("SeaLevel") != status_strings.end())
        {
          std::istringstream is (status_strings.find("SeaLevel")->second);
          aspect::iarchive ia (is);
          ia >> (*this);
        }
    }
  }
}



// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(SeaLevel,
                                  "sea level",
                                  "A postprocessor that computes the sea level for glacial isostatic adjustment"
                                  "modeling. When ice melts and enters the ocean, the ocean water needs to be"
                                  "redistributed in a gravitationally consistent way. With the updated surface"
                                  "loading (ocean and ice) the free surface deformation needs to be computed"
                                  "iteratively before moving to the next time step. "
                                  "A postprocessor intended for use with a deforming top surface. After every step "
                                  "it computes the sea level based on the topography, ocean basin, ice melt, "
                                  "perturbed gravitational potential of the Earth model and gravitational potential "
                                  "of the ice load, relative to a reference datum (initial "
                                  "radius for a spherical shell geometry model). "
                                  "The sea level computation is based on \\cite{Martinec2018}. "
                                  "If 'SeaLevel.Output to file' is set to true, also outputs sea level "
                                  "into text files named `sea_level.NNNNN' in the output directory, "
                                  "where NNNNN is the number of the time step. "
                                  "\n\n"
                                  "The file format then consists of lines with Euclidean coordinates "
                                  "followed by the corresponding sea level value. "
                                  "Sea level is printed/written in meters. ")
  }
}
