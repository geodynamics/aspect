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


#include <aspect/postprocess/sealevel.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/simulator.h>
#include <aspect/global.h>
#include <aspect/postprocess/topography.h>
#include <aspect/postprocess/geoid.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

#include <boost/lexical_cast.hpp>

#include <cmath>
#include <limits>

// TODO:
// Get geoid postprocessor to work with a free surface
// Obtain geoid anomaly of model from the geoid postprocessor
// Obtain geoid anomaly of surface load (from new postprocessor, from geoid postprocessor when implemented in there, or compute here)
// Load an ascii file with the initial ocean basin (equal for all time steps, fixed ocean basin for now)
// Compute ice mass change: follows from ice input data (does not change over iterations)
// Compute the sea level change (from ice mass change, surface displacement, change in gravitational potential from model and ice load)
// Compute new surface gravity, dependent on surface topography and Earth model
// Make new ocean loading parameter: sea level * sea water density * new surface gravity

namespace aspect
{
  // namespace InitialTemperature
  // {
  //   template <int dim>
  //   TOPOGRAPHY<dim>::TOPOGRAPHY ()
  //   {}
    
  //   template <int dim>
  //   void
  //   TOPOGRAPHY<dim>::initialize ()
  //   {
  //     Utilities::AsciiDataInitial<dim>::initialize(1);
  //   }

  //   template <int dim>
  //   void
  //   TOPOGRAPHY<dim>::declare_parameters (ParameterHandler &prm)
  //   {
  //     prm.enter_subsection ("Initial temperature model");
  //     {
  //       Utilities::AsciiDataBase<dim>::declare_parameters(prm,
  //                                                         "$ASPECT_SOURCE_DIR/../input_files/",
  //                                                         "topo.txt");
  //     }
  //     prm.leave_subsection();
  //   }

  //   template <int dim>
  //   void
  //   TOPOGRAPHY<dim>::parse_parameters (ParameterHandler &prm)
  //   {
  //     prm.enter_subsection ("Initial temperature model");
  //     {
  //       Utilities::AsciiDataBase<dim>::parse_parameters(prm);
  //     }
  //     prm.leave_subsection();
  //   }

  // }
  namespace Postprocess
  {
    template <int dim>
    double
    SeaLevel<dim>::sea_level_equation(const Point<dim> &position)
    {
      double c = 0.0;
      double gravity_perturbation = 0.0;
      const double elevation = this->get_geometry_model().height_above_reference_surface(position);
      int ocean_mask = 1;
      double sea_level = (gravity_perturbation - elevation + c)*ocean_mask;
      return sea_level;
    }
    
    template <int dim>
    void
    SeaLevel<dim>::initialize ()
    {
      Utilities::AsciiDataInitial<dim>::initialize(1);
    }

    template <int dim>
    std::pair<std::string,std::string>
    SeaLevel<dim>::execute (TableHandler &statistics)
    {
      //Disallow use of the plugin for any other than a 3D spherical shell geometry
      AssertThrow (Plugins::plugin_type_matches<const GeometryModel::SphericalShell<dim>>(this->get_geometry_model())
                   &&
                   dim == 3,
                   ExcMessage("The sea level postprocessor is currently only implemented for the 3D spherical shell geometry model."));

      const types::boundary_id relevant_boundary = this->get_geometry_model().translate_symbolic_boundary_name_to_id ("top");

//      template <int dim>
//      QGauss<dim>::QGauss(const unsigned int n)
//        : Quadrature<dim>(QGauss<dim - 1>(n), QGauss<1>(n))
//      {}
//
//      template <int dim>
//      QTrapez<dim>::QTrapez()
//        : Quadrature<dim>(QTrapez<dim - 1>(), QTrapez<1>())
//      {}

      // Get a quadrature rule that exists for the Gaussian integration points
      const unsigned int quadrature_degree = this->get_fe().base_element(this->introspection().base_elements.velocities).degree+1;
      QGauss<dim-1> face_points(quadrature_degree);
      FEFaceValues<dim> face_vals (this->get_mapping(), this->get_fe(), face_points, update_quadrature_points);

      // have a stream into which we write the data. the text stream is then
      // later sent to processor 0
      std::ostringstream output_stats;
      std::ostringstream output_file;

      // Choose stupidly large values for initialization
      double local_max_height = -std::numeric_limits<double>::max();
      double local_min_height = std::numeric_limits<double>::max();

      // loop over all of the surface cells and save the elevation to stored_value
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned() && cell->at_boundary())
          for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
            if (cell->face(face_no)->at_boundary())
              {
                if ( cell->face(face_no)->boundary_id() != relevant_boundary)
                  continue;

                face_vals.reinit( cell, face_no);

                for (unsigned int gausspoint = 0; gausspoint < face_points.size(); ++gausspoint)
                  {
                    const Point<dim> vertex = face_vals.quadrature_point(gausspoint);
                    const double elevation = sea_level_equation(vertex);
                    if (write_to_file)
                      output_file << vertex << ' '<< elevation << std::endl;
                    if ( elevation > local_max_height)
                      local_max_height = elevation;
                    if ( elevation < local_min_height)
                      local_min_height = elevation;
                  }
              }

      //Calculate min/max sea level across all processes
      const double max_sealevel = Utilities::MPI::max(local_max_height, this->get_mpi_communicator());
      const double min_sealevel = Utilities::MPI::min(local_min_height, this->get_mpi_communicator());

      //Write results to statistics file
      statistics.add_value ("Minimum sea level (m)",
                            min_sealevel);
      statistics.add_value ("Maximum sea level (m)",
                            max_sealevel);
      const char *columns[] = { "Minimum sealevel (m)",
                                "Maximum sealevel (m)"
                              };
      for (unsigned int i=0; i<sizeof(columns)/sizeof(columns[0]); ++i)
        {
          statistics.set_precision (columns[i], 8);
          statistics.set_scientific (columns[i], true);
        }

      output_stats.precision(4);
      output_stats << min_sealevel << " m, "
                   << max_sealevel << " m";

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
        return std::pair<std::string, std::string> ("Sea level min/max:",
                                                    output_stats.str());

      std::string filename = this->get_output_directory() +
                             "sealevel." +
                             Utilities::int_to_string(this->get_timestep_number(), 5);
      if (this->get_parameters().run_postprocessors_on_nonlinear_iterations)
        filename.append("." + Utilities::int_to_string (this->get_nonlinear_iteration(), 4));

      const unsigned int max_data_length = Utilities::MPI::max (output_file.str().size()+1,
                                                                this->get_mpi_communicator());

      const unsigned int mpi_tag = 777;

      // on processor 0, collect all of the data the individual processors send
      // and concatenate them into one file
      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {
          std::ofstream file (filename.c_str());

          file << "# "
               << "x y z"
               << " sealevel" << std::endl;

          // first write out the data we have created locally
          file << output_file.str();

          std::string tmp;
          tmp.resize (max_data_length, '\0');

          // then loop through all of the other processors and collect
          // data, then write it to the file
          for (unsigned int p=1; p<Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()); ++p)
            {
              MPI_Status status;
              // get the data. note that MPI says that an MPI_Recv may receive
              // less data than the length specified here. since we have already
              // determined the maximal message length, we use this feature here
              // rather than trying to find out the exact message length with
              // a call to MPI_Probe.
              const int ierr = MPI_Recv (&tmp[0], max_data_length, MPI_CHAR, p, mpi_tag,
                                         this->get_mpi_communicator(), &status);
              AssertThrowMPI(ierr);

              // output the string. note that 'tmp' has length max_data_length,
              // but we only wrote a certain piece of it in the MPI_Recv, ended
              // by a \0 character. write only this part by outputting it as a
              // C string object, rather than as a std::string
              file << tmp.c_str();
            }
        }
      else
        // on other processors, send the data to processor zero. include the \0
        // character at the end of the string
        {
          output_file << "\0";
          const int ierr = MPI_Send (&output_file.str()[0], output_file.str().size()+1, MPI_CHAR, 0, mpi_tag,
                                     this->get_mpi_communicator());
          AssertThrowMPI(ierr);
        }

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

      return std::pair<std::string, std::string> ("Sea level min/max:",
                                                  output_stats.str());
    }


    template <int dim>
    std::list<std::string>
    SeaLevel<dim>::required_other_postprocessors() const
    {
      std::list<std::string> deps;
      //deps.emplace_back("topography");
      deps.emplace_back("geoid");
      return deps;
    }


    template <int dim>
    void
    SeaLevel<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Sea level");
        {
          Utilities::AsciiDataBase<dim>::declare_parameters(prm,
                                                          "$ASPECT_SOURCE_DIR/data/geometry-model/initial-topography-model/ascii-data/test/",
                                                          "shell_3d_outer.0.txt");
          prm.declare_entry ("Output to file", "false",
                             Patterns::List(Patterns::Bool()),
                             "Whether or not to write sea level to a text file named named "
                             "'sealevel.NNNNN' in the output directory");

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



    // template <int dim>
    // void
    // TOPOGRAPHY<dim>::declare_parameters (ParameterHandler &prm)
    // {
    //   prm.enter_subsection ("Initial temperature model");
    //   {
    //     Utilities::AsciiDataBase<dim>::declare_parameters(prm,
    //                                                       "$ASPECT_SOURCE_DIR/../input_files/",
    //                                                       "topo.txt");
    //   }
    //   prm.leave_subsection();
    // }

    // template <int dim>
    // void
    // TOPOGRAPHY<dim>::parse_parameters (ParameterHandler &prm)
    // {
    //   prm.enter_subsection ("Initial temperature model");
    //   {
    //     Utilities::AsciiDataBase<dim>::parse_parameters(prm);
    //   }
    //   prm.leave_subsection();
    // }

  // }
    // template <int dim>
    // void
    // TOPOGRAPHY<dim>::declare_parameters (ParameterHandler &prm)
    // {
    //   prm.enter_subsection ("Initial temperature model");
    //   {
    //     Utilities::AsciiDataBase<dim>::declare_parameters(prm,
    //                                                       "$ASPECT_SOURCE_DIR/../input_files/",
    //                                                       "topo.txt");
    //   }
    //   prm.leave_subsection();
    // }



    template <int dim>
    void
    SeaLevel<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Sea level");
        {
          Utilities::AsciiDataBase<dim>::parse_parameters(prm);
          write_to_file = prm.get_bool ("Output to file");
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
    ASPECT_REGISTER_POSTPROCESSOR(SeaLevel,
                                  "sea level",
                                  "A postprocessor intended for use with a deforming top surface. After every step "
                                  "it computes the sea level based on the topography, ocean basin, ice melt, "
                                  "perturbed gravitational potential of the Earth model and gravitational potential "
                                  "of the ice load, relative to a reference datum (initial "
                                  "radius for a spherical shell geometry model). "
                                  "If 'SeaLevel.Output to file' is set to true, also outputs sea level "
                                  "into text files named `sealevel.NNNNN' in the output directory, "
                                  "where NNNNN is the number of the time step.\n"
                                  "The file format then consists of lines with Euclidean coordinates "
                                  "followed by the corresponding sea level value."
                                  "Sea level is printed/written in meters.")
  }
}
