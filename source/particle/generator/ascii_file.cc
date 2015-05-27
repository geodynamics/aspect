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

#include <aspect/particle/generator/ascii_file.h>

#include <deal.II/grid/grid_tools.h>


namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      // Generate random uniform distribution of particles over entire simulation domain

          /**
           * Constructor.
           *
           * @param[in] The MPI communicator for synchronizing particle generation.
           */
        template <int dim>
        AsciiFile<dim>::AsciiFile() {}

          /**
           * TODO: Update comments
           * Generate a uniformly randomly distributed set of particles in the current triangulation.
           */
          // TODO: fix the particle system so it works even with processors assigned 0 cells
        template <int dim>
        void
        AsciiFile<dim>::generate_particles(Particle::World<dim> &world,
                                            const double total_num_particles)
          {
          const std::string filename = data_directory+data_filename;
          std::ifstream in(filename.c_str(), std::ios::in);
          AssertThrow (in,
                       ExcMessage (std::string("Couldn't open data file <"
                                               +
                                               filename
                                               +
                                               ">.")));

          // Read header lines and table size
          while (in.peek() == '#')
            {
            }

          char sep;
          // Read data lines
          unsigned int id = 0;
          Point<dim> coordinates;

          if (dim == 2)
            while (in >> coordinates[0] >> sep >> coordinates[1])
              {
                generate_particle(world,coordinates,id);
                id++;
              }
          else if (dim == 3)
            while (in >> coordinates[0] >> sep >> coordinates[1] >> sep >> coordinates[2])
              {
                generate_particle(world,coordinates,id);
                id++;
              }
          else
            Assert(false,ExcNotImplemented());
          }



        template <int dim>
        void
        AsciiFile<dim>::generate_particle(Particle::World<dim> &world,
                                           const Point<dim> &position,
                                           const unsigned int id)
            {
              typename parallel::distributed::Triangulation<dim>::active_cell_iterator it =
                  (GridTools::find_active_cell_around_point<> (this->get_mapping(), this->get_triangulation(), position)).first;;

              if (it->is_locally_owned())
                {
                  //Only try to add the point if the cell it is in, is on this processor
                  BaseParticle<dim> new_particle(position, id);
                  world.add_particle(new_particle, std::make_pair(it->level(), it->index()));
                }
            }


        template <int dim>
        void
        AsciiFile<dim>::declare_parameters (ParameterHandler &prm)
        {
          prm.enter_subsection("Postprocess");
          {
            prm.enter_subsection("Tracers");
            {
              prm.enter_subsection("Generator");
              {
                prm.enter_subsection("Ascii file");
                {
                  prm.declare_entry ("Data directory",
                                     "$ASPECT_SOURCE_DIR/data/particle/generator/ascii/",
                                     Patterns::DirectoryName (),
                                     "The name of a directory that contains the tracer data. This path "
                                     "may either be absolute (if starting with a '/') or relative to "
                                     "the current directory. The path may also include the special "
                                     "text '$ASPECT_SOURCE_DIR' which will be interpreted as the path "
                                     "in which the ASPECT source files were located when ASPECT was "
                                     "compiled. This interpretation allows, for example, to reference "
                                     "files located in the 'data/' subdirectory of ASPECT. ");
                  prm.declare_entry ("Data file name", "tracer.dat",
                                     Patterns::Anything (),
                                     "The name of the tracer file.");
                prm.leave_subsection();
              }
              prm.leave_subsection();
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        }


        template <int dim>
        void
        AsciiFile<dim>::parse_parameters (ParameterHandler &prm)
        {
          prm.enter_subsection("Postprocess");
          {
            prm.enter_subsection("Tracers");
            {
              prm.enter_subsection("Generator");
              {
                prm.enter_subsection("Ascii file");
                {
                  // Get the path to the data files. If it contains a reference
                  // to $ASPECT_SOURCE_DIR, replace it by what CMake has given us
                  // as a #define
                  data_directory        = prm.get ("Data directory");
                  {
                    const std::string      subst_text = "$ASPECT_SOURCE_DIR";
                    std::string::size_type position;
                    while (position = data_directory.find (subst_text),  position!=std::string::npos)
                      data_directory.replace (data_directory.begin()+position,
                                              data_directory.begin()+position+subst_text.size(),
                                              ASPECT_SOURCE_DIR);
                  }

                  data_filename    = prm.get ("Data file name");
                }
                prm.leave_subsection();
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
  namespace Particle
  {
    namespace Generator
    {
    ASPECT_REGISTER_PARTICLE_GENERATOR(AsciiFile,
                                       "ascii file",
                                       "Generates a distribution of tracers from coordinates "
                                       "specified in an Ascii data file. The file format is "
                                       "identical to the one from other ascii data plugins, "
                                       "except that it only needs to contain as many columns as "
                                       "dimensions to specify the position of each tracer. Note "
                                       "that this plugin does ignore the number of tracers set in "
                                       "the input file.")
    }
  }
}

