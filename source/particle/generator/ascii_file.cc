// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2015-2024 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#include <aspect/particle/generator/ascii_file.h>
#include <aspect/utilities.h>


namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      template <int dim>
      void
      AsciiFile<dim>::generate_particles(Particles::ParticleHandler<dim> &particle_handler)
      {
        const std::string filename = data_directory+data_filename;

        // Read data from disk and distribute among processes
        std::istringstream in(Utilities::read_and_distribute_file_content(filename, this->get_mpi_communicator()));

        // Skip header lines
        while (in.peek() == '#')
          {
            std::string temp;
            std::getline(in,temp);
          }

        // Read data lines
        types::particle_index particle_index = 0;
        Point<dim> particle_position;

        while (in >> particle_position)
          {
            this->insert_particle_at_position(particle_position, particle_index, particle_handler);
            ++particle_index;
          }
        particle_handler.update_cached_numbers();
      }


      template <int dim>
      void
      AsciiFile<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Generator");
        {
          prm.enter_subsection("Ascii file");
          {
            prm.declare_entry ("Data directory",
                               "$ASPECT_SOURCE_DIR/data/particle/generator/ascii/",
                               Patterns::DirectoryName (),
                               "The name of a directory that contains the particle data. This path "
                               "may either be absolute (if starting with a '/') or relative to "
                               "the current directory. The path may also include the special "
                               "text '$ASPECT_SOURCE_DIR' which will be interpreted as the path "
                               "in which the ASPECT source files were located when ASPECT was "
                               "compiled. This interpretation allows, for example, to reference "
                               "files located in the `data/' subdirectory of ASPECT. ");
            prm.declare_entry ("Data file name", "particle.dat",
                               Patterns::Anything (),
                               "The name of the particle file.");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }


      template <int dim>
      void
      AsciiFile<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Generator");
        {
          prm.enter_subsection("Ascii file");
          {
            data_directory = Utilities::expand_ASPECT_SOURCE_DIR(prm.get ("Data directory"));

            data_filename    = prm.get ("Data file name");
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
                                         "Generates a distribution of particles from coordinates "
                                         "specified in an Ascii data file. The file format is "
                                         "a simple text file, with as many columns as spatial "
                                         "dimensions and as many lines as particles to be generated. "
                                         "Initial comment lines starting with `#' will be discarded. "
                                         "Note that this plugin always generates as many particles "
                                         "as there are coordinates in the data file, the "
                                         "``Particles/Number of particles'' parameter "
                                         "has no effect on this plugin. "
                                         "All of the values that define this generator are read "
                                         "from a section ``Particles/Generator/Ascii file'' in the "
                                         "input file, see "
                                         "Section~\\ref{parameters:Particles/Generator/Ascii_20file}.")
    }
  }
}
