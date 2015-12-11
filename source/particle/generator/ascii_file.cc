/*
  Copyright (C) 2015 by the authors of the ASPECT code.

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


namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      template <int dim>
      void
      AsciiFile<dim>::generate_particles(std::multimap<types::LevelInd, Particle<dim> > &particles)
      {
        const std::string filename = data_directory+data_filename;
        std::ifstream in(filename.c_str(), std::ios::in);
        AssertThrow (in,
                     ExcMessage (std::string("Could not open data file <"
                                             +
                                             filename
                                             +
                                             ">.")));

        // Skip header lines
        while (in.peek() == '#')
          {
            std::string temp;
            getline(in,temp);
          }

        // Read data lines
        types::particle_index particle_index = 0;
        Point<dim> particle_position;

        while (in >> particle_position)
          {
            // Try to add the particle. If it is not in this domain, do not
            // worry about it and move on to next point.
            try
              {
                particles.insert(this->generate_particle(particle_position,particle_index++));
              }
            catch (ExcParticlePointNotInDomain &)
              {}
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
                                         "a simple text file, with as many columns as spatial "
                                         "dimensions and as many lines as tracers to be generated. "
                                         "Initial comment lines starting with `\\#' will be discarded."
                                         "All of the values that define this generator are read "
                                         "from a section ``Particle generator/Ascii file'' in the "
                                         "input file, see "
                                         "Section~\\ref{parameters:Particle_20generator/Ascii_20file}.")
    }
  }
}

