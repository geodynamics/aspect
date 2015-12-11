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

#ifndef __aspect__particle_generator_ascii_file_h
#define __aspect__particle_generator_ascii_file_h

#include <aspect/particle/generator/interface.h>

namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      /**
       *  Generate a distribution of particles that is determined by the
       *  coordinates given in an ascii data file.
       *
       * @ingroup ParticleGenerators
       */
      template <int dim>
      class AsciiFile : public Interface<dim>
      {
        public:
          /**
           * Reads in a file and generate a set of particles at the prescribed
           * positions.
           *
           * @param [in,out] particles A multimap between cells and their
           * particles. This map will be filled in this function.
           */
          virtual
          void
          generate_particles(std::multimap<types::LevelInd, Particle<dim> > &particles);

          /**
           * Declare the parameters this class takes through input files.
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Read the parameters this class declares from the parameter file.
           */
          virtual
          void
          parse_parameters (ParameterHandler &prm);

        private:
          std::string data_directory;
          std::string data_filename;
      };

    }
  }
}

#endif
