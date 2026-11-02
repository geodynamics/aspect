// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2015-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#ifndef _aspect_particle_generator_ascii_file_h
#define _aspect_particle_generator_ascii_file_h

#include <aspect/particle/generator/interface.h>

namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      /**
       * Generate a distribution of particles that is determined by the
       * coordinates given in an ascii data file.
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
           * @param [in,out] particle_handler The particle handler into which
           * the generated particles should be inserted.
           */
          void
          generate_particles(Particles::ParticleHandler<dim> &particle_handler) override;

          /**
           * Declare the parameters this class takes through input files.
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Read the parameters this class declares from the parameter file.
           */
          void
          parse_parameters (ParameterHandler &prm) override;

        private:
          /**
           * This variable is read from the parameter file through a parameter
           * called 'Data directory'.
           */
          std::string data_directory;

          /**
           * This variable is read from the parameter file through a parameter
           * called 'Data file name'.
           */
          std::string data_filename;
      };

    }
  }
}

#endif
