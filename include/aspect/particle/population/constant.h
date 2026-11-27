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

#ifndef _aspect_particle_population_constant_h
#define _aspect_particle_population_constant_h

#include <aspect/particle/population/interface.h>

namespace aspect
{
  namespace Particle
  {
    namespace Population
    {

      /**
       * A particle population model that prescribes the same minimum
       * and maximum numbers of particles in each cell.
       *
       * @ingroup ParticlePopulations
       */
      template <int dim>
      class Constant : public Interface<dim>
      {
        public:
          /**
           * @copydoc aspect::Particle::Population::Interface::particle_count_bounds()
           */
          std::pair<unsigned int, unsigned int>
          particle_count_bounds(const typename DoFHandler<dim>::active_cell_iterator &cell) const override;
      };
    }
  }
}

#endif
