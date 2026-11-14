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

#ifndef _aspect_particle_population_presence_of_melt_h
#define _aspect_particle_population_presence_of_melt_h

#include <aspect/particle/population/interface.h>

namespace aspect
{
  namespace Particle
  {
    namespace Population
    {

      /**
       * A particle population model that prescribes a (constant) minimum
       * and maximum number of particles in cells where melt is present,
       * indicated by the porosity exceeding a threshold. Cells without melt
       * have a minimum and maximum particle count of zero.
       *
       * @ingroup ParticlePopulations
       */
      template <int dim>
      class PresenceOfMelt : public Interface<dim>
      {
        public:
          /**
           * Initialize the plugin and determine if the compositional field
           * porosity exists.
           */
          void
          initialize() override;

          /**
           * @copydoc aspect::Particle::Population::Interface::particle_count_bounds()
           */
          std::pair<unsigned int, unsigned int>
          particle_count_bounds(const typename DoFHandler<dim>::active_cell_iterator &cell) const override;

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
           * Threshold above which a cell is considered to contain melt.
           */
          double porosity_threshold;
      };
    }
  }
}

#endif
