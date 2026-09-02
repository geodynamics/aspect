/*
 Copyright (C) 2015 - 2026 by the authors of the ASPECT code.

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

#ifndef _aspect_particle_property_particle_generation_time_h
#define _aspect_particle_property_particle_generation_time_h

#include <aspect/particle/property/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      /**
       * A class that initializes particle properties based on the
       * generation time of the particles.
       *
       * @ingroup ParticleProperties
       */
      template <int dim>
      class ParticleGenerationTime : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * @copydoc aspect::Particle::Property::Interface::initialize_one_particle_property()
           */
          void
          initialize_one_particle_property (const Point<dim> &position,
                                            std::vector<double> &particle_properties) const override;

          /**
           * Returns an enum, which determines how this particle property is
           * initialized for particles that are created later than the initial
           * particle generation. For this property the value of
           * generated particles is set to the initialized value.
           */
          InitializationModeForLateParticles
          late_initialization_mode () const override;

          /**
          * @copydoc aspect::Particle::Property::Interface::get_property_information()
          */
          std::vector<std::pair<std::string, unsigned int>>
          get_property_information() const override;
      };
    }
  }
}

#endif
