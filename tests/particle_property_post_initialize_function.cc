/*
  Copyright (C) 2022 - 2024 by the authors of the ASPECT code.

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
#include <aspect/simulator_access.h>
#include <aspect/particle/property/interface.h>

#include <iostream>

unsigned int counter_without = 0, counter_with = 0;
bool quiet = true;


namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      template <int dim>
      class PostInitializeParticleProperty : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
      {
        public:

          virtual
          void
          initialize ()
          {
            std::cout << "initialize" << std::endl;
            const Particle::Property::Manager<dim> &manager = this->get_particle_manager(this->get_particle_manager_index()).get_property_manager();
            post_initialized_info = manager.get_data_info().get_field_index_by_name("initial position");
            std::cout << "initial position: post_initialized_info = " << post_initialized_info << std::endl;

            post_initialized_info = manager.get_data_info().get_field_index_by_name("initial C_1");
            std::cout << "initial composition C_1: post_initialized_info = " << post_initialized_info << std::endl;

            post_initialized_info = manager.get_data_info().get_field_index_by_name("velocity");
            std::cout << "velocity: post_initialized_info = " << post_initialized_info << std::endl;
            exit(0);
          }

          virtual
          std::vector<std::pair<std::string, unsigned int>>
          get_property_information() const
          {
            return std::vector<std::pair<std::string, unsigned int>>();
          }

        private:
          unsigned int post_initialized_info;

      };
    }
  }
}


namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      ASPECT_REGISTER_PARTICLE_PROPERTY(PostInitializeParticleProperty,
                                        "PostInitializeParticleProperty",
                                        "")
    }
  }
}
