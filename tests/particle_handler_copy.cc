/*
  Copyright (C) 2020 by the authors of the ASPECT code.

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

#include <aspect/simulator_signals.h>
#include <aspect/particle/world.h>

#include <deal.II/particles/particle_handler.h>

#include <iostream>

using namespace aspect;

template <int dim>
void duplicate_particle_handler(const SimulatorAccess<dim> &simulator_access,
                                const bool,
                                const unsigned int,
                                const SolverControl &)
{
  const auto &particle_world = simulator_access.get_particle_world();
  dealii::Particles::ParticleHandler<dim> particle_handler;
  std::cout << "duplicating particle handler" << std::endl;

  particle_world.copy_particle_handler(particle_world.get_particle_handler(),
                                       particle_handler);

  auto copied_particle = particle_handler.begin();
  for (const auto &particle: particle_world.get_particle_handler())
    {
      std::cout << "Original position: " << particle.get_location()
                << ". Copied position: " << copied_particle->get_location() << std::endl;
      std::cout << "Original properties: " << particle.get_properties()[0]
                << ". Copied properties: " << copied_particle->get_properties()[0] << std::endl;
      ++copied_particle;
    }

  particle_handler.begin()->get_properties()[0] = 3.14;

  std::cout << "resetting changed particle handler" << std::endl;

  auto &non_const_particle_handler = const_cast<dealii::Particles::ParticleHandler<dim> &>(particle_world.get_particle_handler());

  particle_world.copy_particle_handler(particle_handler,
                                       non_const_particle_handler);

  copied_particle = particle_handler.begin();
  for (const auto &particle: particle_world.get_particle_handler())
    {
      std::cout << "Original position: " << particle.get_location()
                << ". Copied position: " << copied_particle->get_location() << std::endl;
      std::cout << "Original properties: " << particle.get_properties()[0]
                << ". Copied properties: " << copied_particle->get_properties()[0] << std::endl;
      ++copied_particle;
    }
}


template <int dim>
void signal_connector (SimulatorSignals<dim> &signals)
{
  std::cout << "Connecting signals" << std::endl;
  signals.post_advection_solver.connect (&duplicate_particle_handler<dim>);
}


ASPECT_REGISTER_SIGNALS_CONNECTOR(signal_connector<2>,
                                  signal_connector<3>)

