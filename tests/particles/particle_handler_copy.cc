// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2020-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#include <aspect/simulator_signals.h>

// On purpose use the old file name here to test that backwards compatibility works.
#include <aspect/particle/world.h>

#include <iostream>

using namespace aspect;

template <int dim>
void duplicate_particle_handler(const SimulatorAccess<dim> &simulator_access,
                                const bool,
                                const unsigned int,
                                const SolverControl &)
{
  // On purpose use the old class name here to test that backwards compatibility works.
  const Particle::World<dim> &particle_manager = simulator_access.get_particle_manager(0);
  dealii::Particles::ParticleHandler<dim> particle_handler;
  std::cout << "duplicating particle handler" << std::endl;

  particle_manager.copy_particle_handler(particle_manager.get_particle_handler(),
                                         particle_handler);

  auto copied_particle = particle_handler.begin();
  for (const auto &particle: particle_manager.get_particle_handler())
    {
      std::cout << "Original position: " << particle.get_location()
                << ". Copied position: " << copied_particle->get_location() << std::endl;
      std::cout << "Original properties: " << particle.get_properties()[0]
                << ". Copied properties: " << copied_particle->get_properties()[0] << std::endl;
      ++copied_particle;
    }

  particle_handler.begin()->get_properties()[0] = 3.14;

  std::cout << "resetting changed particle handler" << std::endl;

  auto &non_const_particle_handler = const_cast<dealii::Particles::ParticleHandler<dim> &>(particle_manager.get_particle_handler());

  particle_manager.copy_particle_handler(particle_handler,
                                         non_const_particle_handler);

  copied_particle = particle_handler.begin();
  for (const auto &particle: particle_manager.get_particle_handler())
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
