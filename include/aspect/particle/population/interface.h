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

#ifndef _aspect_particle_population_interface_h
#define _aspect_particle_population_interface_h

#include <aspect/particle/interface.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/parameter_handler.h>

namespace aspect
{
  namespace Particle
  {
    /**
     * A namespace in which we define everything that has to do with defining
     * the particle population.
     *
     * @ingroup ParticlePopulations
     */
    namespace Population
    {
      /**
       * Abstract base class for plugins that determine a minimum and maximum
       * particle count in each cell.
       *
       * @ingroup ParticlePopulations
       */
      template <int dim>
      class Interface : public SimulatorAccess<dim>, public ParticleInterfaceBase
      {
        public:
          /**
           * Return the minimum and maximum number of particles allowed in a
           * given cell. The first value of the returned pair is the minimum
           * and the second is the maximum. Every derived class has to determine
           * these numbers, for example using input parameters declared in their
           * declare_parameters and parse_parameters functions.
           *
           * @param cell the target cell in which to determine the particle count bounds.
           */
          virtual
          std::pair<unsigned int, unsigned int>
          particle_count_bounds(const typename DoFHandler<dim>::active_cell_iterator &cell) const = 0;
      };

      /**
       * Register a particle population so that it can be selected from
       * the parameter file.
       *
       * @param name A string that identifies the particle population
       * @param description A text description of what this population does and that
       * will be listed in the documentation of the parameter file.
       * @param declare_parameters_function A pointer to a function that can be
       * used to declare the parameters that this particle population wants to read
       * from input files.
       * @param factory_function A pointer to a function that can create an
       * object of this particle population.
       *
       * @ingroup ParticlePopulations
       */
      template <int dim>
      void
      register_particle_population (const std::string &name,
                                    const std::string &description,
                                    void (*declare_parameters_function) (ParameterHandler &),
                                    std::unique_ptr<Interface<dim>> (*factory_function) ());

      /**
       * A function that given the name of a model returns a pointer to an
       * object that describes it. Ownership of the pointer is transferred to
       * the caller.
       *
       * The model object returned is not yet initialized and has not
       * read its runtime parameters yet.
       *
       * @ingroup ParticlePopulations
       */
      template <int dim>
      std::unique_ptr<Interface<dim>>
      create_particle_population (ParameterHandler &prm);

      /**
       * Declare the runtime parameters of the registered particle populations.
       *
       * @ingroup ParticlePopulations
       */
      template <int dim>
      void
      declare_parameters (ParameterHandler &prm);


      /**
       * For the current plugin subsystem, write a connection graph of all of the
       * plugins we know about, in the format that the
       * programs dot and neato understand. This allows for a visualization of
       * how all of the plugins that ASPECT knows about are interconnected, and
       * connect to other parts of the ASPECT code.
       *
       * @param output_stream The stream to write the output to.
       */
      template <int dim>
      void
      write_plugin_graph (std::ostream &output_stream);


      /**
       * Given a class name, a name, and a description for the parameter file
       * for a particle population, register it with the functions that
       * can declare their parameters and create these objects.
       *
       * @ingroup ParticlePopulations
       */
#define ASPECT_REGISTER_PARTICLE_POPULATION(classname, name, description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_PARTICLE_POPULATION_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::Particle::Population::Interface<2>,classname<2 >> \
    dummy_ ## classname ## _2d (&aspect::Particle::Population::register_particle_population<2>, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::Particle::Population::Interface<3>,classname<3>> \
    dummy_ ## classname ## _3d (&aspect::Particle::Population::register_particle_population<3>, \
                                name, description); \
  }
    }
  }
}

#endif
