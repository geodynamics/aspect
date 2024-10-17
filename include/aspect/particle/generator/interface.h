/*
 Copyright (C) 2015 - 2024 by the authors of the ASPECT code.

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

#ifndef _aspect_particle_generator_interface_h
#define _aspect_particle_generator_interface_h

#include <aspect/particle/interface.h>
#include <aspect/simulator_access.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/generators.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/exceptions.h>

#include <random>

#include <map>

namespace aspect
{
  namespace Particle
  {
    /**
     * A namespace in which we define everything that has to do with defining
     * the particle generation.
     *
     * @ingroup ParticleGenerators
     */
    namespace Generator
    {
      using namespace dealii::Particles;
      using dealii::Particles::Particle;

      /**
       * Exception denoting a division by zero.
       */
      DeclExceptionMsg (ExcParticlePointNotInDomain,
                        "You requested to generate a particle at a position that "
                        "is not owned by this process, therefore the "
                        "Particle::Generator::Interface::generate_particle() function "
                        "refused to create it. You can circumvent this error message "
                        "by catching the ExcParticlePointNotInDomain exception and "
                        "do whatever you think is appropriate in this case.");

      /**
       * Abstract base class used for classes that generate particles.
       *
       * @ingroup ParticleGenerators
       */
      template <int dim>
      class Interface : public SimulatorAccess<dim>, public ParticleInterfaceBase
      {
        public:
          virtual
          void
          initialize () override;

          /**
           * Generate particles. Every derived class
           * has to decide on the method and number of particles to generate,
           * for example using input parameters declared in their
           * declare_parameters and parse_parameters functions. This function
           * should generate the particles and insert them into @p particle_handler
           * by calling its respective functions.
           *
           * @param [in,out] particle_handler The particle handler into which
           * the generated particles should be inserted.
           */
          virtual
          void
          generate_particles(Particles::ParticleHandler<dim> &particle_handler) = 0;

          /**
           * Generate one particle in the given cell. This function's main purpose
           * is to provide functionality to fill up cells with too few particles
           * after refinement. Of course it can also be utilized by derived classes
           * to generate the initial particle distribution.
           */
          std::pair<Particles::internal::LevelInd,Particle<dim>>
          generate_particle (const typename parallel::distributed::Triangulation<dim>::active_cell_iterator &cell,
                             const types::particle_index id);

        protected:
          /**
           * Generate a particle at the specified position and with the
           * specified id. Many derived classes use this functionality,
           * therefore it is implemented here to avoid duplication.
           * In case the position is not in the local domain this function
           * throws an exception of type ExcParticlePointNotInDomain, which
           * can be caught in the calling plugin.
           *
           * @deprecated: This function uses an old return type and is deprecated.
           * Use the insert_particle_at_position() function below instead.
           */
          DEAL_II_DEPRECATED
          std::pair<Particles::internal::LevelInd,Particle<dim>>
          generate_particle(const Point<dim> &position,
                            const types::particle_index id) const;

          /**
           * Generate a particle at the specified position and with the
           * specified id and insert it into the @p particle_handler.
           * Many derived classes use this functionality,
           * therefore it is implemented here to avoid duplication.
           * In case the position is not in the local domain this function
           * throws an exception of type ExcParticlePointNotInDomain, which
           * can be caught in the calling plugin. Note that since the cell
           * in which the particle is generated is not known, it has to be
           * found, which is an expensive operation.
           *
           * @param position Position of the particle.
           * @param id The id of the particle.
           * @param particle_handler The particle handler into which the particle
           * should be inserted.
           *
           * @return An iterator to the inserted particle. If the position was
           * not in the local domain, an iterator to particle_handler.end() is
           * returned.
           */
          Particles::ParticleIterator<dim>
          insert_particle_at_position(const Point<dim> &position,
                                      const types::particle_index id,
                                      Particles::ParticleHandler<dim> &particle_handler) const;

          /**
           * Random number generator. For reproducibility of tests it is
           * initialized in the constructor with a constant.
           */
          std::mt19937            random_number_generator;
      };

      /**
       * Register a particle generator so that it can be selected from
       * the parameter file.
       *
       * @param name A string that identifies the particle generator
       * @param description A text description of what this generator does and that
       * will be listed in the documentation of the parameter file.
       * @param declare_parameters_function A pointer to a function that can be
       * used to declare the parameters that this particle generator wants to read
       * from input files.
       * @param factory_function A pointer to a function that can create an
       * object of this particle generator.
       *
       * @ingroup ParticleGenerators
       */
      template <int dim>
      void
      register_particle_generator (const std::string &name,
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
       * @ingroup ParticleGenerators
       */
      template <int dim>
      std::unique_ptr<Interface<dim>>
      create_particle_generator (ParameterHandler &prm);

      /**
       * Declare the runtime parameters of the registered particle generators.
       *
       * @ingroup ParticleGenerators
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
       * for a particle generator, register it with the functions that
       * can declare their parameters and create these objects.
       *
       * @ingroup ParticleGenerators
       */
#define ASPECT_REGISTER_PARTICLE_GENERATOR(classname, name, description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_PARTICLE_GENERATOR_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::Particle::Generator::Interface<2>,classname<2 >> \
    dummy_ ## classname ## _2d (&aspect::Particle::Generator::register_particle_generator<2>, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::Particle::Generator::Interface<3>,classname<3>> \
    dummy_ ## classname ## _3d (&aspect::Particle::Generator::register_particle_generator<3>, \
                                name, description); \
  }
    }
  }
}

#endif
