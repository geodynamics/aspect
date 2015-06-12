/*
 Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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

#ifndef __aspect__particle_generator_interface_h
#define __aspect__particle_generator_interface_h

#include <aspect/particle/world.h>
#include <aspect/particle/particle.h>
#include <aspect/plugins.h>

#include <deal.II/base/parameter_handler.h>

namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      /**
       * Abstract base class used for classes that generate particles
       */
      template <int dim>
      class Interface
      {
        public:
          /**
           * Constructor.
           */
          Interface() {}

          /**
           * Destructor. Made virtual so that derived classes can be created
           * and destroyed through pointers to the base class.
           */
          virtual ~Interface () {}

          /**
           * Generate a specified number of particles in the specified world
           * using the type of generation function implemented by this
           * Generator.
           *
           * @param [in] total_num_particles Total number of particles to
           * generate. The actual number of generated particles may differ,
           * for example if the generator reads particles from a file this
           * parameter may be ignored.
           * @param [inout] world The particle world the particles will exist in
           *
           */
          virtual
          void
          generate_particles(World<dim> &world) = 0;


          /**
           * Declare the parameters this class takes through input files. The
           * default implementation of this function does not describe any
           * parameters. Consequently, derived classes do not have to overload
           * this function if they do not take any runtime parameters.
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Read the parameters this class declares from the parameter file.
           * The default implementation of this function does not read any
           * parameters. Consequently, derived classes do not have to overload
           * this function if they do not take any runtime parameters.
           */
          virtual
          void
          parse_parameters (ParameterHandler &prm);
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
                                   Interface<dim> *(*factory_function) ());

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
      Interface<dim> *
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
    aspect::internal::Plugins::RegisterHelper<aspect::Particle::Generator::Interface<2>,classname<2 > > \
    dummy_ ## classname ## _2d (&aspect::Particle::Generator::register_particle_generator<2>, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::Particle::Generator::Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&aspect::Particle::Generator::register_particle_generator<3>, \
                                name, description); \
  }
    }
  }
}

#endif
