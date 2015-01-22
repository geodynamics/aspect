/*
 Copyright (C) 2011, 2012 by the authors of the ASPECT code.

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

#ifndef __aspect__particle_generator_h
#define __aspect__particle_generator_h

#include <aspect/particle/particle.h>
#include <aspect/particle/world.h>


namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      /**
       * Abstract base class used for classes that generate particles
       */
      template <int dim, class T>
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
           * @param [in] world The particle world the particles will exist in
           * @param [in] total_num_particles Total number of particles to
           * generate. The actual number of generated particles may differ,
           * for example if the generator reads particles from a file this
           * parameter may be ignored.
           */
          virtual
          void
          generate_particles(Particle::World<dim, T> &world,
                             const double total_num_particles) = 0;
      };


      /**
       * Create a generator object.
       *
       * @param[in] generator_type Name of the type of generator to create
       * @return pointer to the generator. Caller needs to delete this
       * pointer.
       */
      template <int dim, class T>
      Interface<dim, T> *
      create_generator_object (const std::string &generator_type);


      /**
       * Return a list of names (separated by '|') of possible particle
       * generators.
       */
      std::string
      generator_object_names ();

    }
  }
}

#endif
