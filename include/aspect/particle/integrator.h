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

#ifndef __aspect__particle_integrator_h
#define __aspect__particle_integrator_h

#include <aspect/particle/world.h>
#include <deal.II/numerics/fe_field_function.h>

namespace aspect
{
  namespace Particle
  {
    namespace Integrator
    {
      /**
       * An abstract class defining virtual methods for performing integration
       * of particle paths through the simulation velocity field.
       */
      template <int dim, class T>
      class Interface
      {
        public:
          /**
           * Constructor.
           */
          Interface(void) {}

          /**
           * Destructor. Made virtual so that derived classes can be created
           * and destroyed through pointers to the base class.
           */
          virtual ~Interface () {}

          /**
           * Perform an integration step of moving the particles by the
           * specified timestep dt. Implementations of this function must
           * update the particle location. If the integrator requires multiple
           * internal steps, this function must return true until all internal
           * steps are finished. Between calls to this function the velocity
           * at the updated particle positions is evaluated and passed to
           * integrate_step during the next call.
           *
           * @param [in,out] world The world to integrate particles in. The
           * particle positions will be changed in this function based on the
           * integration scheme.
           * @param [in] dt The timestep length to perform the integration.
           * @return Whether this function needs to be called again (true) for
           * additional integration steps or if all internal steps are
           * complete (false).
           */
          virtual bool integrate_step(Particle::World<dim, T> *world, const double dt) = 0;

          /**
           * Secify the MPI types and data sizes involved in transferring
           * integration related information between processes. If the
           * integrator samples velocities at different locations and the
           * particle moves between processes during the integration step, the
           * sampled velocities must be transferred with the particle.
           *
           * @param [in,out] data_info Adds MPI data info to the specified
           * vector indicating the quantity and type of values the integrator
           * needs saved for this particle.
           */
          virtual void add_mpi_types(std::vector<MPIDataInfo> &data_info) = 0;

          /**
           * Return data length of the integration related data required for
           * communication in terms of number of doubles.
           *
           * @return The number of doubles required to store the relevant
           * integrator data.
           */
          virtual unsigned int data_len() const = 0;

          /**
           * Read integration related data for a particle specified by id_num
           * from the data vector.
           *
           * @param [in] data The vector of double data to read from.
           * @param [in] pos The position in the data vector to start reading
           * from.
           * @param [in] id_num The id number of the particle to read the data
           * for.
           * @return The position in the vector of the next unread double.
           */
          virtual unsigned int read_data(const std::vector<double> &data, const unsigned int &pos, const double &id_num) = 0;

          /**
           * Write integration related data to a vector for a particle
           * specified by id_num.
           *
           * @param [in,out] data The vector of doubles to write integrator
           * data into.
           * @param [in] id_num The id number of the particle to read the data
           * for.
           */
          virtual void write_data(std::vector<double> &data, const double &id_num) const = 0;
      };


      /**
       * Create an integrator object.
       *
       * @param[in] integrator_name Name of the type of integrator.
       * @return Pointer to instantiated generator object
       */
      template <int dim, class T>
      Interface<dim, T> *
      create_integrator_object (const std::string &integrator_name);


      /**
       * Return a list of names (separated by '|') of possible integrator
       * classes for particles.
       */
      std::string
      integrator_object_names ();

    }
  }
}

#endif
