/*
 Copyright (C) 2015 by the authors of the ASPECT code.

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

#ifndef __aspect__particle_integrator_interface_h
#define __aspect__particle_integrator_interface_h

#include <aspect/particle/particle.h>
#include <aspect/plugins.h>
#include <aspect/global.h>

#include <deal.II/base/parameter_handler.h>

namespace aspect
{
  namespace Particle
  {
    namespace Integrator
    {
      using namespace dealii;

      /**
       * An abstract class defining virtual methods for performing integration
       * of particle paths through the simulation velocity field.
       *
       * @ingroup ParticleIntegrators
       */
      template <int dim>
      class Interface
      {
        public:
          /**
           * Destructor. Made virtual so that derived classes can be created
           * and destroyed through pointers to the base class.
           */
          virtual ~Interface ();

          /**
           * Perform an integration step of moving the particles of one cell
           * by the specified timestep dt. Implementations of this function
           * must update the particle location. Between calls to this function
           * the velocity at the updated particle positions is evaluated and
           * passed as input argument during the next call.
           *
           * @param [in] begin_particle An iterator to the first particle to be moved.
           * @param [in] end_particle An iterator to the last particle to be moved.
           * @param [in] old_velocities The velocities at t_n, i.e. before the
           * particle movement, for all particles between @p begin_particle
           * and @p end_particle at their current position.
           * @param [in] velocities The velocities at the particle positions
           * at t_{n+1}, i.e. after the particle movement. Note that this is
           * the velocity at the old positions, but at the new time. It is the
           * responsibility of this function to compute the new location of
           * the particles.
           * @param [in] dt The length of the integration timestep.
           */
          virtual
          void
          local_integrate_step(const typename std::multimap<types::LevelInd, Particle<dim> >::iterator &begin_particle,
                               const typename std::multimap<types::LevelInd, Particle<dim> >::iterator &end_particle,
                               const std::vector<Tensor<1,dim> > &old_velocities,
                               const std::vector<Tensor<1,dim> > &velocities,
                               const double dt) = 0;

          /**
           * This function is called at the end of every integration step.
           * In case of multi-step integrators it signals the beginning of a
           * new integration step. The default implementation always returns
           * false, which is ok for single-step integration methods.
           *
           * @return This function returns true if the integrator requires
           * another integration step. The particle world will continue
           * to start new integration steps until this function returns false.
           */
          virtual bool new_integration_step();

          /**
           * Return data length of the integration related data required for
           * communication in terms of number of bytes. When data about
           * particles is transported from one processor to another, or stored
           * on disk for snapshots, integrators get the chance to store
           * whatever information they need with each particle. This function
           * returns how many pieces of additional information a concrete
           * integrator class needs to store for each particle.
           *
           * @return The number of bytes required to store the relevant
           * integrator data for one particle.
           */
          virtual unsigned int get_data_size() const;

          /**
           * Read integration related data for a particle specified by particle_id
           * from the data array. This function is called after transferring
           * a particle to the local domain during an integration step.
           *
           * @param [in] data A pointer into the data array. The pointer
           * marks the position where this function starts reading.
           * @param [in] particle_id The id number of the particle to read the data
           * for.
           * @return The updated position of the pointer into the data array.
           * The return value is @p data advanced by get_data_size() bytes.
           */
          virtual
          const void *
          read_data(const void *data,
                    const types::particle_index particle_id);

          /**
           * Write integration related data to a vector for a particle
           * specified by particle_id. This function is called in cases where
           * particles leave the local domain during an integration step to
           * transfer this data to another process.
           *
           * @param [in] data A pointer into the array of integrator data.
           * @param [in] particle_id The id number of the particle to write the data
           * for.
           * @return The updated position of the pointer into the data array.
           * The return value is @p data advanced by get_data_size() bytes.
           */
          virtual
          void *
          write_data(void *data,
                     const types::particle_index particle_id) const;


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
       * Return a list of names (separated by '|') of possible integrator
       * classes for particles.
       */
      std::string
      integrator_object_names ();


      /**
       * Register a particle integrator so that it can be selected from
       * the parameter file.
       *
       * @param name A string that identifies the particle integrator
       * @param description A text description of what this integrator does and that
       * will be listed in the documentation of the parameter file.
       * @param declare_parameters_function A pointer to a function that can be
       * used to declare the parameters that this particle integrator wants to read
       * from input files.
       * @param factory_function A pointer to a function that can create an
       * object of this particle integrator.
       *
       * @ingroup ParticleIntegrators
       */
      template <int dim>
      void
      register_particle_integrator (const std::string &name,
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
       * @ingroup ParticleIntegrators
       */
      template <int dim>
      Interface<dim> *
      create_particle_integrator (ParameterHandler &prm);


      /**
       * Declare the runtime parameters of the registered particle integrators.
       *
       * @ingroup ParticleIntegrators
       */
      template <int dim>
      void
      declare_parameters (ParameterHandler &prm);

      /**
       * Given a class name, a name, and a description for the parameter file
       * for a particle integrator, register it with the functions that
       * can declare their parameters and create these objects.
       *
       * @ingroup ParticleIntegrators
       */
#define ASPECT_REGISTER_PARTICLE_INTEGRATOR(classname, name, description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_PARTICLE_INTEGRATOR_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::Particle::Integrator::Interface<2>,classname<2> > \
    dummy_ ## classname ## _2d (&aspect::Particle::Integrator::register_particle_integrator<2>, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::Particle::Integrator::Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&aspect::Particle::Integrator::register_particle_integrator<3>, \
                                name, description); \
  }
    }
  }
}


#endif
