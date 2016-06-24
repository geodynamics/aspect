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

#ifndef __aspect__particle_interpolator_interface_h
#define __aspect__particle_interpolator_interface_h

#include <aspect/particle/particle.h>
#include <aspect/plugins.h>
#include <aspect/global.h>

#include <deal.II/base/point.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/distributed/tria.h>

namespace aspect
{
  namespace Particle
  {
    namespace Interpolator
    {
      using namespace dealii;

      /**
       * An abstract class defining virtual methods for performing
       * interpolation of particle properties to arbitrary points.
       *
       * @ingroup ParticleInterpolators
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
           * Perform an interpolation of the properties of the particles in
           * this cell onto a vector of positions in this cell.
           * Implementations of this function must return a vector of a vector
           * of doubles which contains a somehow computed
           * value of all particle properties at all given positions.
           *
           * @param [in] particles Reference to the particle map.
           * @param [in] positions The vector of positions where the properties
           * should be evaluated.
           * @param [in] cell An optional iterator to the cell containing the
           * particles. Not all callers will know the cell of the particles,
           * but providing the cell when known speeds up the interpolation
           * significantly.
           * @return A vector with as many entries as @p positions. Every entry
           * is a vector of interpolated tracer properties at this position.
           */
          virtual
          std::vector<std::vector<double> >
          properties_at_points(const std::multimap<types::LevelInd, Particle<dim> > &particles,
                               const std::vector<Point<dim> > &positions,
                               const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell = typename parallel::distributed::Triangulation<dim>::cell_iterator()) const = 0;

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
       * Return a list of names (separated by '|') of possible interpolator
       * classes for particles.
       */
      std::string
      interpolator_object_names ();


      /**
       * Register a particle interpolator so that it can be selected from
       * the parameter file.
       *
       * @param name A string that identifies the particle interpolator
       * @param description A text description of what this interpolator does and that
       * will be listed in the documentation of the parameter file.
       * @param declare_parameters_function A pointer to a function that can be
       * used to declare the parameters that this particle interpolator wants to read
       * from input files.
       * @param factory_function A pointer to a function that can create an
       * object of this particle interpolator.
       *
       * @ingroup ParticleInterpolators
       */
      template <int dim>
      void
      register_particle_interpolator (const std::string &name,
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
       * @ingroup ParticleInterpolators
       */
      template <int dim>
      Interface<dim> *
      create_particle_interpolator (ParameterHandler &prm);


      /**
       * Declare the runtime parameters of the registered particle interpolators.
       *
       * @ingroup ParticleInterpolators
       */
      template <int dim>
      void
      declare_parameters (ParameterHandler &prm);

      /**
       * Given a class name, a name, and a description for the parameter file
       * for a particle interpolator, register it with the functions that
       * can declare their parameters and create these objects.
       *
       * @ingroup ParticleInterpolators
       */
#define ASPECT_REGISTER_PARTICLE_INTERPOLATOR(classname, name, description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_PARTICLE_INTERPOLATOR_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::Particle::Interpolator::Interface<2>,classname<2> > \
    dummy_ ## classname ## _2d (&aspect::Particle::Interpolator::register_particle_interpolator<2>, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::Particle::Interpolator::Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&aspect::Particle::Interpolator::register_particle_interpolator<3>, \
                                name, description); \
  }
    }
  }
}


#endif
