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

#ifndef __aspect__particle_property_interface_h
#define __aspect__particle_property_interface_h

#include <aspect/particle/particle.h>
#include <aspect/simulator_access.h>
#include <aspect/plugins.h>

#include <deal.II/base/std_cxx1x/shared_ptr.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {

      enum UpdateTimeFlags
      {
        /**
         * Never update the initially set properties. This is the default
         * behaviour, which is sufficient for particle properties that are
         * set at the beginning of the model and constant for the whole
         * simulation time.
         */
        update_never,
        /**
         * Update the tracer properties before every output. This is
         * sufficient for all passive tracer properties that depend on the
         * current solution, like the current velocity or pressure.
         */
        update_output_step,
        /**
         * Update the tracer properties every timestep. This is only necessary
         * if the properties at the output time depend on some sort of time
         * integration of solution properties or time varying particle
         * properties are used while solving the model problem.
         */
        update_time_step
      };

      /**
        * Interface provides an example of how to extend the Particle
        * class to include related particle data. This allows users to attach
        * scalars/vectors/tensors/etc to particles and ensure they are
        * transmitted correctly over MPI and written to output files.
        *
        * @ingroup ParticleProperties
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
           * Initialization function. This function is called once at the
           * beginning of the program after parse_parameters is run.
           */
          virtual
          void
          initialize ();

          /**
           * Initialization function. This function is called once at the
           * creation of every particle for every property to initialize its
           * value.
           *
           * @param [in] position The current particle position.
           *
           * @param [in] solution The values of the solution variables at the
           * current particle position.
           *
           * @param [in] gradients The gradients of the solution variables at
           * the current particle position.
           *
           * @param [in,out] particle_properties The properties of the particle
           * that is initialized within the call of this function. The purpose
           * of this function should be to extend this vector by a number of
           * properties.
           */
          virtual
          void
          initialize_one_particle_property (const Point<dim> &position,
                                            const Vector<double> &solution,
                                            const std::vector<Tensor<1,dim> > &gradients,
                                            std::vector<double> &particle_properties) const;

          /**
           * Update function. This function is called every time an update is
           * request by need_update() for every particle for every property.
           * It is obvious that
           * this function is called a lot, so its code should be efficient.
           * The interface provides a default implementation that does nothing,
           * therefore derived plugins that do not require an update do not
           * need to implement this function.
           *
           * @param [in] data_position An unsigned integer that denotes which
           * component of the particle property vector is associated with the
           * current property. For properties that own several components it
           * denotes the first component of this property, all other components
           * fill consecutive entries in the @p particle_properties vector.
           *
           * @param [in] position The current particle position.
           *
           * @param [in] solution The values of the solution variables at the
           * current particle position.
           *
           * @param [in] gradients The gradients of the solution variables at
           * the current particle position.
           *
           * @param [in,out] particle_properties The properties of the particle
           * that is updated within the call of this function.
           */
          virtual
          void
          update_one_particle_property (const unsigned int data_position,
                                        const Point<dim> &position,
                                        const Vector<double> &solution,
                                        const std::vector<Tensor<1,dim> > &gradients,
                                        std::vector<double> &particle_properties) const;

          /**
           * Returns an enum, which determines at what times particle properties
           * are updated. The default implementation returns update_never, which
           * signals that particle properties should never be updated.
           * See the documentation of UpdateTimeFlags for a list of possible
           * values and examples for their use. Every
           * plugin that implements this function should return the value
           * appropriate for its purpose, unless it does not need any update,
           * which is the default. This option saves considerable computation
           * time in cases, when no plugin needs to update tracer properties
           * over time.
           */
          virtual
          UpdateTimeFlags
          need_update () const;

          /**
           * Set up the information about the names and number of components
           * this property requires. Derived classes need to implement this
           * function.
           *
           * The purpose of this function is to return a vector of pairs of a
           * property name and the number of components associated with this
           * name (e.g. 1 for scalar properties,
           * n for n-dimensional vectors).
           *
           * @return A vector that contains pairs of the property names and the
           * number of components this property plugin defines.
           */
          virtual
          std::vector<std::pair<std::string, unsigned int> >
          get_property_information() const = 0;


          /**
           * Declare the parameters this class takes through input files.
           * Derived classes should overload this function if they actually do
           * take parameters; this class declares a fall-back function that
           * does nothing, so that property classes that do not take any
           * parameters do not have to do anything at all.
           *
           * This function is static (and needs to be static in derived
           * classes) so that it can be called without creating actual objects
           * (because declaring parameters happens before we read the input
           * file and thus at a time when we don't even know yet which
           * property objects we need).
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Read the parameters this class declares from the parameter file.
           * The default implementation in this class does nothing, so that
           * derived classes that do not need any parameters do not need to
           * implement it.
           */
          virtual
          void
          parse_parameters (ParameterHandler &prm);
      };



      /**
       * Manager class of properties - This class sets the data of the
       * collection of particles and updates it over time if requested by the
       * user selected properties.
       */
      template <int dim>
      class Manager : public SimulatorAccess<dim>
      {
        public:
          /**
           * Empty constructor for Manager
           */
          Manager ();

          /**
           * Destructor for Manager
           */
          ~Manager ();

          /**
           * Initialization function. This function is called once at the
           * beginning of the program after parse_parameters is run.
           */
          void
          initialize ();

          /**
           * Initialization function for particle properties. This function is
           * called once for each of the particles of a particle
           * collection after it was created.
           */
          void
          initialize_one_particle (Particle<dim> &particle,
                                   const Vector<double> &solution,
                                   const std::vector<Tensor<1,dim> > &gradients) const;

          /**
           * Update function for particle properties. This function is
           * called once every time step for every particle.
           */
          void
          update_one_particle (Particle<dim> &particle,
                               const Vector<double> &solution,
                               const std::vector<Tensor<1,dim> > &gradients) const;

          /**
           * Returns an enum, which denotes at what time this class needs to
           * update tracer properties. The result of this class is a
           * combination of the need_update() functions of all individual
           * properties that are selected. More precise, it will choose to
           * update the tracer properties as often as the plugin that needs the
           * most frequent update option requires. This saves considerable
           * computation time, e.g. in cases when no plugin needs to update tracer
           * properties over time, because the solution does not need to be
           * evaluated in this case.
           */
          UpdateTimeFlags
          need_update () const;

          /**
           * Get the number of components required to represent this particle's
           * properties.
           *
           * @return Number of doubles required to represent this particle's
           * additional properties.
           */
          unsigned int
          get_n_property_components () const;

          /**
           * Get the size in number of bytes required to represent this
           * particle's properties for communication. This is essentially the
           * space needed for the property components plus the space for the
           * particle position and the space needed for its ID.
           *
           * @return Number of bytes required to represent this particle.
           */
          std::size_t
          get_particle_size () const;

          /**
           * Get the names and number of components of particle properties.
           *
           * @return A vector of pairs for each property name and the
           * corresponding number of components attached to particles.
           */
          const std::vector<std::pair<std::string,unsigned int> > &
          get_data_info() const;

          /**
           * Get the position of the property specified by name in the property
           * vector of the particles.
           */
          unsigned int
          get_property_component_by_name(const std::string &name) const;

          /**
           * A function that is used to register particle property
           * objects in such a way that the Manager can deal with all of them
           * without having to know them by name. This allows the files in which
           * individual properties are implemented to register these
           * properties, rather than also having to modify the Manager class
           * by adding the new properties class.
           *
           * @param name The name under which this particle property
           * is to be called in parameter files.
           * @param description A text description of what this particle property
           *  does and that will be listed in the documentation of the
           * parameter file.
           * @param declare_parameters_function A pointer to a function that
           * declares the parameters for this property.
           * @param factory_function A pointer to a function that creates such a
           * property object and returns a pointer to it.
           */
          static
          void
          register_particle_property (const std::string &name,
                                      const std::string &description,
                                      void (*declare_parameters_function) (ParameterHandler &),
                                      Property::Interface<dim> *(*factory_function) ());


          /**
           * Declare the parameters this class takes through input files.
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Read the parameters this class declares from the parameter file.
           */
          void
          parse_parameters (ParameterHandler &prm);

        private:
          /**
           * A list of property objects that have been requested in the
           * parameter file.
           */
          std::list<std_cxx1x::shared_ptr<Interface<dim> > > property_list;

          /**
           * A map between names of properties and the position of their
           * first data component in the particle property vector.
           */
          std::map<std::string,unsigned int> property_position_map;

          /**
           * The number of doubles needed to represent a tracer's
           * additional properties.
           */
          unsigned int n_property_components;

          /**
           * Vector of the names and number of components of the properties
           * that are selected in this
           * model. This vector has as many entries as individually named
           * fields in the properties, which does not need to be the size of
           * the property_list, because a single property plugin can define
           * several properties (scalar or vector).
           */
          std::vector<std::pair<std::string,unsigned int> > property_component_list;

          /**
           * Vector of the data positions of individual property plugins.
           * This vector has as many components as property plugins selected.
           * It can be different from property_position_map, because
           * single plugins can define multiple data fields.
           */
          std::vector<unsigned int> positions;
      };


      /**
       * Given a class name, a name, and a description for the parameter file for
       * a tracer property, register it with the aspect::Particle:: class.
       *
       * @ingroup Particle
       */
#define ASPECT_REGISTER_PARTICLE_PROPERTY(classname,name,description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_PARTICLE_PROPERTY_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::Particle::Property::Interface<2>,classname<2> > \
    dummy_ ## classname ## _2d (&aspect::Particle::Property::Manager<2>::register_particle_property, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::Particle::Property::Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&aspect::Particle::Property::Manager<3>::register_particle_property, \
                                name, description); \
  }

    }
  }
}

#endif

