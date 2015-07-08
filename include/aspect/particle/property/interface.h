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

#ifndef __aspect__particle_property_interface_h
#define __aspect__particle_property_interface_h

#include <aspect/particle/particle.h>
#include <aspect/particle/definitions.h>
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
         * Never update the initially set properties.
         */
        update_never,
        /**
         * Update the tracer properties before outputting them. This is
         * sufficient for all tracer properties that depend on the current
         * solution, like current velocity or pressure.
         */
        update_output_step,
        /**
         * Update the tracer properties every timestep. This is only necessary
         * if the properties at the output time depend on some sort of time
         * integration of solution properties.
         */
        update_time_step
      };

      /**
        * Interface provides an example of how to extend the Particle
        * class to include related particle data. This allows users to attach
        * scalars/vectors/tensors/etc to particles and ensure they are
        * transmitted correctly over MPI and written to output files.
        */
      template <int dim>
      class Interface
      {
        public:

          /**
           * Initialization function. This function is called once at the
           * beginning of the program after parse_parameters is run.
           */
          virtual
          void
          initialize ();

          /**
           * Initialization function. This function is called once at the
           * creation of every particle to set up it's properties.
           */
          virtual
          void
          initialize_particle (std::vector<double> &data,
                               const Point<dim> &position,
                               const Vector<double> &solution,
                               const std::vector<Tensor<1,dim> > &gradients);

          /**
           * Update function. This function is called every timestep for
           * every particle to update it's properties. It is obvious that
           * this function is called a lot, so its code should be efficient.
           */
          virtual
          void
          update_particle (unsigned int &data_position,
                           std::vector<double> &particle_properties,
                           const Point<dim> &position,
                           const Vector<double> &solution,
                           const std::vector<Tensor<1,dim> > &gradients);

          /**
           * Returns a bool, which is false in the default implementation,
           * telling the property manager that no update is needed. Every
           * plugin that implements this function should return true. This
           * saves considerable computation time in cases, when no plugin needs
           * to update tracer properties over time.
           */
          virtual
          UpdateTimeFlags
          need_update ();

          virtual
          void
          data_length(std::vector<unsigned int> &length) const = 0;

          /**
           * Set up the name information for the particle property
           *
           * @param [in,out] names Vector that contains the property name
           */
          virtual
          void
          data_names(std::vector<std::string> &names) const = 0;


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
       * Manager class of properties - This class sets the data of the particles
       * and updates it over time if requested by the user selected properties
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
          virtual
          ~Manager ();

          /**
           * Initialization function. This function is called once at the
           * beginning of the program after parse_parameters is run.
           */
          virtual
          void
          initialize ();

          /**
           * Initialization function for particle properties. This function is
           * called once at the creation of a particle
           */
          virtual
          void
          initialize_particle (Particle<dim> &particle,
                               const Vector<double> &solution,
                               const std::vector<Tensor<1,dim> > &gradients);

          /**
           * Update function for particle properties. This function is
           * called once every timestep for every particle
           */
          virtual
          void
          update_particle (Particle<dim> &particle,
                           const Vector<double> &solution,
                           const std::vector<Tensor<1,dim> > &gradients);

          /**
           * Returns a bool, which is false if no selected plugin needs to
           * update tracer properties over time. This saves considerable
           * computation time in cases, when no plugin needs to update tracer
           * properties over time, because the solution does not need to be
           * evaluated at tracer positions in this case.
           */
          virtual
          UpdateTimeFlags
          need_update ();

          /**
           * Get the number of doubles required to represent this particle's
           * properties for communication.
           *
           * @return Number of doubles required to represent this particle
           */
          unsigned int
          get_data_len () const;

          /**
           * Add the MPI data description for this particle type to the vector.
           *
           * @param return This function returns a vector of names of the
           * particle properties.
           */
          void
          get_data_info (std::vector<std::string> &names,
                         std::vector<unsigned int> &length) const;

          /**
           * Get the position of the property specified by name in the property
           * vector of the particles.
           */
          unsigned int
          get_property_component_by_name(const std::string &name) const;

          /**
           * A function that is used to register visualization postprocessor
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
          virtual
          void
          parse_parameters (ParameterHandler &prm);

        private:
          /**
           * A list of property objects that have been requested in the
           * parameter file.
           */
          std::list<std_cxx1x::shared_ptr<Interface<dim> > > property_list;

          /**
           * A map between names of properties and their data component.
           */
          std::map<std::string,unsigned int> property_component_map;

          /**
           * The number of doubles needed to represent a typical tracer
           */
          unsigned int data_len;

          /**
           * Vector of the names of the properties that are selected in this
           * model. This vector has as many components as individually named
           * fields in the properties.
           */
          std::vector<std::string> names;

          /**
           * Vector of the data length of the properties that are selected in
           * this model. This vector has as many components as individually named
           * fields in the properties.
           */
          std::vector<unsigned int> length;

          /**
           * Vector of the data positions of individual property plugins.
           * This vector has as many components as property plugins selected.
           * It can be different from length of names or length, because
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

