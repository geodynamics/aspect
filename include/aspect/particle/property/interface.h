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

#ifndef _aspect_particle_property_interface_h
#define _aspect_particle_property_interface_h

#include <aspect/particle/interface.h>
#include <aspect/global.h>

#include <aspect/particle/interpolator/interface.h>
#include <aspect/simulator_access.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/property_pool.h>
#include <deal.II/fe/fe_update_flags.h>

#include <memory>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      using namespace dealii::Particles;

      /**
       * A data structure with all inputs for the
       * Particle::Property::update_particle_properties() method.
       */
      template <int dim>
      struct ParticleUpdateInputs
      {
        public:
          /**
           * The solution vector at each particle position. This vector is
           * only filled if the update function requires the solution values.
           */
          std::vector<small_vector<double,50>> solution;

          /**
           * The solution gradients at each particle position.
           * This vector is only filled if the update function requires the
           * gradients of the solution values.
           */
          std::vector<small_vector<Tensor<1,dim>,50>> gradients;

          /**
           * Cell iterator of the cell that is currently being updated.
           * This allows for evaluating additional properties at the cell vertices,
           * or to query the cell for material ids, neighbors, or other
           * information that is not available solely from the particles.
           */
          typename DoFHandler<dim>::active_cell_iterator current_cell;
      };

      /**
       * This class is used to store all the necessary information to translate
       * between the data structure of the particle properties (a flat vector of
       * doubles) and the semantic meaning of these properties. It contains
       * information about the three layers of particle property information:
       *
       * By 'property plugins' we mean each separate class that is derived from
       * aspect::Particle::Property::Interface<dim>, and that is selected in
       * the input file. This means in any model there are as many property
       * plugins as entries in the 'List of particle properties' input parameter.
       *
       * Each plugin can create one or more 'property fields'. Property fields
       * are interpreted as distinctly named particle properties. Most plugins
       * contain only one field, but some group distinctly named properties
       * into groups, such as the 'InitialComposition' plugin, which creates
       * one property field per compositional field and names the property
       * fields according to the compositional field names. When writing particle
       * data to output files each 'property field' will be written into a
       * separate output field.
       *
       * Last each field can contain several 'property components'
       * if it represents a vector or tensor property. These components can
       * not be named individually, but it is still important to be able to
       * know how many components belong to a particular field.
       *
       * Information that is often required by other algorithms is for example
       * the number of components (= number of doubles in the property vector)
       * of a particle field or plugin, and its position within the particle
       * property vector. All of this information might be required within
       * loops over all fields, or for a specific field either identified by
       * its index, or by its name.
       *
       * @ingroup ParticleProperties
       */
      class ParticlePropertyInformation
      {
        public:
          /**
           * Empty default constructor.
           */
          ParticlePropertyInformation();

          /**
           * Constructor. Initialize the various arrays of this structure with the
           * given property information collected from the individual plugins.
           *
           * @p property_information A vector that contains one vector per
           * property plugin. Each of these vectors contains one or more pairs
           * that represent a property field name and the number of components
           * for this field. The input argument can be constructed by
           * concatenating the output of the
           * Particle::Property::Interface<dim>::get_property_information()
           * functions of all property plugins.
           */
          ParticlePropertyInformation(const std::vector<std::vector<std::pair<std::string,unsigned int>>> &property_information);

          /**
           * Checks if the particle property specified by @p name exists
           * in this model.
           */
          bool
          fieldname_exists(const std::string &name) const;

          /**
           * Get the field index of the particle property specified by @p name.
           */
          unsigned int
          get_field_index_by_name(const std::string &name) const;

          /**
           * Get the field index of the particle property specified by @p name.
           */
          std::string
          get_field_name_by_index(const unsigned int field_index) const;

          /**
           * Get the data position of the first component of the particle
           * property specified by @p name in the property vector of a particle.
           */
          unsigned int
          get_position_by_field_name(const std::string &name) const;

          /**
           * Get the number of components of the particle property specified
           * by @p name.
           */
          unsigned int
          get_components_by_field_name(const std::string &name) const;

          /**
           * Get the data position of the first component of the particle
           * property specified by @p field_index in the property vector
           * of a particle.
           */
          unsigned int
          get_position_by_field_index(const unsigned int field_index) const;

          /**
           * Get the number of components of the particle property specified
           * by @p field_index.
           */
          unsigned int
          get_components_by_field_index(const unsigned int field_index) const;

          /**
           * Get the data position of the first component of the particle
           * property specified by @p plugin_index in the property vector
           * of a particle.
           */
          unsigned int
          get_position_by_plugin_index(const unsigned int plugin_index) const;

          /**
           * Get the number of components of the particle property specified
           * by @p plugin_index.
           */
          unsigned int
          get_components_by_plugin_index(const unsigned int plugin_index) const;

          /**
           * Get the number of fields of the particle property specified
           * by @p plugin_index.
           */
          unsigned int
          get_fields_by_plugin_index(const unsigned int plugin_index) const;


          /**
           * Return the number of active particle property plugins.
           */
          unsigned int
          n_plugins() const;

          /**
           * Return the number of active particle property fields.
           */
          unsigned int
          n_fields() const;

          /**
           * Return the number of active particle property components.
           */
          unsigned int
          n_components() const;

        private:
          /**
           * A vector of all property field names.
           */
          std::vector<std::string> field_names;

          /**
           * A vector containing the number of components per property field.
           */
          std::vector<unsigned int> components_per_field;

          /**
           * A vector containing the position index of the first data component
           * of each field in the property vector of every particle.
           */
          std::vector<unsigned int> position_per_field;

          /**
           * A vector containing the number of property fields per property
           * plugin.
           */
          std::vector<unsigned int> fields_per_plugin;

          /**
           * A vector containing the number of components per property plugin.
           */
          std::vector<unsigned int> components_per_plugin;

          /**
           * A vector containing the position index of the first data component
           * of each plugin in the property vector of every particle.
           */
          std::vector<unsigned int> position_per_plugin;

          /**
           * The number of doubles needed to represent a particle's
           * additional properties.
           */
          unsigned int number_of_components;

          /**
           * The number of distinctly named particle property fields.
           */
          unsigned int number_of_fields;

          /**
           * The number of active particle property plugins.
           */
          unsigned int number_of_plugins;
      };

      enum UpdateTimeFlags
      {
        /**
         * Never update the initially set properties. This is the default
         * behavior, which is sufficient for particle properties that are
         * set at the beginning of the model and constant for the whole
         * simulation time.
         */
        update_never,
        /**
         * Update the particle properties before every output. This is
         * sufficient for all passive particle properties that depend on the
         * current solution, like the current velocity or pressure.
         */
        update_output_step,
        /**
         * Update the particle properties every nonlinear iteration. This is only necessary
         * if the properties at the output time depend on some sort of time
         * integration of solution properties or time varying particle
         * properties are used while solving the model problem.
         */
        update_time_step
      };

      /**
       * This enum controls how to initialize the properties of particles that
       * have been added later than the initial particle creation, e.g. to
       * improve the load balance or to prevent empty cells.
       */
      enum InitializationModeForLateParticles
      {
        /**
         * Initialize the particle as if it were created at the beginning of
         * the model at its current position with the current solution.
         */
        initialize,
        /**
         * Use the interpolated properties of the surrounding particles as
         * calculated by the selected interpolator.
         */
        interpolate,
        /**
         * Use the interpolated properties of the surrounding particles as
         * calculated by the selected interpolator except for particles in
         * boundary cells. These will use the boundary condition of the
         * compositional fields instead. This mode only makes sense for
         * properties that are associated with compositional fields through
         * the parameter 'Compositional fields/Mapped particle properties'.
         */
        interpolate_respect_boundary,
        /**
         * Initialize the particle properties to zero. If the property is
         * updated over time its update function is called as usual, if not
         * the property will remain zero throughout the model run.
         */
        initialize_to_zero
      };

      /**
       * Interface provides an example of how to extend the Particle class to
       * include related particle data. This allows users to attach
       * scalars/vectors/tensors/etc to particles and ensure they are
       * transmitted correctly over MPI and written to output files.
       *
       * @ingroup ParticleProperties
       */
      template <int dim>
      class Interface : public ParticleInterfaceBase
      {
        public:
          /**
           * Initialization function. This function is called once at the
           * creation of every particle for every property to initialize its
           * value.
           *
           * @param [in] position The current particle position.
           * @param [in,out] particle_properties The properties of the particle
           * that is initialized within the call of this function. The purpose
           * of this function should be to extend this vector by a number of
           * properties.
           */
          virtual
          void
          initialize_one_particle_property (const Point<dim> &position,
                                            std::vector<double> &particle_properties) const;

          /**
           * Update function. This function is called every time an update is
           * requested by need_update() for every cell for every property.
           * It is expected to update the properties of all particles in the
           * given range @p particles, which are all in one cell.
           * It is obvious that
           * this function is called a lot, so its code should be efficient.
           * The interface provides a default implementation that does nothing,
           * therefore derived plugins that do not require an update do not
           * need to implement this function.
           *
           * @param [in] inputs A struct of type ParticleUpdateInputs that contains
           * all necessary inputs to compute the particle updates. See
           * the documentation of this struct in
           * include/aspect/particle/property/interface.h for a list of all
           * available inputs.
           *
           * @param [in,out] particles The particles that are to be updated
           * within this function.
           */
          virtual
          void
          update_particle_properties (const ParticleUpdateInputs<dim> &inputs,
                                      typename ParticleHandler<dim>::particle_iterator_range &particles) const;

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
           * @param [in] solution The values of the solution variables at the
           * current particle position.
           *
           * @param [in] gradients The gradients of the solution variables at
           * the current particle position.
           *
           * @param [in,out] particle The particle that is updated within
           * the call of this function. The particle location can be accessed
           * using particle->get_location() and its properties using
           * particle->get_properties().
           *
           * @deprecated This version of the function is deprecated.
           * Use update_particle_properties() instead, which allows to
           * update all particles of a cell in one function call.
           */
          DEAL_II_DEPRECATED
          virtual
          void
          update_particle_property (const unsigned int data_position,
                                    const Vector<double> &solution,
                                    const std::vector<Tensor<1,dim>> &gradients,
                                    typename ParticleHandler<dim>::particle_iterator &particle) const;


          /**
           * Returns an enum, which determines at what times particle properties
           * are updated. The default implementation returns update_never, which
           * signals that particle properties should never be updated.
           * See the documentation of UpdateTimeFlags for a list of possible
           * values and examples for their use. Every
           * plugin that implements this function should return the value
           * appropriate for its purpose, unless it does not need any update,
           * which is the default. This option saves considerable computation
           * time in cases, when no plugin needs to update particle properties
           * over time.
           */
          virtual
          UpdateTimeFlags
          need_update () const;

          /**
           * Return which data of the solution component @p component
           * has to be provided to update the current particle property.
           *
           * Note that particle properties can only ask for update_default
           * (no data), update_values (solution values), and update_gradients
           * (solution gradients). All other update flags will have no effect.
           *
           * As an example consider a particle property that depends on the
           * solution values and gradients of the velocity field. In this case
           * the function should return update_values | update_gradients if the
           * @p component is one of the velocity components, and update_default
           * otherwise.
           *
           * @param component The component of the solution which is to be
           * evaluated.
           *
           * @return The necessary update flags for the solution component
           * @p component that is required for this particle property.
           */
          virtual
          UpdateFlags
          get_update_flags (const unsigned int component) const;

          /**
           * Return which data has to be provided to update all properties.
           * Note that particle properties can only ask for update_default
           * (no data), update_values (solution values), and update_gradients
           * (solution gradients). All other update flags will have no effect.
           *
           * @return The necessary update flags for this particle property.
           *
           * @deprecated This function is deprecated. Use the above version of
           * get_update_flags() instead.
           */
          DEAL_II_DEPRECATED
          virtual
          UpdateFlags
          get_needed_update_flags () const;

          /**
           * Returns an enum, which determines how this particle property is
           * initialized for particles that are created later than the initial
           * particle generation, e.g. to balance the particle load or prevent
           * empty cells. The default implementation returns
           * interpolate, which will use the particle interpolator to
           * set the new particle properties to a value that is interpolated
           * from the other particles in the cell.
           * See the documentation of InitializationModeForLateParticles for a
           * list of possible values and examples for their use. Every
           * plugin that implements this function should return the value
           * appropriate for its purpose, unless it wants to use the default
           * value.
           */
          virtual
          InitializationModeForLateParticles
          late_initialization_mode () const;

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
          std::vector<std::pair<std::string, unsigned int>>
          get_property_information() const = 0;

          /**
           * Set the position of this property in the particle property vector.
           */
          virtual
          void
          set_data_position (const unsigned int data_position);

          /**
           * Get the position of this property in the particle property vector.
           */
          virtual
          unsigned int
          get_data_position () const;

        protected:
          /**
           * Store the position of the particle property in the particle property vector.
           * If the property has multiple components, the first component is stored
           * and all other components are stored consecutively after the first one.
           */
          unsigned int data_position;
      };

      /**
       * A particle property that provides storage space for
       * the properties that particle integrators need to
       * store. This is an internal property that is not
       * intended for use outside of the particle integrators
       * and that will not be written to output files.
       *
       * @ingroup ParticleProperties
       */
      template <int dim>
      class IntegratorProperties : public Interface<dim>
      {
        public:
          /**
           * Initialization function. Since these properties are set and used
           * by the integrator this function only resizes them to the correct
           * size, but does not need to do any initialization.
           */
          void
          initialize_one_particle_property (const Point<dim> &position,
                                            std::vector<double> &particle_properties) const override;

          /**
           * Set up the information about the names and number of components
           * this property requires. This depends on the chosen integration scheme.
           *
           * @return A vector that contains pairs of the property names and the
           * number of components this property plugin defines.
           */
          std::vector<std::pair<std::string, unsigned int>>
          get_property_information() const override;

          /**
           * Read the parameters this class needs to determine which integrator is used,
           * and therefore how many properties to reserve.
           */
          void
          parse_parameters (ParameterHandler &prm) override;

        private:
          /**
           * The number of integrator properties to store. This variable is initialized in
           * parse_parameters().
           */
          unsigned int n_integrator_properties;
      };



      /**
       * Manager class of properties - This class sets the data of the
       * collection of particles and updates it over time if requested by the
       * user selected properties.
       */
      template <int dim>
      class Manager : public Plugins::ManagerBase<Interface<dim>>, public SimulatorAccess<dim>
      {
        public:
          /**
           * Initialization function. This function is called once at the
           * beginning of the program after parse_parameters is run.
           */
          void
          initialize () override;

          /**
           * Initialization function for particle properties. This function is
           * called once for each of the particles of a particle
           * collection after it was created.
           */
          void
          initialize_one_particle (typename ParticleHandler<dim>::particle_iterator &particle) const;

          /**
           * Initialization function for particle properties. This function is
           * called once for each of the particles of a particle
           * collection that were created later than the initial particle
           * generation.
           */
          std::vector<double>
          initialize_late_particle (const Point<dim> &particle_location,
                                    const ParticleHandler<dim> &particle_handler,
                                    const Interpolator::Interface<dim> &interpolator,
                                    const typename parallel::distributed::Triangulation<dim>::active_cell_iterator &cell = typename parallel::distributed::Triangulation<dim>::active_cell_iterator()) const;

          /**
           * Update function for particle properties. This function is
           * called once every time step for every cell.
           *
           * @param inputs A struct of type ParticleUpdateInputs that contains
           * all necessary inputs to compute the particle updates. See
           * the documentation of this struct in
           * include/aspect/particle/property/interface.h for a list of all
           * available inputs.
           * @param particles The particles that are to be updated within
           * this function.
           */
          void
          update_particles (ParticleUpdateInputs<dim> &inputs,
                            typename ParticleHandler<dim>::particle_iterator_range &particles) const;

          /**
           * Returns an enum, which denotes at what time this class needs to
           * update particle properties. The result of this class is a
           * combination of the need_update() functions of all individual
           * properties that are selected. More precise, it will choose to
           * update the particle properties as often as the plugin that needs the
           * most frequent update option requires. This saves considerable
           * computation time, e.g. in cases when no plugin needs to update particle
           * properties over time, because the solution does not need to be
           * evaluated in this case.
           */
          UpdateTimeFlags
          need_update () const;

          /**
           * Return which data has to be provided to update all properties.
           * Note that particle properties can only ask for update_default
           * (no data), update_values (solution values), and update_gradients
           * (solution gradients). All other update flags will have no effect.
           * The result of this function is a combination of the
           * get_update_flags() functions of all individual properties
           * that are selected.
           *
           * @return A vector that contains the update flags that are
           * required to update all particle properties. The vector has as many entries
           * as there solution components.
           */
          std::vector<UpdateFlags>
          get_update_flags () const;

          /**
           * Checks if the particle plugin specified by @p name exists
           * in this model.
           */
          bool
          plugin_name_exists(const std::string &name) const;

          /**
           * Checks if the particle property plugin specified by @p first
           * is executed before another particle property plugin specified
           * by @p second.
           *
           * Throws an assert when one of the plugin names does not
           * exist. You can use the function plugin_name_exists() to
           * check in advance whether a plugin exists
           */
          bool
          check_plugin_order(const std::string &first, const std::string &second) const;

          /**
           * Get the plugin index of the particle plugin specified by @p name.
           */
          unsigned int get_plugin_index_by_name(const std::string &name) const;

          /**
           * Go through the list of all particle properties that have been selected
           * in the input file (and are consequently currently active) and return
           * true if one of them has the desired type specified by the template
           * argument.
           *
           * This function can only be called if the given template type (the first template
           * argument) is a class derived from the Interface class in this namespace.
           *
           * @deprecated Instead of this function, use the
           *   Plugins::ManagerBase::has_matching_active_plugin() and
           *   Plugins::ManagerBase::get_matching_active_plugin() functions of the base
           *   class of the current class.
           */
          template <typename ParticlePropertyType,
                    typename = typename std::enable_if_t<std::is_base_of<Interface<dim>,ParticlePropertyType>::value>>
          DEAL_II_DEPRECATED
          bool
          has_matching_property () const;

          /**
           * Go through the list of all particle properties that have been selected
           * in the input file (and are consequently currently active) and see
           * if one of them has the type specified by the template
           * argument or can be cast to that type. If so, return a reference
           * to it. If no property is active that matches the given type,
           * throw an exception.
           *
           * This function can only be called if the given template type (the first template
           * argument) is a class derived from the Interface class in this namespace.
           *
           * @deprecated Instead of this function, use the
           *   Plugins::ManagerBase::has_matching_active_plugin() and
           *   Plugins::ManagerBase::get_matching_active_plugin() functions of the base
           *   class of the current class.
           */
          template <typename ParticlePropertyType,
                    typename = typename std::enable_if_t<std::is_base_of<Interface<dim>,ParticlePropertyType>::value>>
          DEAL_II_DEPRECATED
          const ParticlePropertyType &
          get_matching_property () const;

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
          const ParticlePropertyInformation &
          get_data_info() const;

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
           * does and that will be listed in the documentation of the
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
                                      std::unique_ptr<Property::Interface<dim>> (*factory_function) ());


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
          parse_parameters (ParameterHandler &prm) override;

          /**
           * @brief Set the particle manager index for all particle properties
           *
           * @param particle_manager_index The index of the particle manager.
           */
          void set_particle_manager_index(unsigned int particle_manager_index);

          /**
           * For the current plugin subsystem, write a connection graph of all of the
           * plugins we know about, in the format that the
           * programs dot and neato understand. This allows for a visualization of
           * how all of the plugins that ASPECT knows about are interconnected, and
           * connect to other parts of the ASPECT code.
           *
           * @param output_stream The stream to write the output to.
           */
          static
          void
          write_plugin_graph (std::ostream &output_stream);

        private:
          /**
           * Stores the index to the particle manager, to which this manager belongs.
           */
          unsigned int particle_manager_index;

          /**
           * A class that stores all information about the particle properties,
           * their association with property plugins and their storage pattern.
           */
          ParticlePropertyInformation property_information;
      };

      /* -------------------------- inline and template functions ---------------------- */


      template <int dim>
      template <typename ParticlePropertyType, typename>
      inline
      bool
      Manager<dim>::has_matching_property () const
      {
        return this->template has_matching_active_plugin<ParticlePropertyType>();
      }


      template <int dim>
      template <typename ParticlePropertyType, typename>
      inline
      const ParticlePropertyType &
      Manager<dim>::get_matching_property () const
      {
        return this->template get_matching_active_plugin<ParticlePropertyType>();
      }


      /**
       * Given a class name, a name, and a description for the parameter file for
       * a particle property, register it with the aspect::Particle:: class.
       *
       * @ingroup Particle
       */
#define ASPECT_REGISTER_PARTICLE_PROPERTY(classname,name,description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_PARTICLE_PROPERTY_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::Particle::Property::Interface<2>,classname<2>> \
    dummy_ ## classname ## _2d (&aspect::Particle::Property::Manager<2>::register_particle_property, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::Particle::Property::Interface<3>,classname<3>> \
    dummy_ ## classname ## _3d (&aspect::Particle::Property::Manager<3>::register_particle_property, \
                                name, description); \
  }

    }
  }
}

#endif
