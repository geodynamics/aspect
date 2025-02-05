/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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


#ifndef _aspect_boundary_velocity_interface_h
#define _aspect_boundary_velocity_interface_h

#include <aspect/plugins.h>
#include <aspect/simulator_access.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/utilities.h>

#include <deal.II/base/point.h>
#include <deal.II/base/parameter_handler.h>

#include <boost/core/demangle.hpp>
#include <typeinfo>


namespace aspect
{
  /**
   * A namespace in which we define everything that has to do with defining
   * the velocity boundary conditions.
   *
   * @ingroup BoundaryVelocities
   */
  namespace BoundaryVelocity
  {
    /**
     * A base class for parameterizations of velocity boundary conditions.
     *
     * @ingroup BoundaryVelocities
     */
    template <int dim>
    class Interface : public Plugins::InterfaceBase
    {
      public:
        /**
         * Return the velocity that is to hold at a particular position on
         * the boundary of the domain.
         *
         * @param boundary_indicator The boundary indicator of the part of the
         * boundary of the domain on which the point is located at which we
         * are requesting the velocity.
         * @param position The position of the point at which we ask for the
         * velocity.
         *
         * @return Boundary velocity at position @p position.
         */
        virtual
        Tensor<1,dim>
        boundary_velocity (const types::boundary_id boundary_indicator,
                           const Point<dim> &position) const = 0;
    };

    /**
     * A class that manages all boundary velocity objects.
     *
     * @ingroup BoundaryVelocities
     */
    template <int dim>
    class Manager : public Plugins::ManagerBase<Interface<dim>>, public SimulatorAccess<dim>
    {
      public:
        /**
         * A function that calls the boundary_velocity functions of all the
         * individual boundary velocity objects and uses the stored operators
         * to combine them.
         */
        Tensor<1,dim>
        boundary_velocity (const types::boundary_id boundary_indicator,
                           const Point<dim> &position) const;

        /**
         * A function that is used to register boundary velocity objects in such
         * a way that the Manager can deal with all of them without having to
         * know them by name. This allows the files in which individual
         * plugins are implemented to register these plugins, rather than also
         * having to modify the Manager class by adding the new boundary
         * velocity plugin class.
         *
         * @param name A string that identifies the boundary velocity model
         * @param description A text description of what this model does and that
         * will be listed in the documentation of the parameter file.
         * @param declare_parameters_function A pointer to a function that can be
         * used to declare the parameters that this boundary velocity model
         * wants to read from input files.
         * @param factory_function A pointer to a function that can create an
         * object of this boundary velocity model.
         */
        static
        void
        register_boundary_velocity (const std::string &name,
                                    const std::string &description,
                                    void (*declare_parameters_function) (ParameterHandler &),
                                    std::unique_ptr<Interface<dim>> (*factory_function) ());


        /**
         * Return the names of all prescribed boundary velocity models currently
         * used in the computation as specified in the input file. The function
         * returns a map between a boundary identifier and a pair. The
         * first part of the pair is a string that represents the prescribed
         * velocity components on this boundary (e.g. y, xz, or xyz) and the
         * second part is a vector of strings that represent the names of
         * boundary velocity plugins for this boundary.
         * If there are no prescribed boundary velocity plugins
         * for a particular boundary, this boundary identifier will not appear
         * in the map.
         *
         * @deprecated This function will be removed. Use the function
         * get_active_plugin_names() of the base class ManagerBase instead.
         */
        DEAL_II_DEPRECATED
        const std::map<types::boundary_id, std::pair<std::string,std::vector<std::string>>> &
        get_active_boundary_velocity_names () const;

        /**
         * Return pointers to all boundary velocity models
         * currently used in the computation, as specified in the input file.
         * The function returns a map between a boundary identifier and a vector
         * of unique pointers that represent the names of prescribed velocity
         * boundary models for this boundary. If there are no prescribed
         * boundary velocity plugins for a particular boundary this boundary
         * identifier will not appear in the map.
         *
         * @deprecated This function has been removed. Use the function
         * get_active_plugins() of the base class ManagerBase instead.
         */
        DEAL_II_DEPRECATED
        const std::map<types::boundary_id,std::vector<std::unique_ptr<BoundaryVelocity::Interface<dim>>>> &
        get_active_boundary_velocity_conditions () const;

        /**
         * Return a list of boundary indicators that indicate for
         * each active plugin which boundary id
         * it is responsible for. The list of active plugins can be
         * requested by calling get_active_plugins().
         */
        const std::vector<types::boundary_id> &
        get_active_plugin_boundary_indicators() const;

        /**
         * Return a component mask that indicates for the given
         * @p boundary_id which velocity components are prescribed by
         * this manager class. All plugins that are responsible
         * for this boundary use the same component mask.
         * The list of plugin objects can be
         * requested by calling get_active_plugins() and the
         * list of boundaries they are responsible for is
         * returned by get_active_plugin_boundary_indicators().
         */
        ComponentMask
        get_component_mask(const types::boundary_id boundary_id) const;

        /**
         * Return a set of boundary indicators for which boundary
         * velocities are prescribed.
         */
        const std::set<types::boundary_id> &
        get_prescribed_boundary_velocity_indicators () const;

        /**
         * Return a list of boundary ids on which the velocity is prescribed
         * to be zero (no-slip).
         */
        const std::set<types::boundary_id> &
        get_zero_boundary_velocity_indicators () const;

        /**
         * Return a list of boundary ids on which the velocity is prescribed
         * to be tangential (free-slip).
         */
        const std::set<types::boundary_id> &
        get_tangential_boundary_velocity_indicators () const;

        /**
         * Declare the parameters of all known boundary velocity plugins, as
         * well as the ones this class has itself.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         * This determines which boundary velocity objects will be created;
         * then let these objects read their parameters as well.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

        /**
         * Go through the list of all boundary velocity models that have been selected
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
        template <typename BoundaryVelocityType,
                  typename = typename std::enable_if_t<std::is_base_of<Interface<dim>,BoundaryVelocityType>::value>>
        DEAL_II_DEPRECATED
        bool
        has_matching_boundary_velocity_model () const;

        /**
         * Go through the list of all boundary velocity models that have been selected
         * in the input file (and are consequently currently active) and see
         * if one of them has the type specified by the template
         * argument or can be cast to that type. If so, return a reference
         * to it. If no boundary velocity model is active that matches the given type,
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
        template <typename BoundaryVelocityType,
                  typename = typename std::enable_if_t<std::is_base_of<Interface<dim>,BoundaryVelocityType>::value>>
        DEAL_II_DEPRECATED
        const BoundaryVelocityType &
        get_matching_boundary_velocity_model () const;

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


        /**
         * Exception.
         */
        DeclException1 (ExcBoundaryVelocityNameNotFound,
                        std::string,
                        << "Could not find entry <"
                        << arg1
                        << "> among the names of registered boundary velocity objects.");
      private:
        /**
         * A list of boundary indicators that indicate for
         * each plugin in the list of plugin_objects which boundary id
         * it is responsible for. By default each plugin
         * is active for all boundaries, but this list
         * can be modified by derived classes to limit the application
         * of plugins to specific boundaries.
         */
        std::vector<types::boundary_id> boundary_indicators;

        /**
         * A list of boundary indicators that indicate for
         * each plugin in the list of plugin_objects which components
         * it is responsible for. By default each plugin
         * is active for all components, but this list
         * can be modified by derived classes to limit the application
         * of plugins to specific boundaries.
         */
        std::vector<ComponentMask> component_masks;

        /**
         * A list of boundary velocity objects that have been requested in the
         * parameter file.
         *
         * @deprecated This variable is no longer used, but needed to issue a proper
         * error message in the function get_active_boundary_velocity_conditions().
         */
        std::map<types::boundary_id,std::vector<std::unique_ptr<BoundaryVelocity::Interface<dim>>>> boundary_velocity_objects;

        /**
         * Map from boundary id to a pair
         * ("components", list of "velocity boundary type"),
         * where components is of the format "[x][y][z]" and the velocity type is
         * mapped to one of the plugins of velocity boundary conditions (e.g.
         * "function"). If the components string is empty, it is assumed the
         * plugins are used for all components.
         *
         * @deprecated Remove this variable when the deprecated functions
         * get_active_boundary_velocity_names and
         * get_active_boundary_velocity_conditions are removed. Use the base class
         * variable plugin_names instead.
         */
        std::map<types::boundary_id, std::pair<std::string,std::vector<std::string>>> boundary_velocity_indicators;

        /**
         * A set of boundary indicators, on which velocities are prescribed.
         */
        std::set<types::boundary_id> prescribed_velocity_boundary_indicators;

        /**
         * A set of boundary indicators, on which velocities are prescribed to
         * zero (no-slip).
         */
        std::set<types::boundary_id> zero_velocity_boundary_indicators;

        /**
         * A set of boundary indicators, on which velocities are prescribed to
         * be tangential (free-slip).
         */
        std::set<types::boundary_id> tangential_velocity_boundary_indicators;
    };



    template <int dim>
    template <typename BoundaryVelocityType, typename>
    inline
    bool
    Manager<dim>::has_matching_boundary_velocity_model () const
    {
      for (const auto &boundary : boundary_velocity_objects)
        for (const auto &p : boundary.second)
          if (Plugins::plugin_type_matches<BoundaryVelocityType>(*p))
            return true;
      return false;
    }


    template <int dim>
    template <typename BoundaryVelocityType, typename>
    inline
    const BoundaryVelocityType &
    Manager<dim>::get_matching_boundary_velocity_model () const
    {
      AssertThrow(has_matching_boundary_velocity_model<BoundaryVelocityType> (),
                  ExcMessage("You asked BoundaryVelocity::Manager::get_boundary_velocity_model() for a "
                             "boundary velocity model of type <" + boost::core::demangle(typeid(BoundaryVelocityType).name()) + "> "
                             "that could not be found in the current model. Activate this "
                             "boundary velocity model in the input file."));

      for (const auto &boundary : boundary_velocity_objects)
        for (const auto &p : boundary)
          if (Plugins::plugin_type_matches<BoundaryVelocityType>(*p))
            return Plugins::get_plugin_as_type<BoundaryVelocityType>(*p);

      // We will never get here, because we had the Assert above. Just to avoid warnings.
      return Plugins::get_plugin_as_type<BoundaryVelocityType>(**(boundary_velocity_objects.begin()));
    }


    /**
     * Return a string that consists of the names of boundary velocity models that can
     * be selected. These names are separated by a vertical line '|' so
     * that the string can be an input to the deal.II classes
     * Patterns::Selection or Patterns::MultipleSelection.
     */
    template <int dim>
    std::string
    get_valid_model_names_pattern ();


    /**
     * Given a class name, a name, and a description for the parameter file
     * for a boundary velocity model, register it with the functions that
     * can declare their parameters and create these objects.
     *
     * @ingroup BoundaryVelocities
     */
#define ASPECT_REGISTER_BOUNDARY_VELOCITY_MODEL(classname, name, description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_BOUNDARY_VELOCITY_MODEL_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::BoundaryVelocity::Interface<2>,classname<2>> \
    dummy_ ## classname ## _2d (&aspect::BoundaryVelocity::Manager<2>::register_boundary_velocity, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::BoundaryVelocity::Interface<3>,classname<3>> \
    dummy_ ## classname ## _3d (&aspect::BoundaryVelocity::Manager<3>::register_boundary_velocity, \
                                name, description); \
  }
  }
}


#endif
