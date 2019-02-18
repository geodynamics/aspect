/*
  Copyright (C) 2011 - 2018 by the authors of the ASPECT code.

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
    using namespace dealii;

    /**
     * A base class for parameterizations of velocity boundary conditions.
     *
     * @ingroup BoundaryVelocities
     */
    template <int dim>
    class Interface
    {
      public:
        /**
         * Destructor. Made virtual to enforce that derived classes also have
         * virtual destructors.
         */
        virtual ~Interface();

        /**
         * Initialization function. This function is called once at the
         * beginning of the program after parse_parameters is run and after
         * the SimulatorAccess (if applicable) is initialized.
         */
        virtual
        void
        initialize ();

        /**
         * A function that is called at the beginning of each time step. The
         * default implementation of the function does nothing, but derived
         * classes that need more elaborate setups for a given time step may
         * overload the function.
         *
         * The point of this function is to allow complex boundary velocity
         * models to do an initialization step once at the beginning of each
         * time step. An example would be a model that needs to call an
         * external program to compute positions for a set of plates.
         */
        virtual
        void
        update ();

        /**
         * Return the boundary velocity as a function of position.
         */
        virtual
        Tensor<1,dim>
        boundary_velocity (const types::boundary_id boundary_indicator,
                           const Point<dim> &position) const = 0;

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
     * A class that manages all boundary velocity objects.
     *
     * @ingroup BoundaryVelocities
     */
    template <int dim>
    class Manager : public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Destructor. Made virtual since this class has virtual member
         * functions.
         */
        virtual ~Manager ();

        /**
         * A function that is called at the beginning of each time step and
         * calls the corresponding functions of all created plugins.
         *
         * The point of this function is to allow complex boundary velocity
         * models to do an initialization step once at the beginning of each
         * time step. An example would be a model that needs to call an
         * external program to compute the velocity change at a boundary.
         */
        virtual
        void
        update ();

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
                                    Interface<dim> *(*factory_function) ());


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
         */
        const std::map<types::boundary_id, std::pair<std::string,std::vector<std::string> > > &
        get_active_boundary_velocity_names () const;

        /**
         * Return pointers to all boundary velocity models
         * currently used in the computation, as specified in the input file.
         * The function returns a map between a boundary identifier and a vector
         * of shared pointers that represent the names of prescribed velocity
         * boundary models for this boundary. If there are no prescribed
         * boundary velocity plugins for a particular boundary this boundary
         * identifier will not appear in the map.
         */
        const std::map<types::boundary_id,std::vector<std::shared_ptr<BoundaryVelocity::Interface<dim> > > > &
        get_active_boundary_velocity_conditions () const;

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
        parse_parameters (ParameterHandler &prm);

        /**
         * Go through the list of all boundary velocity models that have been selected in
         * the input file (and are consequently currently active) and see if one
         * of them has the desired type specified by the template argument. If so,
         * return a pointer to it. If no boundary velocity model is active
         * that matches the given type, return a nullptr.
         */
        template <typename BoundaryVelocityType>
        DEAL_II_DEPRECATED
        BoundaryVelocityType *
        find_boundary_velocity_model () const;

        /**
         * Go through the list of all boundary velocity models that have been selected
         * in the input file (and are consequently currently active) and return
         * true if one of them has the desired type specified by the template
         * argument.
         */
        template <typename BoundaryVelocityType>
        bool
        has_matching_boundary_velocity_model () const;

        /**
         * Go through the list of all boundary velocity models that have been selected
         * in the input file (and are consequently currently active) and see
         * if one of them has the type specified by the template
         * argument or can be casted to that type. If so, return a reference
         * to it. If no boundary velocity model is active that matches the given type,
         * throw an exception.
         */
        template <typename BoundaryVelocityType>
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
         * A list of boundary velocity objects that have been requested in the
         * parameter file.
         */
        std::map<types::boundary_id,std::vector<std::shared_ptr<BoundaryVelocity::Interface<dim> > > > boundary_velocity_objects;

        /**
         * Map from boundary id to a pair
         * ("components", list of "velocity boundary type"),
         * where components is of the format "[x][y][z]" and the velocity type is
         * mapped to one of the plugins of velocity boundary conditions (e.g.
         * "function"). If the components string is empty, it is assumed the
         * plugins are used for all components.
         */
        std::map<types::boundary_id, std::pair<std::string,std::vector<std::string> > > boundary_velocity_indicators;

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
    template <typename BoundaryVelocityType>
    inline
    BoundaryVelocityType *
    Manager<dim>::find_boundary_velocity_model () const
    {
      for (typename std::map<types::boundary_id,std::vector<std::shared_ptr<BoundaryVelocity::Interface<dim> > > >::const_iterator
           boundary = boundary_velocity_objects.begin();
           boundary != boundary_velocity_objects.end(); ++boundary)
        for (typename std::vector<std::shared_ptr<BoundaryVelocity::Interface<dim>>>::const_iterator
             p = boundary->second.begin();
             p != boundary->second.end(); ++p)
          if (BoundaryVelocityType *x = dynamic_cast<BoundaryVelocityType *> ( (*p).get()) )
            return x;
      return nullptr;
    }


    template <int dim>
    template <typename BoundaryVelocityType>
    inline
    bool
    Manager<dim>::has_matching_boundary_velocity_model () const
    {
      for (typename std::map<types::boundary_id,std::vector<std::shared_ptr<BoundaryVelocity::Interface<dim> > > >::const_iterator
           boundary = boundary_velocity_objects.begin();
           boundary != boundary_velocity_objects.end(); ++boundary)
        for (typename std::vector<std::shared_ptr<BoundaryVelocity::Interface<dim>>>::const_iterator
             p = boundary->second.begin();
             p != boundary->second.end(); ++p)
          if (Plugins::plugin_type_matches<BoundaryVelocityType>(*(*p)))
            return true;
      return false;
    }


    template <int dim>
    template <typename BoundaryVelocityType>
    inline
    const BoundaryVelocityType &
    Manager<dim>::get_matching_boundary_velocity_model () const
    {
      AssertThrow(has_matching_boundary_velocity_model<BoundaryVelocityType> (),
                  ExcMessage("You asked BoundaryVelocity::Manager::get_boundary_velocity_model() for a "
                             "boundary velocity model of type <" + boost::core::demangle(typeid(BoundaryVelocityType).name()) + "> "
                             "that could not be found in the current model. Activate this "
                             "boundary velocity model in the input file."));

      typename std::map<types::boundary_id,std::vector<std::shared_ptr<BoundaryVelocity::Interface<dim> > > >::const_iterator boundary_velocity_model;
      for (typename std::map<types::boundary_id,std::vector<std::shared_ptr<BoundaryVelocity::Interface<dim> > > >::const_iterator
           boundary = boundary_velocity_objects.begin();
           boundary != boundary_velocity_objects.end(); ++boundary)
        for (typename std::vector<std::shared_ptr<BoundaryVelocity::Interface<dim>>>::const_iterator
             p = boundary->second.begin();
             p != boundary->second.end(); ++p)
          if (Plugins::plugin_type_matches<BoundaryVelocityType>(*(*p)))
            return Plugins::get_plugin_as_type<BoundaryVelocityType>(*(*p));

      // We will never get here, because we had the Assert above. Just to avoid warnings.
      return Plugins::get_plugin_as_type<BoundaryVelocityType>(*(*boundary_velocity_model));
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
    aspect::internal::Plugins::RegisterHelper<aspect::BoundaryVelocity::Interface<2>,classname<2> > \
    dummy_ ## classname ## _2d (&aspect::BoundaryVelocity::Manager<2>::register_boundary_velocity, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::BoundaryVelocity::Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&aspect::BoundaryVelocity::Manager<3>::register_boundary_velocity, \
                                name, description); \
  }
  }
}


#endif
