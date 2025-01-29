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


#ifndef _aspect_boundary_traction_interface_h
#define _aspect_boundary_traction_interface_h

#include <aspect/plugins.h>
#include <aspect/geometry_model/interface.h>

#include <deal.II/base/point.h>
#include <deal.II/base/parameter_handler.h>

namespace aspect
{
  /**
   * A namespace in which we define everything that has to do with defining
   * traction boundary conditions for the Stokes equations.
   *
   * @ingroup BoundaryTractions
   */
  namespace BoundaryTraction
  {
    /**
     * A base class for parameterizations of traction boundary conditions.
     *
     * @ingroup BoundaryTractions
     */
    template <int dim>
    class Interface : public Plugins::InterfaceBase
    {
      public:
        /**
         * Return the traction that is to hold at a particular position on
         * the boundary of the domain.
         *
         * @param boundary_indicator The boundary indicator of the part of the
         * boundary of the domain on which the point is located at which we
         * are requesting the traction.
         * @param position The position of the point at which we ask for the
         * traction.
         * @param normal_vector The (outward) normal vector to the boundary
         * of the domain.
         *
         * @return Boundary traction at position @p position.
         */
        virtual
        Tensor<1,dim>
        boundary_traction (const types::boundary_id boundary_indicator,
                           const Point<dim> &position,
                           const Tensor<1,dim> &normal_vector) const = 0;
    };

    template <int dim>
    class Manager : public SimulatorAccess<dim>
    {
      public:
        /**
         * Destructor. Made virtual since this class has virtual member
         * functions.
         */
        ~Manager () override;

        /**
         * A function that is called at the beginning of each time step and
         * calls the corresponding functions of all created plugins.
         *
         * The point of this function is to allow complex boundary traction
         * models to do an initialization step once at the beginning of each
         * time step. An example would be a model that needs to call an
         * external program to compute the traction change at a boundary.
         */
        virtual
        void
        update ();

        /**
         * A function that calls the boundary_traction functions of all the
         * individual boundary traction objects and uses the stored operators
         * to combine them.
         */
        Tensor<1,dim>
        boundary_traction (const types::boundary_id boundary_indicator,
                           const Point<dim> &position,
                           const Tensor<1,dim> &normal_vector) const;

        /**
         * Return the names of all prescribed boundary traction models currently
         * used in the computation as specified in the input file. The function
         * returns a map between a boundary identifier and a pair. The
         * first part of the pair is a string that represents the prescribed
         * traction components on this boundary (e.g. y, xz, or xyz) and the
         * second part is a vector of strings that represent the names of
         * boundary traction plugins for this boundary.
         * If there are no prescribed boundary traction plugins
         * for a particular boundary, this boundary identifier will not appear
         * in the map.
         */
        const std::map<types::boundary_id, std::pair<std::string,std::vector<std::string>>> &
        get_active_boundary_traction_names () const;

        /**
         * Return pointers to all boundary traction models
         * currently used in the computation, as specified in the input file.
         * The function returns a map between a boundary identifier and a vector
         * of unique pointers that represent the names of prescribed traction
         * boundary models for this boundary. If there are no prescribed
         * boundary traction plugins for a particular boundary this boundary
         * identifier will not appear in the map.
         */
        const std::map<types::boundary_id,std::vector<std::unique_ptr<BoundaryTraction::Interface<dim>>>> &
        get_active_boundary_traction_conditions () const;

        /**
         * Declare the parameters of all known boundary traction plugins, as
         * well as the ones this class has itself.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         * This determines which boundary traction objects will be created;
         * then let these objects read their parameters as well.
         */
        void
        parse_parameters (ParameterHandler &prm);

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
        DeclException1 (ExcBoundaryTractionNameNotFound,
                        std::string,
                        << "Could not find entry <"
                        << arg1
                        << "> among the names of registered boundary traction objects.");

        /**
         * Register a traction boundary conditions model so that it can be
         * selected from the parameter file.
         *
         * @param name A string that identifies the traction boundary conditions
         * model
         * @param description A text description of what this model does and that
         * will be listed in the documentation of the parameter file.
         * @param declare_parameters_function A pointer to a function that can be
         * used to declare the parameters that this traction boundary conditions
         * model wants to read from input files.
         * @param factory_function A pointer to a function that can create an
         * object of this traction boundary conditions model.
         *
         * @ingroup BoundaryTractions
         */
        static
        void
        register_boundary_traction (const std::string &name,
                                    const std::string &description,
                                    void (*declare_parameters_function) (ParameterHandler &),
                                    std::unique_ptr<Interface<dim>> (*factory_function) ());

      private:
        /**
         * A list of boundary traction objects that have been requested in the
         * parameter file.
         */
        std::map<types::boundary_id,std::vector<std::unique_ptr<BoundaryTraction::Interface<dim>>>> boundary_traction_objects;

        /**
         * Map from boundary id to a pair
         * ("components", list of "traction boundary type"),
         * where components is of the format "[x][y][z]" and the traction type is
         * mapped to one of the plugins of traction boundary conditions (e.g.
         * "function"). If the components string is empty, it is assumed the
         * plugins are used for all components.
         */
        std::map<types::boundary_id, std::pair<std::string,std::vector<std::string>>> boundary_traction_indicators;

    };

    /**
     * A function that given the name of a model returns a pointer to an
     * object that describes it. Ownership of the pointer is transferred to
     * the caller.
     *
     * The model object returned is not yet initialized and has not
     * read its runtime parameters yet.
     *
     * @ingroup BoundaryTractions
     */
    template <int dim>
    std::unique_ptr<Interface<dim>>
    create_boundary_traction (const std::string &name);

    /**
     * Return a list of names of all implemented boundary traction models,
     * separated by '|' so that it can be used in an object of type
     * Patterns::Selection.
     */
    template <int dim>
    std::string
    get_names ();

    /**
     * Declare the runtime parameters of the registered traction boundary
     * conditions models.
     *
     * @ingroup BoundaryTractions
     */
    template <int dim>
    void
    declare_parameters (ParameterHandler &prm);


    /**
     * For the current plugin subsystem, write a connection graph of all of the
     * plugins we know about, in the format that the
     * programs dot and neato understand. This allows for a visualization of
     * how all of the plugins that ASPECT knows about are interconnected, and
     * connect to other parts of the ASPECT code.
     *
     * @param output_stream The stream to write the output to.
     */
    template <int dim>
    void
    write_plugin_graph (std::ostream &output_stream);



    /**
     * Given a class name, a name, and a description for the parameter file
     * for a traction boundary conditions model, register it with the
     * functions that can declare their parameters and create these objects.
     *
     * @ingroup BoundaryTractions
     */
#define ASPECT_REGISTER_BOUNDARY_TRACTION_MODEL(classname,name,description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_BOUNDARY_TRACTION_MODEL_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::BoundaryTraction::Interface<2>,classname<2>> \
    dummy_ ## classname ## _2d (&aspect::BoundaryTraction::Manager<2>::register_boundary_traction, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::BoundaryTraction::Interface<3>,classname<3>> \
    dummy_ ## classname ## _3d (&aspect::BoundaryTraction::Manager<3>::register_boundary_traction, \
                                name, description); \
  }
  }
}


#endif
