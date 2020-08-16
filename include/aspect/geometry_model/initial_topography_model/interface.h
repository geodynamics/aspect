/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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


#ifndef _aspect_geometry_model_initial_topography_model_interface_h
#define _aspect_geometry_model_initial_topography_model_interface_h

#include <aspect/plugins.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>

#include <set>


namespace aspect
{
  /**
   * A namespace for the definition of properties of the initial topography.
   * This includes mainly the storage and retrieval of the initial topography.
   * The retrieval is done through the value function, which requires a point
   * of size dim-1, and it returns a double which represents the elevation.
   *
   * @ingroup InitialTopographyModels
   */
  namespace InitialTopographyModel
  {
    using namespace dealii;

    /**
     * Base class for classes that describe particular initial topographies
     * for the domain.
     *
     * @ingroup InitialTopographyModels
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
        virtual void initialize ();

        /**
         * Return the value of the elevation at the given point.
         */
        virtual
        double value (const Point<dim-1> &p) const = 0;

        /**
         * Return the maximum value of the elevation.
         */
        virtual
        double max_topography () const = 0;

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
     * Register a initial topography model so that it can be selected from the parameter
     * file.
     *
     * @param name A string that identifies the initial topography model
     * @param description A text description of what this model does and that
     * will be listed in the documentation of the parameter file.
     * @param declare_parameters_function A pointer to a function that can be
     * used to declare the parameters that this initial topography model wants
     * to read from input files.
     * @param factory_function A pointer to a function that can create an
     * object of this initial topography model.
     *
     * @ingroup InitialTopographyModels
     */
    template <int dim>
    void
    register_initial_topography_model (const std::string &name,
                                       const std::string &description,
                                       void (*declare_parameters_function) (ParameterHandler &),
                                       Interface<dim> *(*factory_function) ());

    /**
     * A function that given the name of a model returns a pointer to an
     * object that describes it. Ownership of the pointer is transferred to
     * the caller.
     *
     * The initial topography model will also be asked to read its runtime
     * parameters already.
     *
     * @ingroup InitialTopographyModels
     */
    template <int dim>
    Interface<dim> *
    create_initial_topography_model (ParameterHandler &prm);


    /**
     * Declare the runtime parameters of the registered initial topography
     * models.
     *
     * @ingroup InitialTopographyModels
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
     * for a initial topography model, register it with the functions that
     * can declare their parameters and create these objects.
     *
     * @ingroup InitialTopographyModels
     */
#define ASPECT_REGISTER_INITIAL_TOPOGRAPHY_MODEL(classname,name,description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_INITIAL_TOPOGRAPHY_MODEL_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::InitialTopographyModel::Interface<2>,classname<2> > \
    dummy_ ## classname ## _2d (&aspect::InitialTopographyModel::register_initial_topography_model<2>, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::InitialTopographyModel::Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&aspect::InitialTopographyModel::register_initial_topography_model<3>, \
                                name, description); \
  }
  }
}


#endif
