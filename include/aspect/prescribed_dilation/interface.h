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


#ifndef _aspect_prescibed_dilation_interface_h
#define _aspect_prescibed_dilation_interface_h

#include <aspect/plugins.h>
#include <deal.II/base/point.h>
#include <deal.II/base/parameter_handler.h>

namespace aspect
{
  /**
   * A namespace in which we define everything that has to do with
   * prescribing dilation.
   *
   * @ingroup PrescribedDilation
   */
  namespace PrescribedDilation
  {
    /**
     * A base class for parameterizations of prescribed dilation models.
     *
     * @ingroup PrescribedDilation
     */
    template <int dim>
    class Interface : public Plugins::InterfaceBase
    {
      public:
        /**
         * Return the dilation as a function of position.
         */
        virtual double dilation (const Point<dim> &position) const = 0;
    };




    /**
     * Register a prescribed dilation model so that it can be selected from the parameter
     * file.
     *
     * @param name A string that identifies the prescribed dilation model
     * @param description A text description of what this model does and that
     * will be listed in the documentation of the parameter file.
     * @param declare_parameters_function A pointer to a function that can be
     * used to declare the parameters that this prescribed dilation model wants to read
     * from input files.
     * @param factory_function A pointer to a function that can create an
     * object of this prescribed dilation model.
     *
     * @ingroup PrescribedDilation
     */
    template <int dim>
    void
    register_prescribed_dilation_model (const std::string &name,
                                        const std::string &description,
                                        void (*declare_parameters_function) (ParameterHandler &),
                                        std::unique_ptr<Interface<dim>> (*factory_function) ());

    /**
     * A function that given the name of a model returns a pointer to an
     * object that describes it. Ownership of the pointer is transferred to
     * the caller.
     *
     * The model object returned is not yet initialized and has not read its
     * runtime parameters yet.
     *
     * @ingroup PrescribedDilation
     */
    template <int dim>
    std::unique_ptr<Interface<dim>>
    create_prescribed_dilation_model (ParameterHandler &prm);


    /**
     * Declare the runtime parameters of the registered prescribed dilation models.
     *
     * @ingroup PrescribedDilation
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
     * for a prescribed dilation model, register it with the functions that can declare
     * their parameters and create these objects.
     *
     * @ingroup PrescribedDilation
     */
#define ASPECT_REGISTER_PRESCRIBED_DILATION_MODEL(classname,name,description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_PRESCRIBED_DILATION_MODEL_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::PrescribedDilation::Interface<2>,classname<2>> \
    dummy_ ## classname ## _2d (&aspect::PrescribedDilation::register_prescribed_dilation_model<2>, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::PrescribedDilation::Interface<3>,classname<3>> \
    dummy_ ## classname ## _3d (&aspect::PrescribedDilation::register_prescribed_dilation_model<3>, \
                                name, description); \
  }
  }
}


#endif
