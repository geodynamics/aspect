/*
  Copyright (C) 2011, 2012 by the authors of the ASPECT code.

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


#ifndef __aspect__initial_conditions_interface_h
#define __aspect__initial_conditions_interface_h

#include <aspect/plugins.h>

#include <deal.II/base/point.h>
#include <deal.II/base/parameter_handler.h>

namespace aspect
{
  /**
   * A namespace in which we define everything that has to do with defining
   * the initial conditions.
   *
   * @ingroup InitialConditionsModels
   */
  namespace InitialConditions
  {
    using namespace dealii;

    /**
     * A base class for parameterizations of initial conditions.
     *
     * @ingroup InitialConditionsModels
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
         * Return the initial temperature as a function of position.
         */
        virtual
        double initial_temperature (const Point<dim> &position) const = 0;


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
     * Register an initial conditions model so that it can be selected from
     * the parameter file.
     *
     * @param name A string that identifies the initial conditions model
     * @param description A text description of what this model does and that
     * will be listed in the documentation of the parameter file.
     * @param declare_parameters_function A pointer to a function that can be
     * used to declare the parameters that this initial conditions model wants
     * to read from input files.
     * @param factory_function A pointer to a function that can create an
     * object of this initial conditions model.
     *
     * @ingroup InitialConditionsModels
     */
    template <int dim>
    void
    register_initial_conditions_model (const std::string &name,
                                       const std::string &description,
                                       void (*declare_parameters_function) (ParameterHandler &),
                                       Interface<dim> *(*factory_function) ());

    /**
     * A function that given the name of a model returns a pointer to an
     * object that describes it. Ownership of the pointer is transferred to
     * the caller.
     *
     * The model object returned is not yet initialized and has not read its
     * runtime parameters yet.
     *
     * @ingroup InitialConditionsModels
     */
    template <int dim>
    Interface<dim> *
    create_initial_conditions (ParameterHandler &prm);


    /**
     * Declare the runtime parameters of the registered initial conditions
     * models.
     *
     * @ingroup InitialConditionsModels
     */
    template <int dim>
    void
    declare_parameters (ParameterHandler &prm);



    /**
     * Given a class name, a name, and a description for the parameter file
     * for a initial conditions model, register it with the functions that can
     * declare their parameters and create these objects.
     *
     * @ingroup InitialConditionsModels
     */
#define ASPECT_REGISTER_INITIAL_CONDITIONS(classname,name,description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_INITIAL_CONDITIONS_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::InitialConditions::Interface<2>,classname<2> > \
    dummy_ ## classname ## _2d (&aspect::InitialConditions::register_initial_conditions_model<2>, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::InitialConditions::Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&aspect::InitialConditions::register_initial_conditions_model<3>, \
                                name, description); \
  }
  }
}


#endif
