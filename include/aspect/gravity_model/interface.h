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


#ifndef __aspect__gravity_model_interface_h
#define __aspect__gravity_model_interface_h

#include <aspect/plugins.h>
#include <deal.II/base/point.h>
#include <deal.II/base/parameter_handler.h>

namespace aspect
{
  /**
   * A namespace in which we define everything that has to do with modeling
   * gravity.
   *
   * @ingroup GravityModels
   */
  namespace GravityModel
  {
    using namespace dealii;

    /**
     * A base class for parameterizations of gravity models.
     *
     * @ingroup GravityModels
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
         * Return the gravity vector as a function of position.
         */
        virtual Tensor<1,dim> gravity_vector (const Point<dim> &position) const = 0;

        /**
         * A function that is called at the beginning of each time step and
         * that allows the implementation to update internal data structures.
         * This is useful, for example, if you have a gravity model that
         * depends on time, or if you have a gravity model that depends on the
         * solution of the previous step.
         *
         * The default implementation of this function does nothing.
         */
        virtual void update();

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
     * Register a gravity model so that it can be selected from the parameter
     * file.
     *
     * @param name A string that identifies the gravity model
     * @param description A text description of what this model does and that
     * will be listed in the documentation of the parameter file.
     * @param declare_parameters_function A pointer to a function that can be
     * used to declare the parameters that this gravity model wants to read
     * from input files.
     * @param factory_function A pointer to a function that can create an
     * object of this gravity model.
     *
     * @ingroup GravityModels
     */
    template <int dim>
    void
    register_gravity_model (const std::string &name,
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
     * @ingroup GravityModels
     */
    template <int dim>
    Interface<dim> *
    create_gravity_model (ParameterHandler &prm);


    /**
     * Declare the runtime parameters of the registered gravity models.
     *
     * @ingroup GravityModels
     */
    template <int dim>
    void
    declare_parameters (ParameterHandler &prm);


    /**
     * Given a class name, a name, and a description for the parameter file
     * for a gravity model, register it with the functions that can declare
     * their parameters and create these objects.
     *
     * @ingroup GravityModels
     */
#define ASPECT_REGISTER_GRAVITY_MODEL(classname,name,description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_GRAVITY_MODEL_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::GravityModel::Interface<2>,classname<2> > \
    dummy_ ## classname ## _2d (&aspect::GravityModel::register_gravity_model<2>, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::GravityModel::Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&aspect::GravityModel::register_gravity_model<3>, \
                                name, description); \
  }
  }
}


#endif
