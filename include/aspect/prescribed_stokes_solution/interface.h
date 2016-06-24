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



#ifndef __aspect__prescribed_stokes_solution_interface_h
#define __aspect__prescribed_stokes_solution_interface_h

#include <aspect/plugins.h>

#include <deal.II/base/point.h>
#include <deal.II/base/parameter_handler.h>

#include <deal.II/base/function.h>

namespace aspect
{
  /**
   * A namespace in which we define everything that has to do with defining
   * the prescribed Stokes solution.
   *
   * @ingroup PrescribedStokesSolution
   */
  namespace PrescribedStokesSolution
  {
    using namespace dealii;

    /**
     * This plugin allows the user to prescribe a Stokes solution and can be
     * thought of as velocity and pressure's equivalent of the initial
     * conditions plugin.
     *
     * Note: Only used if solver type  is "Advection only".
     *
     * @ingroup PrescribedStokesSolution
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
         * The point of this function is to allow complex prescribed Stokes
         * solutions to do an initialization step once at the beginning of each
         * time step. An example would be a time-dependent prescribed velocity.
         */
        virtual
        void
        update ();

        /**
         * Given a position @p p, fill in desired velocity and pressure at
         * that point into @p value, which will have dim+1 components. In @p
         * value, the velocity components come first, followed by the pressure
         * component.
         */
        virtual
        void stokes_solution (const Point<dim> &p, Vector<double> &value) const = 0;

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
     * Register initial prescribed Stokes solution model so that it can be
     * selected from the parameter file.
     *
     * @param name A string that identifies the prescribed Stokes solution
     * model
     * @param description A text description of what this model does and that
     * will be listed in the documentation of the parameter file.
     * @param declare_parameters_function A pointer to a function that can be
     * used to declare the parameters that this prescribed Stokes solution
     * model wants to read from input files.
     * @param factory_function A pointer to a function that can create an
     * object of this prescribed Stokes solution model.
     *
     * @ingroup PrescribedStokesSolution
     */
    template <int dim>
    void
    register_prescribed_stokes_solution_model (const std::string &name,
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
     * @ingroup PrescribedStokesSolution
     */
    template <int dim>
    Interface<dim> *
    create_prescribed_stokes_solution (ParameterHandler &prm);


    /**
     * Declare the runtime parameters of the registered prescribed Stokes
     * solution models.
     *
     * @ingroup PrescribedStokesSolution
     */
    template <int dim>
    void
    declare_parameters (ParameterHandler &prm);



    /**
     * Given a class name, a name, and a description for the parameter file
     * for a prescribed Stokes solution model, register it with the functions
     * that can declare their parameters and create these objects.
     *
     * @ingroup PrescribedStokesSolution
     */
#define ASPECT_REGISTER_PRESCRIBED_STOKES_SOLUTION(classname,name,description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_PRESCRIBED_STOKES_SOLUTION_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::PrescribedStokesSolution::Interface<2>,classname<2> > \
    dummy_ ## classname ## _2d (&aspect::PrescribedStokesSolution::register_prescribed_stokes_solution_model<2>, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::PrescribedStokesSolution::Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&aspect::PrescribedStokesSolution::register_prescribed_stokes_solution_model<3>, \
                                name, description); \
  }
  }
}


#endif
