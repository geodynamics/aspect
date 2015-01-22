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


#ifndef __aspect__boundary_temperature_interface_h
#define __aspect__boundary_temperature_interface_h

#include <aspect/plugins.h>
#include <aspect/geometry_model/interface.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/distributed/tria.h>


namespace aspect
{
  /**
   * A namespace for the definition of things that have to do with describing
   * the boundary values for the temperature.
   *
   * @ingroup BoundaryTemperatures
   */
  namespace BoundaryTemperature
  {
    using namespace dealii;

    /**
     * Base class for classes that describe temperature boundary values.
     *
     * @ingroup BoundaryTemperatures
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
         * Return the temperature that is to hold at a particular location on
         * the boundary of the domain.
         *
         * @param geometry_model The geometry model that describes the domain.
         * This may be used to determine whether the boundary temperature
         * model is implemented for this geometry.
         * @param boundary_indicator The boundary indicator of the part of the
         * boundary of the domain on which the point is located at which we
         * are requesting the temperature.
         * @param location The location of the point at which we ask for the
         * temperature.
         */
        virtual
        double temperature (const GeometryModel::Interface<dim> &geometry_model,
                            const unsigned int                   boundary_indicator,
                            const Point<dim>                    &location) const = 0;

        /**
         * Return the minimal the temperature on that part of the boundary on
         * which Dirichlet conditions are posed.
         *
         * This value is used in computing dimensionless numbers such as the
         * Nusselt number indicating heat flux.
         */
        virtual
        double minimal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids =
                                      std::set<types::boundary_id>()) const = 0;

        /**
         * Return the maximal the temperature on that part of the boundary on
         * which Dirichlet conditions are posed.
         *
         * This value is used in computing dimensionless numbers such as the
         * Nusselt number indicating heat flux.
         */
        virtual
        double maximal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids =
                                      std::set<types::boundary_id>()) const = 0;

        /**
         * A function that is called at the beginning of each time step. The
         * default implementation of the function does nothing, but derived
         * classes that need more elaborate setups for a given time step may
         * overload the function.
         *
         * The point of this function is to allow complex boundary temperature
         * models to do an initialization step once at the beginning of each
         * time step. An example would be a model that needs to call an
         * external program to compute temperature change at bottom.
         */
        virtual
        void
        update ();

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
     * Register a boundary temperature model so that it can be selected from
     * the parameter file.
     *
     * @param name A string that identifies the boundary temperature model
     * @param description A text description of what this model does and that
     * will be listed in the documentation of the parameter file.
     * @param declare_parameters_function A pointer to a function that can be
     * used to declare the parameters that this geometry model wants to read
     * from input files.
     * @param factory_function A pointer to a function that can create an
     * object of this boundary temperature model.
     *
     * @ingroup BoundaryTemperatures
     */
    template <int dim>
    void
    register_boundary_temperature (const std::string &name,
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
     * @ingroup BoundaryTemperatures
     */
    template <int dim>
    Interface<dim> *
    create_boundary_temperature (ParameterHandler &prm);


    /**
     * Declare the runtime parameters of the registered boundary temperature
     * models.
     *
     * @ingroup BoundaryTemperatures
     */
    template <int dim>
    void
    declare_parameters (ParameterHandler &prm);


    /**
     * Given a class name, a name, and a description for the parameter file
     * for a boundary temperature model, register it with the functions that
     * can declare their parameters and create these objects.
     *
     * @ingroup BoundaryTemperatures
     */
#define ASPECT_REGISTER_BOUNDARY_TEMPERATURE_MODEL(classname, name, description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_BOUNDARY_TEMPERATURE_MODEL_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::BoundaryTemperature::Interface<2>,classname<2> > \
    dummy_ ## classname ## _2d (&aspect::BoundaryTemperature::register_boundary_temperature<2>, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::BoundaryTemperature::Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&aspect::BoundaryTemperature::register_boundary_temperature<3>, \
                                name, description); \
  }
  }
}


#endif
