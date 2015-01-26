/*
  Copyright (C) 2013 by the authors of the ASPECT code.

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


#ifndef __aspect__boundary_composition_interface_h
#define __aspect__boundary_composition_interface_h

#include <aspect/plugins.h>
#include <aspect/geometry_model/interface.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/distributed/tria.h>


namespace aspect
{
  /**
   * A namespace for the definition of things that have to do with describing
   * the boundary values for the composition.
   *
   * @ingroup BoundaryCompositions
   */
  namespace BoundaryComposition
  {
    using namespace dealii;

    /**
     * Base class for classes that describe composition boundary values.
     *
     * @ingroup BoundaryCompositions
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
         * A function that is called at the beginning of each time step. The
         * default implementation of the function does nothing, but derived
         * classes that need more elaborate setups for a given time step may
         * overload the function.
         *
         * The point of this function is to allow complex boundary composition
         * models to do an initialization step once at the beginning of each
         * time step. An example would be a model that needs to call an
         * external program to compute composition changes at sides.
         */
        virtual
        void
        update ();

        /**
         * Return the composition that is to hold at a particular location on
         * the boundary of the domain.
         *
         * @param geometry_model The geometry model that describes the domain.
         * This may be used to determine whether the boundary composition
         * model is implemented for this geometry.
         * @param boundary_indicator The boundary indicator of the part of the
         * boundary of the domain on which the point is located at which we
         * are requesting the composition.
         * @param location The location of the point at which we ask for the
         * composition.
         * @param compositional_field The number of the compositional field
         * for which we are requesting the composition.
         * @param compositional_field The index of the compositional field
         * between 0 and @p parameters.n_compositional_fields.
         */
        virtual
        double composition (const GeometryModel::Interface<dim> &geometry_model,
                            const unsigned int                   boundary_indicator,
                            const Point<dim>                    &location,
                            const unsigned int                   compositional_field) const = 0;

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
     * Register a boundary composition model so that it can be selected from
     * the parameter file.
     *
     * @param name A string that identifies the boundary composition model
     * @param description A text description of what this model does and that
     * will be listed in the documentation of the parameter file.
     * @param declare_parameters_function A pointer to a function that can be
     * used to declare the parameters that this geometry model wants to read
     * from input files.
     * @param factory_function A pointer to a function that can create an
     * object of this boundary composition model.
     *
     * @ingroup BoundaryCompositions
     */
    template <int dim>
    void
    register_boundary_composition (const std::string &name,
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
     * @ingroup BoundaryCompositions
     */
    template <int dim>
    Interface<dim> *
    create_boundary_composition (ParameterHandler &prm);


    /**
     * Declare the runtime parameters of the registered boundary composition
     * models.
     *
     * @ingroup BoundaryCompositions
     */
    template <int dim>
    void
    declare_parameters (ParameterHandler &prm);


    /**
     * Given a class name, a name, and a description for the parameter file
     * for a boundary composition model, register it with the functions that
     * can declare their parameters and create these objects.
     *
     * @ingroup BoundaryCompositions
     */
#define ASPECT_REGISTER_BOUNDARY_COMPOSITION_MODEL(classname, name, description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_BOUNDARY_COMPOSITION_MODEL_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::BoundaryComposition::Interface<2>,classname<2> > \
    dummy_ ## classname ## _2d (&aspect::BoundaryComposition::register_boundary_composition<2>, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::BoundaryComposition::Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&aspect::BoundaryComposition::register_boundary_composition<3>, \
                                name, description); \
  }
  }
}


#endif
