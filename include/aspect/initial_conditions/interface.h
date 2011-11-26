//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------
#ifndef __aspect__initial_conditions_interface_h
#define __aspect__initial_conditions_interface_h

#include <aspect/geometry_model/interface.h>
#include <aspect/adiabatic_conditions.h>

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
         * Initialization function. Take references to the geometry model and the
         * adiabatic conditions and store them so that derived classes can access them.
         */
        void
        initialize (const GeometryModel::Interface<dim> &geometry_model,
                    const AdiabaticConditions<dim>      &adiabatic_conditions);

        /**
         * Return the initial temperature as a function of position.
         */
        virtual
        double initial_temperature (const Point<dim> &position) const = 0;

        /**
         * Declare the parameters this class takes through input files.
         * The default implementation of this function does not describe
         * any parameters. Consequently, derived classes do not have to
         * overload this function if they do not take any runtime parameters.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter
         * file. The default implementation of this function does not read
         * any parameters. Consequently, derived classes do not have to
         * overload this function if they do not take any runtime parameters.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

      protected:
        /**
         * Pointer to the geometry object in use.
         */
        const GeometryModel::Interface<dim> *geometry_model;

        /**
         * Pointer to an object that describes adiabatic conditions.
         */
        const AdiabaticConditions<dim>      *adiabatic_conditions;
    };




    /**
     * Register an initial conditions model so that it can be selected from the parameter file.
     *
     * @param name A string that identifies the initial conditions model
     * @param declare_parameters_function A pointer to a function that can be used to
     *   declare the parameters that this initial conditions model wants to read from input files.
     * @param factory_function A pointer to a function that can create an object of
     *   this initial conditions model.
     *
     * @ingroup InitialConditionsModels
     */
    template <int dim>
    void
    register_initial_conditions_model (const std::string &name,
                                       void (*declare_parameters_function) (ParameterHandler &),
                                       Interface<dim> * (*factory_function) ());

    /**
     * A function that given the name of a model returns a pointer to an object
     * that describes it. Ownership of the pointer is transferred to the caller.
     *
     * This function makes the newly created object read its parameters from the
     * input parameter object, and then initializes it with the given geometry
     * model and adiabatic conditions object.
     *
     * @ingroup InitialConditionsModels
     */
    template <int dim>
    Interface<dim> *
    create_initial_conditions (ParameterHandler &prm,
                               const GeometryModel::Interface<dim> &geometry_model,
                               const AdiabaticConditions<dim>      &adiabatic_conditions);


    /**
     * Declare the runtime parameters of the registered initial conditions models.
     *
     * @ingroup InitialConditionsModels
     */
    void
    declare_parameters (ParameterHandler &prm);



    namespace internal
    {
      /**
       * An internal class that is used in the definition of the
       * ASPECT_REGISTER_INITIAL_CONDITIONS macro below. Given a name
       * and a classname, it registers the initial conditions model.
       */
      template <const char **name, class InitialConditionsModelClass>
      struct InitialConditionsModelHelper
      {
        InitialConditionsModelHelper ()
        {
          register_initial_conditions_model
          (*name,
           &InitialConditionsModelClass::declare_parameters,
           &factory);
        }

        static
        Interface<deal_II_dimension> * factory ()
        {
          return new InitialConditionsModelClass();
        }
      };
    }


    /**
     * Given a name and a classname for a initial conditions model, register it with
     * the functions that can declare their parameters and create these objects.
     *
     * @ingroup InitialConditionsModels
     */
#define ASPECT_REGISTER_INITIAL_CONDITIONS(name,classname) \
  namespace ASPECT_REGISTER_INITIAL_CONDITIONS_ ## classname \
  { const char *local_name = name; \
    aspect::InitialConditions::internal::InitialConditionsModelHelper<&local_name,classname<deal_II_dimension> > \
    dummy_ ## classname; }
  }
}


#endif
