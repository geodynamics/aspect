//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------
#ifndef __aspect__gravity_model_interface_h
#define __aspect__gravity_model_interface_h

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
         * Return the gravity vector as a function of position.
         */
        virtual Tensor<1,dim> gravity_vector (const Point<dim> &position) const = 0;

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
    };




    /**
     * Register a gravity model so that it can be selected from the parameter file.
     *
     * @param name A string that identifies the gravity model
     * @param declare_parameters_function A pointer to a function that can be used to
     *   declare the parameters that this gravity model wants to read from input files.
     * @param factory_function A pointer to a function that can create an object of
     *   this gravity model.
     *
     * @ingroup GravityModels
     */
    template <int dim>
    void
    register_gravity_model (const std::string &name,
                            void (*declare_parameters_function) (ParameterHandler &),
                            Interface<dim> * (*factory_function) ());

    /**
     * A function that given the name of a model returns a pointer to an object
     * that describes it. Ownership of the pointer is transferred to the caller.
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
    void
    declare_parameters (ParameterHandler &prm);



    namespace internal
    {
      /**
       * An internal class that is used in the definition of the
       * ASPECT_REGISTER_GRAVITY_MODEL macro below. Given a name
       * and a classname, it registers the gravity model.
       */
      template <const char **name, class GravityModelClass>
      struct GravityModelHelper
      {
        GravityModelHelper ()
        {
          register_gravity_model
          (*name,
           &GravityModelClass::declare_parameters,
           &factory);
        }

        static
        Interface<deal_II_dimension> * factory ()
        {
          return new GravityModelClass();
        }
      };
    }


    /**
     * Given a name and a classname for a gravity model, register it with
     * the functions that can declare their parameters and create these objects.
     *
     * @ingroup GravityModels
     */
#define ASPECT_REGISTER_GRAVITY_MODEL(name,classname) \
  namespace ASPECT_REGISTER_GRAVITY_MODEL_ ## classname \
  { const char *local_name = name; \
    aspect::GravityModel::internal::GravityModelHelper<&local_name,classname<deal_II_dimension> > \
    dummy_ ## classname; }
  }
}


#endif
