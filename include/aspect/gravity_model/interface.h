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
     * @param description A text description of what this model
     * does and that will be listed in the documentation of
     * the parameter file.
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
                            const std::string &description,
                            void (*declare_parameters_function) (ParameterHandler &),
                            Interface<dim> *(*factory_function) ());

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


    /**
     * Given a class name, a name, and a description for the parameter file for a gravity model, register it with
     * the functions that can declare their parameters and create these objects.
     *
     * @ingroup GravityModels
     */
#define ASPECT_REGISTER_GRAVITY_MODEL(classname,name,description) \
  namespace ASPECT_REGISTER_GRAVITY_MODEL_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<Interface<deal_II_dimension>,classname<deal_II_dimension> > \
    dummy_ ## classname (&aspect::GravityModel::register_gravity_model<deal_II_dimension>, \
                         name, description); }
  }
}


#endif
