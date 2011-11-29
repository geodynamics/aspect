//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------
#ifndef __aspect__material_model_interface_h
#define __aspect__material_model_interface_h

#include <deal.II/base/point.h>
#include <deal.II/base/parameter_handler.h>

namespace aspect
{
  /**
   * A namespace in which we define everything that has to do with modeling
   * convecting material, including descriptions of material parameters such
   * as viscosities, densities, etc.
   *
   * @ingroup MaterialModels
   */
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * A base class for parameterizations of material models. Classes derived from
     * this class will need to implement functions that provide material parameters
     * such as the viscosity, density, etc, typically as a function of position,
     * temperature and pressure at that location.
     *
     * @ingroup MaterialModels
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
         * Return the viscosity $\eta$ of the model as a function of temperature,
         * pressure and position.
         */
        virtual double viscosity (const double      temperature,
                                  const double      pressure,
                                  const Point<dim> &position) const = 0;

        /**
         * Return a reference value typical of the viscosities that
         * appear in this model. This value is not actually used in the
         * material description itself, but is used in scaling variables
         * to the same numerical order of magnitude when solving linear
         * systems. Specifically, the reference viscosity appears in
         * the factor scaling the pressure against the velocity.
         */
        virtual double reference_viscosity () const = 0;

        /**
         * Return the specific heat $c_P$ of the model as a function of temperature,
         * pressure and position.
         */
        virtual double specific_heat (const double      temperature,
                                      const double      pressure,
                                      const Point<dim> &position) const = 0;

        /**
         * Return the thermal conductivity $k$ of the model as a function of temperature,
         * pressure and position.
         */
        virtual double thermal_conductivity (const double temperature,
                                             const double pressure,
                                             const Point<dim> &position) const = 0;

        /**
         * Return the density $\rho$ of the model as a function of temperature,
         * pressure and position.
         */
        virtual double density (const double      temperature,
                                const double      pressure,
                                const Point<dim> &position) const = 0;

        /**
         * Return the compressibility coefficient of the model as a function of temperature,
         * pressure and position.
         */
        virtual double compressibility (const double temperature,
                                        const double pressure,
                                        const Point<dim> &position) const = 0;

        /**
         * Return whether the model is compressible or not.  Incompressibility
         * does not necessarily imply that the density is constant; rather, it
         * may still depend on temperature or pressure. In the current
         * context, compressibility means whether we should solve the contuity
         * equation as $\nabla \cdot (\rho \mathbf u)=0$ (compressible Stokes)
         * or as $\nabla \cdot \mathbf{u}=0$ (incompressible Stokes).
         */
        virtual bool is_compressible () const = 0;

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
     * Register a material model so that it can be selected from the parameter file.
     *
     * @param name A string that identifies the material model
     * @param declare_parameters_function A pointer to a function that can be used to
     *   declare the parameters that this material model wants to read from input files.
     * @param factory_function A pointer to a function that can create an object of
     *   this material model.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    void
    register_material_model (const std::string &name,
                             void (*declare_parameters_function) (ParameterHandler &),
                             Interface<dim> * (*factory_function) ());

    /**
     * A function that given the name of a model returns a pointer to an object
     * that describes it. Ownership of the pointer is transferred to the caller.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    Interface<dim> *
    create_material_model (ParameterHandler &prm);


    /**
     * Declare the runtime parameters of the registered material models.
     *
     * @ingroup MaterialModels
     */
    void
    declare_parameters (ParameterHandler &prm);



    namespace internal
    {
      /**
       * An internal class that is used in the definition of the
       * ASPECT_REGISTER_MATERIAL_MODEL macro below. Given a name
       * and a classname, it registers the material model.
       */
      template <const char **name, class MaterialModelClass>
      struct MaterialModelHelper
      {
        MaterialModelHelper ()
        {
          register_material_model
          (*name,
           &MaterialModelClass::declare_parameters,
           &factory);
        }

        static
        Interface<deal_II_dimension> * factory ()
        {
          return new MaterialModelClass();
        }
      };
    }


    /**
     * Given a name and a classname for a material model, register it with
     * the functions that can declare their parameters and create these objects.
     *
     * @ingroup MaterialModels
     */
#define ASPECT_REGISTER_MATERIAL_MODEL(name,classname) \
  namespace ASPECT_REGISTER_MATERIAL_MODEL_ ## classname \
  { const char *local_name = name; \
    aspect::MaterialModel::internal::MaterialModelHelper<&local_name,classname<deal_II_dimension> > \
    dummy_ ## classname; }
  }
}


#endif
