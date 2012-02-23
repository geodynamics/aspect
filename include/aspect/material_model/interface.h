//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011, 2012 by the authors of the ASPECT code
//
//-------------------------------------------------------------
#ifndef __aspect__material_model_interface_h
#define __aspect__material_model_interface_h

#include <aspect/plugins.h>
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
         * @name Physical parameters used in the basic equations
         * @{
         */
        /**
         * Return the viscosity $\eta$ of the model as a function of temperature,
         * pressure and position.
         */
        virtual double viscosity (const double      temperature,
                                  const double      pressure,
                                  const Point<dim> &position) const = 0;

        /**
         * Return the density $\rho$ of the model as a function of temperature,
         * pressure and position.
         */
        virtual double density (const double      temperature,
                                const double      pressure,
                                const Point<dim> &position) const = 0;

        /**
         * Return the compressibility coefficient
         * $\frac{\partial\rho}{\partial p}$ of the model as a
         * function of temperature, pressure and position.
         *
         * A more intuitive way would be to say that this function is
         * in fact the negative compressibility, since it will return
         * a negative value for realistic materials.
         */
        virtual double compressibility (const double temperature,
                                        const double pressure,
                                        const Point<dim> &position) const = 0;

        /**
         * Return the specific heat $C_p$ of the model as a function of temperature,
         * pressure and position.
         */
        virtual double specific_heat (const double      temperature,
                                      const double      pressure,
                                      const Point<dim> &position) const = 0;

        /**
         * Return the thermal conductivity $k$ of the model as a function of temperature,
         * pressure and position. The units of $k$ are $\textrm{W} / \textrm{m} / \textrm{K}$.
        *
        * Note that the thermal <i>conductivity</i> $k$ is related to the thermal
        * <i>diffusivity</i> $\kappa$ as $k = \kappa \rho c_p$. In essence, the conductivity
        * relates to the question of how thermal energy diffuses whereas the diffusivity
        * relates to the question of how the temperature diffuses. $\kappa$ has units
        * $\textrm{m}^2/\textrm{s}$.
         */
        virtual double thermal_conductivity (const double temperature,
                                             const double pressure,
                                             const Point<dim> &position) const = 0;
        /**
         * @}
         */

        /**
         * @name Qualitative properties one can ask a material model
         * @{
         */
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
         * @}
         */

        /**
         * @name Reference quantities
         * @{
         */
        /**
         * Return a reference value typical of the viscosities that
         * appear in this model. This value is not actually used in the
         * material description itself, but is used in scaling variables
         * to the same numerical order of magnitude when solving linear
         * systems. Specifically, the reference viscosity appears in
         * the factor scaling the pressure against the velocity. It is
        * also used in computing dimension-less quantities.
         */
        virtual double reference_viscosity () const = 0;

        /**
        * Return the reference density $\rho$. Like the value returned
        * by reference_viscosity(), this value is not actually used in
        * computations but only in postprocessing such as when computing
        * dimension-less quantities.
        */
        virtual double reference_density () const = 0;
        /**
         * @}
         */

        /**
         * @name Auxiliary material properties used for postprocessing
         * @{
         */
        /**
         * Return the p-wave seismic velocity Vp of the model as a
         * function of temperature and pressure.
        *
        * This function is only called in postprocessing. Derived classes do
         * not need to implement it if no useful information is known to
         * compute this quantity, in which case graphical output will simply
         * show an uninformative field of constant value. By default this
         * function returns -1 to indicate that no useful value is
         * implemented.
         */
        virtual
        double
        seismic_Vp (const double      temperature,
                    const double      pressure) const;

        /**
         * Return the s-wave seismic velocity Vs of the model as a
         * function of temperature and pressure.
        *
        * This function is only called in postprocessing. Derived classes do
         * not need to implement it if no useful information is known to
         * compute this quantity, in which case graphical output will simply
         * show an uninformative field of constant value. By default this
         * function returns -1 to indicate that no useful value is
         * implemented.
         */
        virtual
        double
        seismic_Vs (const double      temperature,
                    const double      pressure) const;
        /**
         * Return the Phase number of the model as a function of
         * temperature and pressure.
        *
        * This function is only called in postprocessing. Derived classes do
         * not need to implement it if no useful information is known to
         * compute this quantity, in which case graphical output will simply
         * show an uninformative field of constant value. By default this
         * function returns 0 to indicate everything is part of the same phase
         */
        virtual
        unsigned int
        thermodynamic_phase (const double      temperature,
                             const double      pressure) const;
        /**
         * @}
         */

        /**
         * @name Functions used in dealing with run-time parameters
         * @{
         */
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
        /**
         * @}
         */
    };




    /**
     * Register a material model so that it can be selected from the parameter file.
     *
     * @param name A string that identifies the material model
     * @param description A text description of what this model
     * does and that will be listed in the documentation of
     * the parameter file.
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
                             const std::string &description,
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



    /**
     * Given a class name, a name, and a description for the parameter file for a material model, register it with
     * the functions that can declare their parameters and create these objects.
     *
     * @ingroup MaterialModels
     */
#define ASPECT_REGISTER_MATERIAL_MODEL(classname,name,description) \
  namespace ASPECT_REGISTER_MATERIAL_MODEL_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<Interface<deal_II_dimension>,classname<deal_II_dimension> > \
    dummy_ ## classname (&aspect::MaterialModel::register_material_model<deal_II_dimension>, \
                         name, description); }
  }
}


#endif
