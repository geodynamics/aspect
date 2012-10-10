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
/*  $Id$  */


#ifndef __aspect__material_model_interface_h
#define __aspect__material_model_interface_h

#include <aspect/plugins.h>
#include <deal.II/base/point.h>
#include <deal.II/base/symmetric_tensor.h>
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
     * An namespace whose enum members are used in querying the
     * nonlinear dependence of physical parameters on other quantities.
     */
    namespace NonlinearDependence
    {
      enum Dependence
      {
        temperature,
        pressure,
        strain_rate
      };
    };


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
         * Called at the beginning of each time step and allows the material model
         * to update internal data structures.
         */
        virtual void update();

        /**
         * @name Physical parameters used in the basic equations
         * @{
         */
        /**
         * Return the viscosity $\eta$ of the model as a function of temperature,
         * pressure, strain rate, and position.
        *
        * @note The strain rate given as the third argument of this function
        * is computed as $\varepsilon(\mathbf u)=\frac 12 (\nabla \mathbf u +
        * \nabla \mathbf u^T)$, regardless of whether the model is
        * compressible or not. This is relevant since in some other contexts,
        * the strain rate in the compressible case is computed as
        * $\varepsilon(\mathbf u)=\frac 12 (\nabla \mathbf u +
        * \nabla \mathbf u^T) - \frac 13 \nabla \cdot \mathbf u \mathbf 1$.
         */
        virtual double viscosity (const double                  temperature,
                                  const double                  pressure,
                                  const SymmetricTensor<2,dim> &strain_rate,
                                  const Point<dim>             &position) const = 0;

        /**
         * Return the viscosity ratio between disclocation creep and diffusion creep
         * in the case of composite rheology
         */
        virtual double viscosity_ratio (const double      temperature,
                                        const double      pressure,
                                        const SymmetricTensor<2,dim> &strainrate,
                                        const Point<dim> &position) const;

        /**
         * Return the density $\rho$ of the model as a function of temperature,
         * pressure and position.
         */
        virtual double density (const double      temperature,
                                const double      pressure,
                                const Point<dim> &position) const = 0;

        /**
         * Return the compressibility coefficient
         * $\frac 1\rho \frac{\partial\rho}{\partial p}$ of the model as a
         * function of temperature, pressure and position.
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
         * Return the thermal expansion coefficient $\alpha$ of the model,
         * possibly as a function of temperature, pressure and position.
        * The thermal expansion coefficient is defined as
        * $\alpha=-\frac{1}{\rho} \frac{d\rho}{dT}$. Since the density
        * <i>decreases</i> with temperature for almost all models,
        * $\alpha$ is usually positive.
        *
        * This function has a default implementation that computes $\alpha$
        * through its definition above, using the density() and density_derivative()
        * functions.
         */
        virtual double thermal_expansion_coefficient (const double      temperature,
                                                      const double      pressure,
                                                      const Point<dim> &position) const;

        /**
         * Return the thermal conductivity $k$ of the model as a function of temperature,
         * pressure and position. The units of $k$ are $\textrm{W} / \textrm{m} / \textrm{K}$
         * in 3d, and $\textrm{W} / \textrm{K}$ in 2d. This is easily see by considering that
         * $k$ is the heat flux density (i.e., Watts per unit area perpendicular to the heat
         * flux direction) per unit temperature gradient (i.e., Kelvin per meter). The unit
         * area has units $m^2$ in 3d, but only $m$ in 2d, yielding the stated units for $k$.
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
        * Return true if the viscosity() function returns something that
        * may depend on the variable identifies by the argument.
        */
        virtual bool
        viscosity_depends_on (const NonlinearDependence::Dependence dependence) const = 0;

        /**
        * Return true if the density() function returns something that
        * may depend on the variable identifies by the argument.
        */
        virtual bool
        density_depends_on (const NonlinearDependence::Dependence dependence) const = 0;

        /**
        * Return true if the compressibility() function returns something that
        * may depend on the variable identifies by the argument.
        *
        * This function must return false for all possible arguments if the
        * is_compressible() function returns false.
        */
        virtual bool
        compressibility_depends_on (const NonlinearDependence::Dependence dependence) const = 0;

        /**
        * Return true if the specific_heat() function returns something that
        * may depend on the variable identifies by the argument.
        */
        virtual bool
        specific_heat_depends_on (const NonlinearDependence::Dependence dependence) const = 0;

        /**
        * Return true if the thermal_conductivity() function returns something that
        * may depend on the variable identifies by the argument.
        */
        virtual bool
        thermal_conductivity_depends_on (const NonlinearDependence::Dependence dependence) const = 0;

        /**
         * Return whether the model is compressible or not.  Incompressibility
         * does not necessarily imply that the density is constant; rather, it
         * may still depend on temperature or pressure. In the current
         * context, compressibility means whether we should solve the continuity
         * equation as $\nabla \cdot (\rho \mathbf u)=0$ (compressible Stokes)
         * or as $\nabla \cdot \mathbf{u}=0$ (incompressible Stokes).
         */
        virtual bool is_compressible () const = 0;
        /**
         * @}
         */


        /**
         * @name Partial derivatives of physical parameters
         * @{
         */

        /**
        * Return the partial derivative of the viscosity function on the
        * variable indicates as last argument.
        *
        * The default implementation of this function returns zero
        * provided viscosity_depends_on() returns false for the given
        * dependence and throws an exception otherwise.
        */
        virtual double
        viscosity_derivative (const double              temperature,
                              const double              pressure,
                              const Point<dim>         &position,
                              const NonlinearDependence::Dependence dependence) const;

        /**
        * Return the partial derivative of the density function on the
        * variable indicates as last argument.
        *
        * The default implementation of this function returns zero
        * provided density_depends_on() returns false for the given
        * dependence and throws an exception otherwise.
        */
        virtual double
        density_derivative (const double              temperature,
                            const double              pressure,
                            const Point<dim>         &position,
                            const NonlinearDependence::Dependence dependence) const;

        /**
        * Return the partial derivative of the compressibility function on the
        * variable indicates as last argument.
        *
        * The default implementation of this function returns zero
        * provided compressibility_depends_on() returns false for the given
        * dependence and throws an exception otherwise.
        */
        virtual double
        compressibility_derivative (const double              temperature,
                                    const double              pressure,
                                    const Point<dim>         &position,
                                    const NonlinearDependence::Dependence dependence) const;

        /**
        * Return the partial derivative of the specific heat function on the
        * variable indicates as last argument.
        *
        * The default implementation of this function returns zero
        * provided specific_heat_depends_on() returns false for the given
        * dependence and throws an exception otherwise.
        */
        virtual double
        specific_heat_derivative (const double              temperature,
                                  const double              pressure,
                                  const Point<dim>         &position,
                                  const NonlinearDependence::Dependence dependence) const;

        /**
        * Return the partial derivative of the thermal conductivity
        * function on the variable indicates as last argument.
        *
        * The default implementation of this function returns zero
        * provided thermal_conductivity_depends_on() returns false
        * for the given dependence and throws an exception otherwise.
        */
        virtual double
        thermal_conductivity_derivative (const double              temperature,
                                         const double              pressure,
                                         const Point<dim>         &position,
                                         const NonlinearDependence::Dependence dependence) const;
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
         * Return a reference value for the thermal expansion
         * coefficient $\alpha$. See the thermal_expansion_coefficient()
         * function for a definition of $\alpha$.
         */
        virtual double reference_thermal_expansion_coefficient () const = 0;
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
                    const double      pressure,
                    const Point<dim> &position) const;

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
                    const double      pressure,
                    const Point<dim> &position) const;
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

        struct MaterialModelInputs
        {
          MaterialModelInputs(unsigned int n_points, unsigned int n_comp);

          std::vector<Point<dim>> position;
          std::vector<double> temperature;
          std::vector<double> pressure;
          std::vector<std::vector<double>> composition;
          std::vector<SymmetricTensor<2,dim>> strain_rate;
       };

        struct MaterialModelOutputs
        {
          MaterialModelOutputs(unsigned int n_points);

          std::vector<double> viscosities;
          std::vector<double> densities;
          std::vector<double> thermal_expansion_coefficients;
          std::vector<double> seismic_Vp;
          std::vector<double> seismic_Vs;
          std::vector<double> specific_heat;
          std::vector<double> thermal_conductivities;
          std::vector<double> compressibilities;
          std::vector<int> thermodynamic_phases;
          bool is_compressible;
        };

        virtual void compute_parameters(MaterialModelInputs & in, MaterialModelOutputs & out);

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
                             Interface<dim> *(*factory_function) ());

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
    template <int dim>
    void
    declare_parameters (ParameterHandler &prm);



    /**
     * Given a class name, a name, and a description for the parameter file for a material model, register it with
     * the functions that can declare their parameters and create these objects.
     *
     * @ingroup MaterialModels
     */
#define ASPECT_REGISTER_MATERIAL_MODEL(classname,name,description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_MATERIAL_MODEL_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<Interface<2>,classname<2> > \
    dummy_ ## classname ## _2d (&aspect::MaterialModel::register_material_model<2>, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&aspect::MaterialModel::register_material_model<3>, \
                                name, description); \
  }
  }
}


#endif
