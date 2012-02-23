//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011, 2012 by the authors of the ASPECT code
//
//-------------------------------------------------------------
#ifndef __aspect__model_table_h
#define __aspect__model_table_h

#include <aspect/material_model/interface.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * A variable viscosity material model that reads the essential values of coefficients from
     * tables in input files.
     *
     * This model is based on the paper Steinberger/Calderwood 2006:
     * "Models of large-scale viscous flow in the Earth's mantle with
     * contraints from mineral physics and surface observations"
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class Steinberger: public MaterialModel::Interface<dim>
    {
      public:
        virtual double viscosity (const double temperature,
                                  const double pressure,
                                  const Point<dim> &position) const;

        virtual double reference_viscosity () const;

        virtual double reference_density () const;

        virtual double specific_heat (const double temperature,
                                      const double pressure,
                                      const Point<dim> &position) const;

        virtual double thermal_conductivity (const double temperature,
                                             const double pressure,
                                             const Point<dim> &position) const;

        virtual double density (const double temperature,
                                const double pressure,
                                const Point<dim> &position) const;

        virtual double compressibility (const double temperature,
                                        const double pressure,
                                        const Point<dim> &position) const;

        /**
         * @name Qualitative properties one can ask a material model
         * @{
         */

        /**
	 * Return true if the viscosity() function returns something that
	 * may depend on the variable identifies by the argument.
	 */
	virtual bool
	viscosity_depends_on (const NonlinearDependence::Dependence dependence) const;

        /**
	 * Return true if the density() function returns something that
	 * may depend on the variable identifies by the argument.
	 */
	virtual bool
	density_depends_on (const NonlinearDependence::Dependence dependence) const;

        /**
	 * Return true if the compressibility() function returns something that
	 * may depend on the variable identifies by the argument.
	 *
	 * This function must return false for all possible arguments if the
	 * is_compressible() function returns false.
	 */
	virtual bool
	compressibility_depends_on (const NonlinearDependence::Dependence dependence) const;

        /**
	 * Return true if the specific_heat() function returns something that
	 * may depend on the variable identifies by the argument.
	 */
	virtual bool
	specific_heat_depends_on (const NonlinearDependence::Dependence dependence) const;

        /**
	 * Return true if the thermal_conductivity() function returns something that
	 * may depend on the variable identifies by the argument.
	 */
	virtual bool
	thermal_conductivity_depends_on (const NonlinearDependence::Dependence dependence) const;

        /**
         * Return whether the model is compressible or not.  Incompressibility
         * does not necessarily imply that the density is constant; rather, it
         * may still depend on temperature or pressure. In the current
         * context, compressibility means whether we should solve the contuity
         * equation as $\nabla \cdot (\rho \mathbf u)=0$ (compressible Stokes)
         * or as $\nabla \cdot \mathbf{u}=0$ (incompressible Stokes).
         */
        virtual bool is_compressible () const;
        /**
         * @}
         */
    };
  }
}

#endif
