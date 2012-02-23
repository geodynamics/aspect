//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011, 2012 by the authors of the ASPECT code
//
//-------------------------------------------------------------
#ifndef __aspect__model_simple_h
#define __aspect__model_simple_h

#include <aspect/material_model/interface.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * A material model that consists of globally constant values for all
     * material parameters except that the density decays linearly with the
     * temperature.
     *
     * The model is considered incompressible, following the definition
     * described in Interface::is_compressible. This is essentially
     * the material model used in the step-32 tutorial program.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class Simple : public MaterialModel::Interface<dim>
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

        double reference_thermal_diffusivity () const;

        double reference_thermal_alpha () const;

        double reference_cp () const;

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

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter
         * file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);
      private:
        double reference_rho;
        double reference_T;
        double eta;
        double thermal_alpha;
        double reference_specific_heat;

        /**
         * The thermal conductivity.
         */
        double k_value;
    };

  }
}

#endif
