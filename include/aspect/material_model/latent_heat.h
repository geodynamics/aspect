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


#ifndef __aspect__model_latent_heat_h
#define __aspect__model_latent_heat_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * A material model that implements a standard approximation of the latent
     * heat terms following Christensen \& Yuen, 1986. The change of entropy
     * is calculated as $Delta S = \gamma \frac{\Delta\rho}{\rho^2}$ with the
     * Clapeyron slope $\gamma$ and the density change $\Delta\rho$ of the
     * phase transition being input parameters. This model employs an analytic
     * phase function in the form $X=0.5 \left( 1 + \tanh \left( \frac{\Delta
     * p}{\Delta p_0} \right) \right)$ with $\Delta p = p - p_transition -
     * \gamma \left( T - T_transition \right)$ and $\Delta p_0$ being the
     * pressure difference over the width of the phase transition (specified
     * as input parameter).
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class LatentHeat : public MaterialModel::InterfaceCompatibility<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * @name Physical parameters used in the basic equations
         * @{
         */
        virtual double viscosity (const double                  temperature,
                                  const double                  pressure,
                                  const std::vector<double>    &compositional_fields,
                                  const SymmetricTensor<2,dim> &strain_rate,
                                  const Point<dim>             &position) const;

        virtual double density (const double temperature,
                                const double pressure,
                                const std::vector<double> &compositional_fields,
                                const Point<dim> &position) const;

        virtual double compressibility (const double temperature,
                                        const double pressure,
                                        const std::vector<double> &compositional_fields,
                                        const Point<dim> &position) const;

        virtual double specific_heat (const double temperature,
                                      const double pressure,
                                      const std::vector<double> &compositional_fields,
                                      const Point<dim> &position) const;

        virtual double thermal_expansion_coefficient (const double      temperature,
                                                      const double      pressure,
                                                      const std::vector<double> &compositional_fields,
                                                      const Point<dim> &position) const;

        virtual double thermal_conductivity (const double temperature,
                                             const double pressure,
                                             const std::vector<double> &compositional_fields,
                                             const Point<dim> &position) const;

        virtual double entropy_derivative (const double temperature,
                                           const double pressure,
                                           const std::vector<double> &compositional_fields,
                                           const Point<dim> &position,
                                           const NonlinearDependence::Dependence dependence) const;
        /**
         * @}
         */

        /**
         * @name Qualitative properties one can ask a material model
         * @{
         */

        /**
         * Return true if the viscosity() function returns something that may
         * depend on the variable identifies by the argument.
         */
        virtual bool
        viscosity_depends_on (const NonlinearDependence::Dependence dependence) const;

        /**
         * Return true if the density() function returns something that may
         * depend on the variable identifies by the argument.
         */
        virtual bool
        density_depends_on (const NonlinearDependence::Dependence dependence) const;

        /**
         * Return true if the compressibility() function returns something
         * that may depend on the variable identifies by the argument.
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
         * Return true if the thermal_conductivity() function returns
         * something that may depend on the variable identifies by the
         * argument.
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
         * @name Reference quantities
         * @{
         */
        virtual double reference_viscosity () const;

        virtual double reference_density () const;

        virtual double reference_thermal_expansion_coefficient () const;

//TODO: should we make this a virtual function as well? where is it used?
        double reference_thermal_diffusivity () const;

        double reference_cp () const;
        /**
         * @}
         */


        /**
         * @name Functions used in dealing with run-time parameters
         * @{
         */
        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);
        /**
         * @}
         */

      private:
        double reference_rho;
        double reference_T;
        double eta;
        double composition_viscosity_prefactor;
        double thermal_viscosity_exponent;
        double thermal_alpha;
        double reference_specific_heat;
        double reference_compressibility;

        /**
         * The thermal conductivity.
         */
        double k_value;

        double compositional_delta_rho;

        /**
         * Percentage of material that has already undergone the phase
         * transition to the higher-pressure material (this is done
         * individually for each transition and summed up in the end)
         */
        virtual
        double
        phase_function (const Point<dim> &position,
                        const double temperature,
                        const double pressure,
                        const int phase) const;

        /**
         * Derivative of the phase function (argument is the pressure
         * deviation).
         */
        virtual
        double
        phase_function_derivative (const Point<dim> &position,
                                   const double temperature,
                                   const double pressure,
                                   const int phase) const;

        // list of depth, width and Clapeyron slopes for the different phase
        // transitions
        std::vector<double> transition_depths;
        std::vector<double> transition_temperatures;
        std::vector<double> transition_widths;
        std::vector<double> transition_slopes;
        std::vector<double> density_jumps;
        std::vector<int> transition_phases;
        std::vector<double> phase_prefactors;
    };

  }
}

#endif
