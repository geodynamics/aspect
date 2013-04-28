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


#ifndef __aspect__model_simple_h
#define __aspect__model_simple_h

#include <aspect/material_model/interface.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * A namespace for the implementation of material models that correspond
     * to the benchmarks defined in the following paper:
     * @code
     *  @Article{DMGT11,
     *    author =       {T. Duretz and D. A. May and T. V. Gerya and P. J. Tackley},
     *    title =        {Discretization errors and free surface stabilization in the
     *                  finite difference and marker-in-cell method for applied
     *                  geodynamics: {A} numerical study},
     *    journal =      {Geochemistry Geophysics Geosystems},
     *    year =         2011,
     *    volume =       12,
     *    pages =        {Q07004/1--26}}
     * @endcode
     *
     * @note While this paper summarizes the benchmarks used here, some
     * of the benchmarks actually originate in earlier papers. For the original
     * references, see the bibliography of the paper above.
     */
    namespace DuretzEtAl
    {
      /**
       * A material model that describes the <i>SolCx</i> benchmark of the paper
       * cited in the documentation of the DuretzEtAl namespace.
       *
       * @note The SolCx benchmark only talks about the flow field, not about
       * a temperature field. All quantities related to the temperature are
       * therefore set to zero in the implementation of this class.
       *
       * @note The analytic solution of this benchmark is implemented in the
       * "SolCx error" postprocessor in aspect::Postprocessor::DuretzEtAl::SolCx
       * class and can be used to assess the accuracy of the computed solution.
       *
       * @ingroup MaterialModels
       */
      template <int dim>
      class SolCx : public MaterialModel::InterfaceCompatibility<dim>
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
            * Returns the viscosity value on the right half of the domain, typically 1 or 1e6
            */
          double get_eta_B() const;

          /**
            * Returns the background density of this model. See the corresponding
            * member variable of this class for more information.
            */
          double get_background_density() const;

        private:
          /**
            * Viscosity value on the right half of the domain, typically 1 or 1e6
            */
          double eta_B;

          /**
           * A constant background density over which the density
           * variations are overlaid. This constant density has no
           * effect on the dynamic pressure and consequently on
           * the flow field, but it contributes to the total pressure
           * via the adiabatic pressure. We use this field to support
           * our claim in the first ASPECT paper that the accuracy
           * of the solutions is guaranteed even if we don't
           * subtract the adiabatic pressure in our computations.
           */
          double background_density;
      };


      /**
       * A material model that describes the <i>SolKz</i> benchmark of the paper
       * cited in the documentation of the DuretzEtAl namespace.
       *
       * @note The SolKz benchmark only talks about the flow field, not about
       * a temperature field. All quantities related to the temperature are
       * therefore set to zero in the implementation of this class.
       *
       * @ingroup MaterialModels
       */
      template <int dim>
      class SolKz : public MaterialModel::InterfaceCompatibility<dim>
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
      };

      /**
       * A material model that describes the "Pure shear/Inclusion" benchmark of the paper
       * cited in the documentation of the DuretzEtAl namespace.
       *
       * @note This benchmark only talks about the flow field, not about
       * a temperature field. All quantities related to the temperature are
       * therefore set to zero in the implementation of this class.
       *
       * @ingroup MaterialModels
       */
      template <int dim>
      class Inclusion : public MaterialModel::InterfaceCompatibility<dim>
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
            * Returns the viscosity value in the inclusion
            */
          double get_eta_B() const;

        private:
          /**
            * viscosity value in the inclusion
            */
          double eta_B;
      };
    }
  }
}

#endif
