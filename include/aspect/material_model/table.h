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


#ifndef __aspect__model_table_h
#define __aspect__model_table_h

#include <aspect/material_model/interface.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * A material model that reads the essential values of coefficients from
     * tables in input files that describe their dependence as a function of
     * pressure and temperature.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class Table: public MaterialModel::InterfaceCompatibility<dim>
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

        virtual double viscosity_ratio (const double temperature,
                                        const double pressure,
                                        const std::vector<double>    &compositional_fields,
                                        const SymmetricTensor<2,dim> &strain_rate,
                                        const Point<dim> &position) const;

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
         * may still depend on temperature or pressure but is only interpreted
         * as such in the rhs of the Stokes equation. This is consistent with
         * the so called Bousinesq formulation. In the current context,
         * compressibility means whether we should solve the contuity equation
         * as $\nabla \cdot (\rho \mathbf u)=0$ (compressible Stokes) or as
         * $\nabla \cdot \mathbf{u}=0$ (incompressible Stokes).
         */
        virtual bool is_compressible () const;
        /**
         * @}
         */

        /**
         * @name Reference quantities
         * @{
         */

        /**
         * A reference viscosity
         *
         * The value here is not used in the computation of things but only in
         * postprocessing the solution when we want dimension-less quantities.
         */
        virtual double reference_viscosity () const;

        /**
         * A reference density
         *
         * The value here is not used in the computation of things but only in
         * postprocessing the solution when we want dimension-less quantities.
         */
        virtual double reference_density () const;

        /**
         * A reference thermal diffusivity $\kappa$. $\kappa$ is related to
         * the thermal conductivity $k$ as $\kappa = k/(rho c_p)$.
         *
         * The value here is not used in the computation of things but only in
         * postprocessing the solution when we want dimension-less quantities.
         */
        double reference_thermal_diffusivity () const;

        /**
         * A reference thermal expansion coefficient $\alpha$.
         *
         * The value here is not used in the computation of things but only in
         * postprocessing the solution when we want dimension-less quantities.
         */
        double reference_thermal_expansion_coefficient () const;

        /**
         * Return the thermal expansion coefficient $\alpha$ of the model,
         * possibly as a function of depth.
         */
        virtual double thermal_expansion_coefficient (const double temperature,
                                                      const double pressure,
                                                      const std::vector<double> &compositional_fields,
                                                      const Point<dim> &position) const;

        /**
         * A reference thermal specific heat $c_p$.
         *
         * The value here is not used in the computation of things but only in
         * postprocessing the solution when we want dimension-less quantities.
         */
        double reference_cp () const;
        /**
         * @}
         */

        /**
         * @name Auxiliary material properties used for postprocessing
         * @{
         */

        /**
         * the seismic pressure wave speed
         *
         * The value here is not used in the computation of things but only in
         * postprocessing the solution when we want dimension-less quantities.
         */
        virtual double seismic_Vp (const double temperature,
                                   const double pressure,
                                   const std::vector<double> &compositional_fields) const;

        /**
         * the seismic shear wave speed
         *
         * The value here is not used in the computation of things but only in
         * postprocessing the solution when we want dimension-less quantities.
         */
        virtual double seismic_Vs (const double temperature,
                                   const double pressure,
                                   const std::vector<double> &compositional_fields) const;

        /**
         * the phase of the composition at given pressure and temperature this
         * returns an integer value that is associated with a specific phase
         * for instance Majorite of PostPerovskite
         *
         * The value here is not used in the computation of things but only in
         * postprocessing the solution when we want dimension-less quantities.
         */
        virtual unsigned int thermodynamic_phase (const double temperature,
                                                  const double pressure,
                                                  const std::vector<double> &compositional_fields) const;


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
        double reference_kappa;
        double reference_specific_heat;
        double reference_alpha;
        std::string composition;
        std::string data_directory;
        bool compute_phases;
        bool model_is_compressible;

        std::string viscosity_model;
        double reference_eta;
        double exponential_T;
        double exponential_P;
        double increase_lower_mantle;
        double activation_energy_diffusion;
        double activation_volume_diffusion;
        double prefactor_diffusion;
        double activation_energy_dislocation;
        double activation_volume_dislocation;
        double prefactor_dislocation;
        double stress_exponent;

        double k_value;
    };
  }
}

#endif
