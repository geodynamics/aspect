/*
  Copyright (C) 2013 - 2023 by the authors of the ASPECT code.

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
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/

#ifndef _aspect_material_model_latent_heat_melt_h
#define _aspect_material_model_latent_heat_melt_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/melt.h>

namespace aspect
{
  namespace MaterialModel
  {
    /**
     * A material model that implements latent heat of melting for two
     * materials: peridotite and pyroxenite. The density and thermal
     * expansivity depend on the melt fraction.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class LatentHeatMelt : public MaterialModel::Interface<dim>,
      public MaterialModel::MeltFractionModel<dim>,
      public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * @name Physical parameters used in the basic equations
         * @{
         */
        void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                      MaterialModel::MaterialModelOutputs<dim> &out) const override;


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
         * context, compressibility means whether we should solve the continuity
         * equation as $\nabla \cdot (\rho \mathbf u)=0$ (compressible Stokes)
         * or as $\nabla \cdot \mathbf{u}=0$ (incompressible Stokes).
         */
        bool is_compressible () const override;
        /**
         * @}
         */

        void melt_fractions (const MaterialModel::MaterialModelInputs<dim> &in,
                             std::vector<double> &melt_fractions,
                             const MaterialModel::MaterialModelOutputs<dim> *out = nullptr) const override;


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
        void
        parse_parameters (ParameterHandler &prm) override;
        /**
         * @}
         */

      private:
        /**
         *  This variable is read from the parameter file through a parameter called 'Reference density'.
         */
        double reference_rho;
        /**
         *  This variable is read from the parameter file through a parameter called 'Reference temperature'.
         */
        double reference_T;
        /**
         *  This variable is read from the parameter file through a parameter called 'Viscosity'.
         */
        double eta;
        /**
         *  This variable is read from the parameter file through a parameter called 'Composition viscosity prefactor'.
         */
        double composition_viscosity_prefactor;
        /**
         *  This variable is read from the parameter file through a parameter called 'Thermal viscosity exponent'.
         */
        double thermal_viscosity_exponent;
        /**
         *  This variable is read from the parameter file through a parameter called 'Thermal expansion coefficient'.
         */
        double thermal_alpha;
        /**
         *  This variable is read from the parameter file through a parameter called 'Thermal expansion coefficient of melt'.
         */
        double melt_thermal_alpha;
        /**
         *  This variable is read from the parameter file through a parameter called 'Reference specific heat'.
         */
        double reference_specific_heat;
        /**
         *  This variable is read from the parameter file through a parameter called 'Compressibility'.
         */
        double reference_compressibility;

        /**
         * The thermal conductivity.
         *
         * This variable is read from the parameter file through a parameter called 'Thermal conductivity'.
         */
        double k_value;

        /**
         *  This variable is read from the parameter file through a parameter called 'Density differential for compositional field 1'.
         */
        double compositional_delta_rho;

        /**
         * Parameters for anhydrous melting of peridotite after Katz, 2003
         */

        // for the solidus temperature
        /**
         *  This variable is read from the parameter file through a parameter called 'A1'.
         */
        double A1;   // °C
        /**
         *  This variable is read from the parameter file through a parameter called 'A2'.
         */
        double A2; // °C/Pa
        /**
         *  This variable is read from the parameter file through a parameter called 'A3'.
         */
        double A3; // °C/(Pa^2)

        // for the lherzolite liquidus temperature
        /**
         *  This variable is read from the parameter file through a parameter called 'B1'.
         */
        double B1;   // °C
        /**
         *  This variable is read from the parameter file through a parameter called 'B2'.
         */
        double B2;   // °C/Pa
        /**
         *  This variable is read from the parameter file through a parameter called 'B3'.
         */
        double B3; // °C/(Pa^2)

        // for the liquidus temperature
        /**
         *  This variable is read from the parameter file through a parameter called 'C1'.
         */
        double C1;   // °C
        /**
         *  This variable is read from the parameter file through a parameter called 'C2'.
         */
        double C2;  // °C/Pa
        /**
         *  This variable is read from the parameter file through a parameter called 'C3'.
         */
        double C3; // °C/(Pa^2)

        // for the reaction coefficient of pyroxene
        /**
         *  This variable is read from the parameter file through a parameter called 'r1'.
         */
        double r1;     // cpx/melt
        /**
         *  This variable is read from the parameter file through a parameter called 'r2'.
         */
        double r2;     // cpx/melt/GPa
        /**
         *  This variable is read from the parameter file through a parameter called 'Mass fraction cpx'.
         */
        double M_cpx;  // mass fraction of pyroxene

        // melt fraction exponent
        /**
         *  This variable is read from the parameter file through a parameter called 'beta'.
         */
        double beta;

        // entropy change upon melting
        /**
         *  This variable is read from the parameter file through a parameter called 'Peridotite melting entropy change'.
         */
        double peridotite_melting_entropy_change;

        /**
         * Parameters for melting of pyroxenite after Sobolev et al., 2011
         */

        // for the melting temperature
        /**
         *  This variable is read from the parameter file through a parameter called 'D1'.
         */
        double D1;    // °C
        /**
         *  This variable is read from the parameter file through a parameter called 'D2'.
         */
        double D2;  // °C/Pa
        /**
         *  This variable is read from the parameter file through a parameter called 'D3'.
         */
        double D3; // °C/(Pa^2)
        // for the melt-fraction dependence of productivity
        /**
         *  This variable is read from the parameter file through a parameter called 'E1'.
         */
        double E1;
        /**
         *  This variable is read from the parameter file through a parameter called 'E2'.
         */
        double E2;

        // for the maximum melt fraction of pyroxenite
        /**
         *  This variable is read from the parameter file through a parameter called 'Maximum pyroxenite melt fraction'.
         */
        double F_px_max;

        // the relative density of molten material (compared to solid)
        /**
         *  This variable is read from the parameter file through a parameter called 'Relative density of melt'.
         */
        double relative_melt_density;

        /**
         *  This variable is read from the parameter file through a parameter called 'Pyroxenite melting entropy change'.
         */
        double pyroxenite_melting_entropy_change;

        /**
         * Percentage of material that is molten. Melting model after Katz,
         * 2003 (for peridotite) and Sobolev et al., 2011 (for pyroxenite)
         */
        virtual
        double
        melt_fraction (const double temperature,
                       const double pressure,
                       const std::vector<double> &compositional_fields,
                       const Point<dim> &position) const;

        virtual
        double
        peridotite_melt_fraction (const double temperature,
                                  const double pressure,
                                  const std::vector<double> &compositional_fields,
                                  const Point<dim> &position) const;

        virtual
        double
        pyroxenite_melt_fraction (const double temperature,
                                  const double pressure,
                                  const std::vector<double> &compositional_fields,
                                  const Point<dim> &position) const;

        double
        entropy_derivative ( const double temperature,
                             const double pressure,
                             const std::vector<double> &compositional_fields,
                             const Point<dim> &position,
                             const NonlinearDependence::Dependence dependence) const;
    };

  }
}

#endif
