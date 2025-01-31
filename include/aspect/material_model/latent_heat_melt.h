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
                             std::vector<double> &melt_fractions) const override;


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
        double reference_rho;
        double reference_T;
        double eta;
        double composition_viscosity_prefactor;
        double thermal_viscosity_exponent;
        double thermal_alpha;
        double melt_thermal_alpha;
        double reference_specific_heat;
        double reference_compressibility;

        /**
         * The thermal conductivity.
         */
        double k_value;

        double compositional_delta_rho;

        /**
         * Parameters for anhydrous melting of peridotite after Katz, 2003
         */

        // for the solidus temperature
        double A1;   // °C
        double A2; // °C/Pa
        double A3; // °C/(Pa^2)

        // for the lherzolite liquidus temperature
        double B1;   // °C
        double B2;   // °C/Pa
        double B3; // °C/(Pa^2)

        // for the liquidus temperature
        double C1;   // °C
        double C2;  // °C/Pa
        double C3; // °C/(Pa^2)

        // for the reaction coefficient of pyroxene
        double r1;     // cpx/melt
        double r2;     // cpx/melt/GPa
        double M_cpx;  // mass fraction of pyroxene

        // melt fraction exponent
        double beta;

        // entropy change upon melting
        double peridotite_melting_entropy_change;

        /**
         * Parameters for melting of pyroxenite after Sobolev et al., 2011
         */

        // for the melting temperature
        double D1;    // °C
        double D2;  // °C/Pa
        double D3; // °C/(Pa^2)
        // for the melt-fraction dependence of productivity
        double E1;
        double E2;

        // for the maximum melt fraction of pyroxenite
        double F_px_max;

        // the relative density of molten material (compared to solid)
        double relative_melt_density;

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
