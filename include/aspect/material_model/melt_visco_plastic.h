/*
  Copyright (C) 2015 - 2017 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_melt_visco_plastic_h
#define _aspect_material_model_melt_visco_plastic_h

#include <aspect/material_model/interface.h>
#include <aspect/material_model/visco_plastic.h>
#include <aspect/postprocess/melt_statistics.h>
#include <aspect/melt.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * A material model that implements a simple formulation of the
     * material parameters required for the modelling of melt transport,
     * including a source term for the porosity according to the melting
     * model for dry peridotite of Katz, 2003. This also includes a
     * computation of the latent heat of melting (if the latent heat
     * heating model is active).
     *
     * Most of the material properties are constant, except for the shear,
     * compaction and melt viscosities and the permeability, which depend on
     * the porosity; and the solid and melt densities, which depend on
     * temperature and pressure.
     *
     * The model is compressible (following the definition described in
     * Interface::is_compressible) only if this is specified in the input file,
     * and contains compressibility for both solid and melt.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class MeltViscoPlastic : public MaterialModel::ViscoPlastic<dim>, public MaterialModel::MeltFractionModel<dim>
    {
      public:
        /**
         * Initialization function. This function is called once at the
         * beginning of the program after parse_parameters is run and after
         * the SimulatorAccess (if applicable) is initialized.
         */
        virtual
        void
        initialize ();

        /**
         * Return whether the model is compressible or not.  Incompressibility
         * does not necessarily imply that the density is constant; rather, it
         * may still depend on temperature or pressure. In the current
         * context, compressibility means whether we should solve the continuity
         * equation as $\nabla \cdot (\rho \mathbf u)=0$ (compressible Stokes)
         * or as $\nabla \cdot \mathbf{u}=0$ (incompressible Stokes).
         */
        virtual bool is_compressible () const;

        /**
         * @name Reference quantities
         * @{
         */
        virtual double reference_viscosity () const;

        virtual void evaluate(const typename Interface<dim>::MaterialModelInputs &in,
                              typename Interface<dim>::MaterialModelOutputs &out) const;

        virtual void melt_fractions (const MaterialModel::MaterialModelInputs<dim> &in,
                                     std::vector<double> &melt_fractions) const;
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
        double melt_density_change;
        double xi_0;
        double eta_f;
        double reference_permeability;
        double alpha_phi;
        double depletion_density_change;
        double freezing_rate;

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

        // strain weakening parameters
        bool use_strain_weakening;
        std::vector<double> start_strain_weakening_intervals;
        std::vector<double> end_strain_weakening_intervals;
        std::vector<double> cohesion_strain_weakening_factors;
        std::vector<double> friction_strain_weakening_factors;

        // plasticity parameters
        std::vector<double> angles_internal_friction;
        std::vector<double> cohesions;
        std::vector<double> strength_reductions;
        std::vector<int> field_used_in_viscosity_averaging;
        double min_strain_rate;
        double min_visc;
        double max_visc;

        // elasticity parameters
        std::vector<double> elastic_shear_moduli;

        std::vector<double> compute_volume_fractions(
          const std::vector<double> &compositional_fields) const;

        /**
         * Percentage of material that is molten for a given @p temperature and
         * @p pressure (assuming equilibrium conditions). Melting model after Katz,
         * 2003, for dry peridotite.
         */
        virtual
        double
        melt_fraction (const double temperature,
                       const double pressure) const;
    };

  }
}

#endif
