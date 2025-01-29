/*
  Copyright (C) 2015 - 2024 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_melt_simple_h
#define _aspect_material_model_melt_simple_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/postprocess/melt_statistics.h>
#include <aspect/melt.h>
#include <aspect/material_model/reaction_model/katz2003_mantle_melting.h>

namespace aspect
{
  namespace MaterialModel
  {
    /**
     * A material model that implements a simple formulation of the
     * material parameters required for the modeling of melt transport,
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
    class MeltSimple : public MaterialModel::MeltInterface<dim>,
      public MaterialModel::MeltFractionModel<dim>,
      public ::aspect::SimulatorAccess<dim>
    {
      public:
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
         * Initialization function. This function is called once at the
         * beginning of the program after parse_parameters is run and after
         * the SimulatorAccess (if applicable) is initialized.
         */
        void
        initialize () override;

        void evaluate(const typename Interface<dim>::MaterialModelInputs &in,
                      typename Interface<dim>::MaterialModelOutputs &out) const override;

        void melt_fractions (const MaterialModel::MaterialModelInputs<dim> &in,
                             std::vector<double> &melt_fractions) const override;

        /**
         * @name Reference quantities
         * @{
         */
        double reference_darcy_coefficient () const override;

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
        void
        parse_parameters (ParameterHandler &prm) override;

        /**
         * @}
         */

        void
        create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const override;


      private:
        bool model_is_compressible;
        double thermal_expansivity;
        double eta_0;
        double reference_specific_heat;
        double thermal_conductivity;
        double compressibility;
        double thermal_viscosity_exponent;
        double reference_T;
        double depletion_density_change;
        double reference_rho_solid;



        /*
        * Object for computing the melt parameters
        */
        ReactionModel::Katz2003MantleMelting<dim> katz2003_model;

    };

  }
}

#endif
