/*
  Copyright (C) 2023 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_reactive_fluid_transport_h
#define _aspect_material_model_reactive_fluid_transport_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/melt.h>
#include <aspect/utilities.h>
#include <aspect/geometry_model/interface.h>

#include <aspect/melt.h>
#include <aspect/utilities.h>
#include <aspect/geometry_model/interface.h>



namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;
    /**
     * A material model that simulates both fluid-rock interactions
     * and the advection of fluids. It is designed to be composited with another material
     * model that computes the solid material properties.
     * @ingroup MaterialModels
     */

    template <int dim>
    class ReactiveFluidTransport : public MaterialModel::MeltInterface<dim>, public ::aspect::SimulatorAccess<dim>, public MaterialModel::MeltFractionModel<dim>
    {
      public:
        /**
         * @copydoc MaterialModel::Interface::is_compressible()
         *
         * Returns value from material model providing compressibility.
         */
        bool is_compressible () const override;

        /**
         * @name Reference quantities
         * @{
         */
        virtual double reference_darcy_coefficient () const override;

        /**
         * Compute the free fluid fraction that can be present in the material based on the
         * fluid content of the material and the fluid solubility for the given input conditions.
         * @p in and @p melt_fractions need to have the same size.
         *
         * @param in Object that contains the current conditions.
         * @param melt_fractions Vector of doubles that is filled with the
         * allowable free fluid fraction for each given input conditions.
         */
        virtual void melt_fractions (const MaterialModel::MaterialModelInputs<dim> &in,
                                     std::vector<double> &melt_fractions) const override;

        /**
         * Initialize the base model at the beginning of the model run
         * @copydoc MaterialModel::Interface::initialize()
         */
        void
        initialize () override;

        /**
         * Update the base model at the beginning of each timestep.
         */
        void update() override;

        /**
         * @copydoc MaterialModel::Interface::evaluate()
         */
        void
        evaluate (const typename Interface<dim>::MaterialModelInputs &in,
                  typename Interface<dim>::MaterialModelOutputs &out) const override;

        /**
         * @copydoc MaterialModel::Interface::declare_parameters()
         */
        static void
        declare_parameters (ParameterHandler &prm);

        /**
         * @copydoc MaterialModel::Interface::parse_parameters()
         */
        void
        parse_parameters (ParameterHandler &prm) override;

        /**
         * If this material model can produce additional named outputs
         * that are derived from NamedAdditionalOutputs, create them in here.
         */
        virtual
        void
        create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const override;

      private:

        /**
         * Pointer to the material model used as the base model
         */
        std::unique_ptr<MaterialModel::Interface<dim>> base_model;

        // Variables that describe the properties of the fluid, i.e. its density,
        // viscosity, and compressibility.
        // Properties of the solid are defined in the base model.
        double reference_rho_f;
        double eta_f;
        double fluid_compressibility;

        // Material properties governing the transport of the fluid with respect
        // to the solid, i.e., the bulk viscosity (relative to the shear viscosity),
        // the permeability, and how much the solid viscosity changes in the presence
        // of fluids.
        double shear_to_bulk_viscosity_ratio;
        double reference_permeability;
        double alpha_phi;

        // Time scale for fluid release and absorption.
        double fluid_reaction_time_scale;

        /**
         * Enumeration for selecting which type of scheme to use for
         * reactions between fluids and solids. The available
         * reaction models are described below.
        *
        * The no reaction model does not include any reactions
        * between the solid and fluid phases. As a result,
        * there is no exchange between the bound fluid and porosity
        * compositional fields. However, the values of each field
        * may vary through the model evolution through advection
        * from their initial configurations.
         *
        * The zero solubility model describes a scenario where the
         * solid cannot accommodate any fluid (i.e., zero solubility).
         * The fluid volume fraction in equilibrium with the solid
         * at any point (stored in the melt_fractions vector) is
         * equal to the sum of the bound fluid content and porosity,
         * with the latter determined by the assigned initial
         * porosity, fluid boundary conditions, and fluid
         * transport through the model. Significantly, this reaction
         * model is thus assuming that the bound water fraction is a
         * volume fraction (i.e., since porosity is always a volume
         * fraction). This latter assumption also requires the selected
         * base model is incompressible, as otherwise the advection
         * equation would only be valid for mass and not volume
         * fractions.
         */
        enum ReactionScheme
        {
          no_reaction,
          zero_solubility
        } fluid_solid_reaction_scheme;
    };
  }
}

#endif
