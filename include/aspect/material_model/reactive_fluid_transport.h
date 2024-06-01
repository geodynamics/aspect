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
#include <aspect/material_model/reaction_model/katz2003_mantle_melting.h>



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
    class ReactiveFluidTransport : public MaterialModel::MeltInterface<dim>, public MaterialModel::MeltFractionModel<dim>, public ::aspect::SimulatorAccess<dim>
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
         * Compute the maximum allowed bound water content at the input
         * pressure and temperature conditions. This is used to determine
         * how free water interacts with the solid phase.
         * @param in Object that contains the current conditions.
         * @param q unsigned int from 0-3 indexing which rock phase the equilbrium
         * bound water content is being calculated for
         */
        std::vector<double> tian_equilibrium_bound_water_content(const MaterialModel::MaterialModelInputs<dim> &in,
                                                                 unsigned int q) const;

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

        // The maximum water content for each of the 4 rock types in the tian approximation
        // method. These are important for keeping the polynomial bounded within reasonable
        // values.
        double tian_max_peridotite_water;
        double tian_max_gabbro_water;
        double tian_max_MORB_water;
        double tian_max_sediment_water;

        // The following coefficients are taken from a publication from Tian et al., 2019, and can be found
        // in Table 3 (Gabbro), Table B1 (MORB), Table B2 (Sediments) and Table B3 (peridotite).
        // LR refers to the effective enthalpy change for devolatilization reactions,
        // csat is the saturated mass fraction of water in the solid, and Td is the
        // onset temperature of devolatilization for water.
        std::vector<double> LR_peridotite_poly_coeffs {-19.0609, 168.983, -630.032, 1281.84, -1543.14, 1111.88, -459.142, 95.4143, 1.97246};
        std::vector<double> csat_peridotite_poly_coeffs {0.00115628, 2.42179};
        std::vector<double> Td_peridotite_poly_coeffs {-15.4627, 94.9716, 636.603};

        std::vector<double> LR_gabbro_poly_coeffs {-1.81745, 7.67198, -10.8507, 5.09329, 8.14519};
        std::vector<double> csat_gabbro_poly_coeffs {-0.0176673, 0.0893044, 1.52732};
        std::vector<double> Td_gabbro_poly_coeffs {-1.72277, 20.5898, 637.517};

        std::vector<double> LR_MORB_poly_coeffs {-1.78177, 7.50871, -10.4840, 5.19725, 7.96365};
        std::vector<double> csat_MORB_poly_coeffs {0.0102725, -0.115390, 0.324452, 1.41588};
        std::vector<double> Td_MORB_poly_coeffs {-3.81280, 22.7809, 638.049};

        std::vector<double> LR_sediment_poly_coeffs {-2.03283, 10.8186, -21.2119, 18.3351, -6.48711, 8.32459};
        std::vector<double> csat_sediment_poly_coeffs {-0.150662, 0.301807, 1.01867};
        std::vector<double> Td_sediment_poly_coeffs {2.83277, -24.7593, 85.9090, 524.898};

        // The polynomials breakdown above certain pressures, 10 GPa for peridotite, 26 GPa for gabbro, 16 GPa for MORB,
        // and 50 GPa for sediment. These cutoff pressures were determined by extending the pressure range in Tian et al. (2019)
        // and observing where the maximum allowed water contents jump towards infinite values.
        const std::vector<double> pressure_cutoffs {10, 26, 16, 50};

        std::vector<std::vector<double>> devolatilization_enthalpy_changes {LR_peridotite_poly_coeffs, LR_gabbro_poly_coeffs, \
                                                                             LR_MORB_poly_coeffs, LR_sediment_poly_coeffs
                                                                            };

        std::vector<std::vector<double>> water_mass_fractions {csat_peridotite_poly_coeffs, csat_gabbro_poly_coeffs, \
                                                                csat_MORB_poly_coeffs, csat_sediment_poly_coeffs
                                                               };

        std::vector<std::vector<double>> devolatilization_onset_temperatures {Td_peridotite_poly_coeffs, Td_gabbro_poly_coeffs, \
                                                                               Td_MORB_poly_coeffs, Td_sediment_poly_coeffs
                                                                              };

        /*
        * Object for computing Katz 2003 melt parameters
        */
        ReactionModel::Katz2003MantleMelting<dim> katz2003_model;

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
         * base model to be incompressible, as otherwise the advection
         * equation would only be valid for mass and not volume
         * fractions.
         *
         * The tian approximation model implements parametrized phase
         * diagrams from Tian et al., 2019 G3, https://doi.org/10.1029/2019GC008488
         * and calculates the fluid-solid reactions for four different rock types:
         * sediments, MORB, gabbro and peridotite. This is achieved by calculating the
         * maximum allowed bound water content for each composition at the current
         * Pressure-Temperature conditions, and releasing bound water as free water if:
         * (maximum bound water content < current bound water content)
         * or incorporating free water (if present) into the solid phase as bound water:
         * maximum bound water content > current bound water content
         * This model requires that 4 compositional fields named after the 4 different rock
         * types exist in the input file.
         *
         * The Katz2003 model implements anhydrous the mantle melting model from
         * Katz et. al., 2003 G3, doi:10.1029/2002GC000433.
         */
        enum ReactionScheme
        {
          no_reaction,
          zero_solubility,
          tian_approximation,
          katz2003
        }
        fluid_solid_reaction_scheme;
    };
  }
}

#endif
