/*
  Copyright (C) 2024 by the authors of the ASPECT code.

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

#ifndef _aspect_material_reaction_model_melt_tian2019_solubility_h
#define _aspect_material_reaction_model_melt_tian2019_solubility_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/melt.h>
#include <aspect/utilities.h>


namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    namespace ReactionModel
    {

      /**
       * A melt model that calculates the solubility of water according to
       * parameterized phase diagrams for four lithologies:
       *  1) sediment
       *  2) mid-ocean ridge basalt (MORB)
       *  3) gabbro
       *  4) peridotite
       * from Tian, 2019 https://doi.org/10.1029/2019GC008488.
       *
       * These functions can be used in the calculation of reactive fluid transport
       * of water.
       *
       * @ingroup ReactionModel
       */
      template <int dim>
      class Tian2019Solubility : public ::aspect::SimulatorAccess<dim>
      {
        public:

          /**
           * Compute the free fluid fraction that is present in the material based on the
           * fluid content of the material and the fluid solubility for the given input conditions.
           * @p in and @p melt_fraction need to have the same size.
           *
           * @param in Object that contains the current conditions.
           * @param porosity_idx the index of the "porosity" composition
           * @param q the quadrature point index
           */
          double
          melt_fraction (const MaterialModel::MaterialModelInputs<dim> &in,
                         const unsigned int porosity_idx,
                         unsigned int q) const;

          /**
           * Compute the maximum allowed bound water content at the input
           * pressure and temperature conditions. This is used to determine
           * how free water interacts with the solid phase.
           * @param in Object that contains the current conditions.
           * @param q the quadrature point index
           */
          std::vector<double> tian_equilibrium_bound_water_content(const MaterialModel::MaterialModelInputs<dim> &in,
                                                                   unsigned int q) const;

          /**
           * Declare the parameters this function takes through input files.
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Read the parameters from the parameter file.
           */
          void
          parse_parameters (ParameterHandler &prm);

        private:

          /**
           * The maximum water content for each of the 4 rock types in the tian approximation
           * method. These are important for keeping the polynomial bounded within reasonable
           * values.
           */
          double tian_max_peridotite_water;
          double tian_max_gabbro_water;
          double tian_max_MORB_water;
          double tian_max_sediment_water;

          /**
           *
           * The following coefficients are taken from a publication from Tian et al., 2019, and can be found
           * in Table 3 (Gabbro), Table B1 (MORB), Table B2 (Sediments) and Table B3 (peridotite).
           * LR refers to the effective enthalpy change for devolatilization reactions,
           * csat is the saturated mass fraction of water in the solid, and Td is the
           * onset temperature of devolatilization for water.
           */
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

          /**
           * The polynomials breakdown above certain pressures, 10 GPa for peridotite, 26 GPa for gabbro, 16 GPa for MORB,
           * and 50 GPa for sediment. These cutoff pressures were determined by extending the pressure range in Tian et al. (2019)
           * and observing where the maximum allowed water contents jump towards infinite values.
           */
          const std::array<double,4 > pressure_cutoffs {{10, 26, 16, 50}};

          std::vector<std::vector<double>> devolatilization_enthalpy_changes {LR_peridotite_poly_coeffs, LR_gabbro_poly_coeffs, \
                                                                               LR_MORB_poly_coeffs, LR_sediment_poly_coeffs
                                                                              };

          std::vector<std::vector<double>> water_mass_fractions {csat_peridotite_poly_coeffs, csat_gabbro_poly_coeffs, \
                                                                  csat_MORB_poly_coeffs, csat_sediment_poly_coeffs
                                                                 };

          std::vector<std::vector<double>> devolatilization_onset_temperatures {Td_peridotite_poly_coeffs, Td_gabbro_poly_coeffs, \
                                                                                 Td_MORB_poly_coeffs, Td_sediment_poly_coeffs
                                                                                };
      };
    }

  }
}

#endif
