/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_rheology_diffusion_dislocation_h
#define _aspect_material_model_rheology_diffusion_dislocation_h

#include <aspect/global.h>
#include <aspect/material_model/interface.h>
#include <aspect/material_model/utilities.h>
#include <aspect/simulator_access.h>
#include <aspect/material_model/rheology/diffusion_creep.h>
#include <aspect/material_model/rheology/dislocation_creep.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace Rheology
    {
      /**
       * A viscous rheology including diffusion and dislocation creep.
       *
       * The effective viscosity is defined as the harmonic mean of the two
       * effective viscosity functions describing diffusion and dislocation creep:
       * @f[
       *   v_\text{eff} = \left( \frac{1}{v_\text{eff}^\text{diff}} + \frac{1}{v_\text{eff}^\text{dis}} \right)^{-1}
       * @f]
       * where
       * @f[
       *   v_\text{eff}^\text{diff} = A_\text{diff}^{-1} \exp\left(\frac{E_\text{diff} + PV_\text{diff}}{RT}\right),
       * @f]
       * and
       * @f[
       *   v_\text{eff}^\text{dis} =  A_\text{dis}^{\frac{-1}{n_\text{dis}}} \dot{\varepsilon}^{\frac{1-n}{n}}
       *                        \exp\left(\frac{E_\text{diff} + PV_\text{diff}}{n_\text{dis}RT}\right)
       * @f]
       * where $\dot{\varepsilon}$ is the second invariant of the strain rate tensor,
       * $A_i$ are prefactors where $i$ corresponds to diffusion or dislocation creep,
       * $E_i$ are the activation energies, $V_i$ are the activation volumes,
       * $\rho_m$ is the mantle density, $R$ is the gas constant,
       * $T$ is temperature, and $P$ is pressure.
       *
       * Rheological parameters can be defined per-compositional field.
       * For each material parameter the user supplies a comma delimited list of
       * length N+1, where N is the number of compositional fields.  The
       * additional field corresponds to the value for background mantle.  They
       * should be ordered ``background, composition1, composition2...''
       *
       * The individual output viscosities for each compositional field are
       * averaged. The user can choose from a range of options for this
       * viscosity averaging. If only one value is given for any of these parameters,
       * all compositions are assigned the same value.
       * The first value in the list is the value assigned to "background mantle"
       * (regions where the sum of the compositional fields is < 1.0).
       */
      template <int dim>
      class DiffusionDislocation : public ::aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * Constructor.
           */
          DiffusionDislocation();

          /**
           * Calculate the viscosities for each compositional field
           * assuming that each composition experiences the same strain rate.
           */
          std::vector<double>
          calculate_isostrain_viscosities ( const double pressure,
                                            const double temperature,
                                            const SymmetricTensor<2,dim> &strain_rate) const;

          /**
           * Compute the viscosity based on the composite viscous creep law,
           * averaging over all compositional fields according to their
           * volume fractions.
           */
          double
          compute_viscosity (const double pressure,
                             const double temperature,
                             const std::vector<double> &volume_fractions,
                             const SymmetricTensor<2,dim> &strain_rate) const;

          static
          void
          declare_parameters (ParameterHandler &prm);

          void
          parse_parameters (ParameterHandler &prm);

        private:
          /**
           * Objects for computing viscous creep viscosities.
           */
          Rheology::DiffusionCreep<dim> diffusion_creep;
          Rheology::DislocationCreep<dim> dislocation_creep;

          /**
           * Defining a minimum strain rate stabilizes the viscosity calculation,
           * which involves a division by the strain rate. Units: 1/s.
           */
          double min_strain_rate;
          double minimum_viscosity;
          double maximum_viscosity;

          double log_strain_rate_residual_threshold;
          unsigned int stress_max_iteration_number;

          double grain_size;
          unsigned int n_chemical_composition_fields;

          MaterialUtilities::CompositionalAveragingOperation viscosity_averaging;

      };

    }
  }
}
#endif
