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

#ifndef _aspect_material_model_diffusion_dislocation_h
#define _aspect_material_model_diffusion_dislocation_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/material_model/rheology/diffusion_dislocation.h>

namespace aspect
{
  namespace MaterialModel
  {
    /**
     * A material model based on a viscous rheology including diffusion and
     * dislocation creep.
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
     * Several model parameters (reference densities, thermal expansivities
     * and rheology parameters) can be defined per-compositional field.
     * For each material parameter the user supplies a comma delimited list of
     * length N+1, where N is the number of compositional fields.  The
     * additional field corresponds to the value for background mantle.  They
     * should be ordered ``background, composition1, composition2...''
     *
     * If a list of values is given for the density and thermal expansivity,
     * the volume weighted sum of the values of each of the compositional fields
     * is used in their place, for example
     * $\rho = \sum \left( \rho_i V_i \right)$
     *
     * The individual output viscosities for each compositional field are
     * also averaged. The user can choose from a range of options for this
     * viscosity averaging. If only one value is given for any of these parameters,
     * all compositions are assigned the same value.
     * The first value in the list is the value assigned to "background mantle"
     * (regions where the sum of the compositional fields is < 1.0).
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class DiffusionDislocation : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:

        void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                      MaterialModel::MaterialModelOutputs<dim> &out) const override;

        /**
         * Return whether the model is compressible or not.  Incompressibility
         * does not necessarily imply that the density is constant; rather, it
         * may still depend on temperature or pressure. In the current
         * context, compressibility means whether we should solve the
         * continuity equation as $\nabla \cdot (\rho \mathbf u)=0$
         * (compressible Stokes) or as $\nabla \cdot \mathbf{u}=0$
         * (incompressible Stokes).
         *
         * This material model is incompressible.
         */
        bool is_compressible () const override;

        static
        void
        declare_parameters (ParameterHandler &prm);

        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        /**
         * Object for computing viscous creep viscosities.
         */
        Rheology::DiffusionDislocation<dim> diffusion_dislocation;

        double reference_T;

        double thermal_diffusivity;
        double heat_capacity;

        std::vector<double> densities;
        std::vector<double> thermal_expansivities;

        MaterialUtilities::CompositionalAveragingOperation viscosity_averaging;

    };

  }
}

#endif
