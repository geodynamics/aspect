/*
  Copyright (C) 2025 - by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_thermal_conductivity_tosi_stackhouse_h
#define _aspect_material_model_thermal_conductivity_tosi_stackhouse_h

#include <aspect/material_model/thermal_conductivity/interface.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace ThermalConductivity
    {
      using namespace dealii;

      /**
       * A class that implements a pressure- and temperature-dependent thermal conductivity
       * following the formulation of
       *
       * Nicola Tosi, David A. Yuen, Nico de Koker, Renata M. Wentzcovitch,
       * Mantle dynamics with pressure- and temperature-dependent thermal expansivity and conductivity,
       * Physics of the Earth and Planetary Interiors, Volume 217,
       * 2013, Pages 48-58, ISSN 0031-9201, https://doi.org/10.1016/j.pepi.2013.02.004,
       *
       * and
       *
       * Stephen Stackhouse, Lars Stixrude, Bijaya B. Karki,
       * First-principles calculations of the lattice thermal conductivity of the lower mantle,
       * Earth and Planetary Science Letters, Volume 427, 2015, Pages 11-17, ISSN 0012-821X,
       * https://doi.org/10.1016/j.epsl.2015.06.050.
       *
       * The thermal conductivity parameter sets can be chosen in such a
       * way that either the Stackhouse or the Tosi relations are used.
       * The conductivity description can consist of several layers with
       * different sets of parameters. Note that the Stackhouse
       * parametrization is only valid for the lower mantle (bridgmanite).
       *
       * The default parameters of this class use the Tosi parametrization in the upper
       * mantle and the Stackhouse parametrization in the lower mantle, which is how
       * it was used in the publication
       *
       * Juliane Dannberg, Rene Gassmöller, Daniele Thallner, Frederick LaCombe, Courtney Sprain,
       * Changes in core–mantle boundary heat flux patterns throughout the supercontinent cycle,
       * Geophysical Journal International, Volume 237, Issue 3, June 2024, Pages 1251–1274, https://doi.org/10.1093/gji/ggae075,
       *
       * which introduced this implementation.
       *
       * @ingroup MaterialModels
       */
      template <int dim>
      class TosiStackhouse : public Interface<dim>, public aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * Function to compute the thermal conductivities in @p out given the
           * inputs in @p in.
           */
          void evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                         MaterialModel::MaterialModelOutputs<dim> &out) const override;

          /**
           * Declare the parameters this plugin takes through input files.
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Read the parameters from the parameter file.
           */
          void
          parse_parameters (ParameterHandler &prm) override;

        private:
          /**
           * Parameters for the temperature and pressure dependence of the
           * thermal conductivity.
           */
          std::vector<double> conductivity_transition_depths;
          std::vector<double> reference_thermal_conductivities;
          std::vector<double> conductivity_pressure_dependencies;
          std::vector<double> conductivity_reference_temperatures;
          std::vector<double> conductivity_exponents;
          std::vector<double> saturation_scaling;

          /**
           * The maximum allowed thermal conductivity.
           */
          double maximum_conductivity;
      };
    }
  }
}

#endif
