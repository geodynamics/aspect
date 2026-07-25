
/*
  Copyright (C) 2019 - 2026 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_pressure_factor_h
#define _aspect_material_model_pressure_factor_h

#include <aspect/global.h>
#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/parsed_function.h>
#include <aspect/utilities.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    namespace Rheology
    {

      /**
       * Enumeration for selecting which pressure factor model to use.
       */
      enum class PressureFactorModel
      {
        none,
        function
      };

      template <int dim>
      class PressureFactor : public ::aspect::SimulatorAccess<dim>
      {
        public:
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
          parse_parameters (ParameterHandler &prm,
                            const std::unique_ptr<std::vector<unsigned int>> &expected_n_phases_per_composition = nullptr);

          /**
           * A function that computes a pressure factor based on the chosen model.
           * This pressure factor can, for example, be considered the pore fluid factor
           * in brittle regimes and applied to the pressure_for_plasticity.
           */
          double
          compute_pressure_factor(const MaterialModel::MaterialModelInputs<dim> &in,
                                  const unsigned int q,
                                  const unsigned int volume_fraction_index) const;

          /**
           * A function that scales the pressure based on the input pressure
           * and pressure factor. This can either scale the pressure as
           * P = P * lambda when computing e.g. the pressure accounting for
           * pore fluid pressure. Or P = P*lambda when computing the fluid
           * pressure in cases where it is not explicitly given in the model.
           */
          double
          compute_scaled_pressure(const MaterialModel::MaterialModelInputs<dim> &in,
                                  const double &pressure_to_scale,
                                  const unsigned int q,
                                  const unsigned int volume_fraction_index) const;


          /**
           * Calculate the pressure factor based on a user-defined function.
           */
          double
          pressure_factor_function(const Point<dim> &position, const unsigned int n_comp) const;

          /**
           * Return which pressure factor model is in use.
           */
          PressureFactorModel
          get_pressure_factor_model() const;

        private:
          /**
           * A function object representing the pressure factor for compositional fields.
           */
          std::unique_ptr<Functions::ParsedFunction<dim>> function;

          /**
           * The coordinate representation to evaluate the function. Possible
           * choices are depth, cartesian and spherical.
           */
          Utilities::Coordinates::CoordinateSystem coordinate_system;

          PressureFactorModel pressure_factor_model;
          bool return_solid_scaling;
      };
    }
  }
}
#endif
