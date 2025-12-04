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

#ifndef _aspect_material_model_rheology_compositional_plasticity_prefactors_h
#define _aspect_material_model_rheology_compositional_plasticity_prefactors_h

#include <aspect/global.h>
#include <aspect/material_model/interface.h>
#include <deal.II/base/parsed_function.h>
#include <aspect/simulator_access.h>
#include <aspect/melt.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace Rheology
    {
      /**
       * A class that handles multiplication of the plastic parameters for a given compositional
       * field. The multiplication factors for each composition (plasticity
       * prefactors) are also declared, parsed, and in some cases calculated in this class.
       */
      template <int dim>
      class CompositionalPlasticityPrefactors : public ::aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * Constructor.
           */
          CompositionalPlasticityPrefactors();

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

          /**
           * Compute the plasticity weakening factor.
           */
          std::pair<double, double>
          compute_weakening_factor (const MaterialModel::MaterialModelInputs<dim> &in,
                                    const unsigned int volume_fraction_index,
                                    const unsigned int i,
                                    const double static_cohesion,
                                    const double static_friction_angle,
                                    const Point<dim> &position) const;

        private:
          /**
           * The plasticity prefactors or terms used to calculate the plasticity
           * prefactors, which are read in from the input file by the
           * parse_parameters() function. Users can choose between different schemes.
           * none: no change in the plasticity parameters
           * porosity: calculate the change in the plasticity parameters due to the amount
           * of porosity at each quadrature point.
           * function: calculate the change in the plasticity parameters by specifying the
           * prefactors via a function.
           * The prefactor for a given compositional field is multiplied with a
           * cohesions and the friction angles provided by the material model, which is then returned
           * to the material model.
           */
          enum WeakeningMechanism
          {
            none,
            porosity,
            function
          } prefactor_mechanism;

          std::vector<double> max_porosities_for_cohesion_prefactors;
          std::vector<double> max_porosities_for_friction_angle_prefactors;

          std::vector<double> minimum_cohesions;
          std::vector<double> minimum_friction_angles;
          /**
           * Parsed functions that specify the plasticity prefactors which must be
           * given in the input file using the function method.
           */
          std::unique_ptr<Functions::ParsedFunction<dim>> prefactor_function;

          /**
           * The coordinate representation to evaluate the function for the plasticity prefactors.
           * Possible choices are depth, cartesian and spherical.
           */
          Utilities::Coordinates::CoordinateSystem coordinate_system_prefactor_function;

      };
    }
  }
}
#endif
