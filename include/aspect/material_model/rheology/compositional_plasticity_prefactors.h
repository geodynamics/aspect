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
       * A class that handles multiplication of viscosity for a given compositional
       * field. The multiplication factors for each composition (viscosity
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
           * Compute the viscosity.
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
           * The viscosity prefactors or terms used to calculate the viscosity
           * prefactors, which are read in from the input file by the
           * parse_parameters() function. Users can choose between different schemes.
           * none: no viscosity change
           * hk04_olivine_hydration: calculate the viscosity change due to hydrogen
           * incorporation into olivine using Hirth & Kohlstaedt 2004 10.1029/138GM06.
           * This method requires a composition called 'bound_fluid' which tracks the wt%
           * water in the solid, which is used to compute an atomic ratio of H/Si ppm
           * assuming 90 mol% forsterite and 10 mol% fayalite, and finally calculates
           * a water fugacity.
           * The prefactor for a given compositional field is multiplied with a
           * base_viscosity value provided by the material model, which is then returned
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
           * Parsed functions that specify the friction angle which must be
           * given in the input file using the function method.
           */
          std::unique_ptr<Functions::ParsedFunction<dim>> prefactor_function;

          /**
           * The coordinate representation to evaluate the function for the friction angle.
           * Possible choices are depth, cartesian and spherical.
           */
          Utilities::Coordinates::CoordinateSystem coordinate_system_prefactor_function;

      };
    }
  }
}
#endif
