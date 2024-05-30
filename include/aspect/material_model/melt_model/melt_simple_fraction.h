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

#ifndef _aspect_material_model_melt_model_melt_simple_fraction_h
#define _aspect_material_model_melt_model_melt_simple_fraction_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/postprocess/melt_statistics.h>
#include <aspect/melt.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    namespace MeltModel
    {
      template <int dim>
      class MeltSimpleFraction : public ::aspect::SimulatorAccess<dim>
      {
        public:
          // constructor
          MeltSimpleFraction();

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
           * Percentage of material that is molten for a given @p temperature and
           * @p pressure (assuming equilibrium conditions). Melting model after Katz,
           * 2003, for dry peridotite.
           */
          double
          melt_fraction (const double temperature,
                         const double pressure) const;
        private:
          /**
          * Parameters for anhydrous melting of peridotite after Katz, 2003
          */

          // for the solidus temperature
          double A1;   // °C
          double A2; // °C/Pa
          double A3; // °C/(Pa^2)

          // for the lherzolite liquidus temperature
          double B1;   // °C
          double B2;   // °C/Pa
          double B3; // °C/(Pa^2)

          // for the liquidus temperature
          double C1;   // °C
          double C2;  // °C/Pa
          double C3; // °C/(Pa^2)

          // for the reaction coefficient of pyroxene
          double r1;     // cpx/melt
          double r2;     // cpx/melt/GPa
          double M_cpx;  // mass fraction of pyroxene

          // melt fraction exponent
          double beta;
      };
    }

  }
}

#endif
