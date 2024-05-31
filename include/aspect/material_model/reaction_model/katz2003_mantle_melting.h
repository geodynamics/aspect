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

#ifndef _aspect_material_reaction_melt_katz2003_mantle_melting_h
#define _aspect_material_reaction_melt_katz2003_mantle_melting_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/postprocess/melt_statistics.h>
#include <aspect/melt.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    namespace ReactionModel
    {

      /**
      * A melt model that calculates melt fraction and entropy change
      * according to the melting model for dry peridotite of Katz, 2003.
      * This also includes a computation of the latent heat of melting (if the latent heat
      * heating model is active).
      *
      * These functions can be used in the calculation of melting and melt transport
      * in the melt_simple material model and can be extended to other material models
      *
      * @ingroup ReactionModel
      */
      template <int dim>
      class Katz2003MantleMelting : public ::aspect::SimulatorAccess<dim>
      {
        public:
          // constructor
          Katz2003MantleMelting();

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

          /**
          * Compute the change in entropy due to melting for a given @p temperature
          * and @p pressure, and under the assumption that a fraction
          * @p maximum_melt_fraction of the material has already been molten
          * previously. The entropy change is computed with respect to temperature
          * or pressure, depending on @p dependence.
          * This is needed to calculate the latent heat of melt.
          */
          double
          entropy_change (const double temperature,
                          const double pressure,
                          const double maximum_melt_fraction,
                          const NonlinearDependence::Dependence dependence) const;
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

          // entropy change upon melting
          double peridotite_melting_entropy_change;
      };
    }

  }
}

#endif
