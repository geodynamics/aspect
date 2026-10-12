// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2020-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#ifndef _aspect_material_model_rheology_frank_kamenetskii_h
#define _aspect_material_model_rheology_frank_kamenetskii_h

#include <aspect/global.h>
#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace Rheology
    {
      /**
       * A class that computes a Frank-Kamenetskii viscosity approximation
       * of the form:
       * viscosity = A * exp(E * 0.5 * (1.0-(T/ref_T)) + F * (P-ref_P)/(rho*g*h))
       * A: prefactor of viscosity, E: adjusted viscosity ratio,
       * ref_T: reference temperature, T: temperature. F: prefactor of pressure,
       * ref_P: reference pressure, rho: density, g: gravity, h, model depth
       *
       * Refer to Noack and Breuer, 2013, GJI. doi: 10.1093/gji/ggt248 Eq. 2.10 for reference.
       */

      template <int dim>
      class FrankKamenetskii : public ::aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * Constructor.
           */
          FrankKamenetskii();

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
           * Compute the viscosity based on the Frank-Kamenetskii approximation.
           */
          double
          compute_viscosity (const double temperature,
                             const unsigned int composition,
                             const double pressure = std::numeric_limits<double>::infinity(),
                             const double density = std::numeric_limits<double>::infinity(),
                             const double gravity = std::numeric_limits<double>::infinity()) const;

        private:
          /**
           * List of Frank-Kamenetskii viscosity ratios (E).
           *
           * This variable is read from the parameter file through a parameter called 'Viscosity ratios for Frank Kamenetskii'.
           */
          std::vector<double> viscosity_ratios_frank_kamenetskii;

          /**
           * List of Frank-Kamenetskii prefactors (A).
           * This variable is read from the parameter file through a parameter called 'Prefactors for Frank Kamenetskii'.
           */
          std::vector<double> prefactors_frank_kamenetskii;

          /**
           * List of Frank-Kamenetskii pressure prefactors (F).
           *
           * This variable is read from the parameter file through a parameter called 'Pressure prefactors for Frank Kamenetskii'.
           */
          std::vector<double> pressure_prefactors_frank_kamenetskii;

          /**
           * This variable is read from the parameter file through a parameter called 'Reference temperatures for Frank Kamenetskii'.
           */
          std::vector<double> reference_temperatures;
          /**
           * This variable is read from the parameter file through a parameter called 'Reference pressures for Frank Kamenetskii'.
           */
          std::vector<double> reference_pressures;
      };
    }
  }
}
#endif
