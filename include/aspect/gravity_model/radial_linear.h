// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2014-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_gravity_model_radial_linear_h
#define _aspect_gravity_model_radial_linear_h

#include <aspect/simulator_access.h>
#include <aspect/gravity_model/interface.h>

namespace aspect
{
  namespace GravityModel
  {
    /**
     * A class that describes gravity as a radial vector of linearly
     * changing magnitude with depth.
     *
     * @ingroup GravityModels
     */
    template <int dim>
    class RadialLinear : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Return the gravity vector as a function of position.
         */
        Tensor<1,dim> gravity_vector (const Point<dim> &position) const override;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        /**
         * Magnitude of the gravity vector at the surface, m/s^2
         *
         * This variable is read from the parameter file through a parameter called 'Magnitude at surface'.
         */
        double magnitude_at_surface;

        /**
         * Magnitude of the gravity vector at the bottom, m/s^2.
         * 'Bottom' means at the maximum depth of the provided geometry, for
         * a full sphere this means the center.
         *
         * This variable is read from the parameter file through a parameter called 'Magnitude at bottom'.
         */
        double magnitude_at_bottom;

    };
  }
}

#endif
