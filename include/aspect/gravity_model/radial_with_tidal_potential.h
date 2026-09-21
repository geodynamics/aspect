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


#ifndef _aspect_gravity_model_radial_with_tidal_potential_h
#define _aspect_gravity_model_radial_with_tidal_potential_h

#include <aspect/simulator_access.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/gravity_model/radial_constant.h>

namespace aspect
{
  namespace GravityModel
  {
    /**
     * A class that describes gravity as a radial vector of linearly
     * changing magnitude, which is modified by a tidal potential from flattening and non-synchronous rotation (NSR).
     *
     * The equation implemented in this gravity model is from Tobie et al. (2025) (https://doi.org/10.1007/s11214-025-01136-y),
     * which is defined as:
     * $\mathbf g = -\text{magnitude} - \nabla (-V_T)$.
     * Tidal potential $V_T$ is positive because the formula follows conventions from geodesy research, where potential is taken as positive. We use
     * $V_T = \frac{3 G M_p}{2 a_s^3} r^2 (T_\ast + T_0)$.
     * $T_\ast = \frac{1}{6} (1-3 \cos(\theta)^2)$ for static flattening of the satellite and
     * $T_0=\frac{1}{2}\sin(\theta)^2 \cos(2\lambda + 2bt)$ for the quasi-static stretching in the direction of the planet varying in nonsynchronous rotation,
     * where $G$ is the gravitational constant, $M_p$ is the mass of the perturbing body, $a_s$ is the semimajor axis of the orbit, $b$ is the angular rate of non-synchronous rotation, and
     * $b = 2  \frac{\pi}{P}$ where $P$ is the period of NSR.
     * $r$, $\theta$ and $\lambda$ are radial distance, polar angle, and azimuthal angle, respectively.
     *
     * @ingroup GravityModels
     */
    template <int dim>
    class RadialWithTidalPotential : public Interface<dim>, public SimulatorAccess<dim>
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
         * The radially constant part of the gravity is handled by an object of type RadialConstant.
         */
        RadialConstant<dim> radialconstant;

        /**
         * Mass of the perturbing body.
         * This variable is read from the parameter file through a parameter called 'Mass of perturbing body'.
         */
        double M_p;

        /**
         * Semimajor axis of the orbit that causes the tidal perturbation.
         * This variable is read from the parameter file through a parameter called 'Semimajor axis of orbit'.
         */
        double a_s;

        /**
         * Period of the non-synchronous rotation.
         * This variable is read from the parameter file through a parameter called 'Period of nonsynchronous rotation'.
         */
        double P;
    };
  }
}

#endif
