// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2013-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_boundary_composition_spherical_constant_h
#define _aspect_boundary_composition_spherical_constant_h

#include <aspect/boundary_composition/interface.h>

#include <aspect/simulator_access.h>


namespace aspect
{
  namespace BoundaryComposition
  {
    /**
     * A class that implements a composition boundary condition for a
     * spherical shell geometry in which the composition at the inner and
     * outer surfaces (i.e. at the core-mantle and the
     * mantle-lithosphere/atmosphere boundaries) are constant.
     * This class works for the sphere, spherical shell,
     * chunk and ellipsoidal chunk geometries.
     *
     * @ingroup BoundaryCompositions
     */
    template <int dim>
    class SphericalConstant : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Initialize some variables.
         */
        void
        initialize() override;

        /**
         * This function returns the constant compositions read from the
         * parameter file for the inner and outer boundaries.
         *
         * @copydoc aspect::BoundaryComposition::Interface::boundary_composition()
         */
        double boundary_composition (const types::boundary_id boundary_indicator,
                                     const Point<dim> &position,
                                     const unsigned int compositional_field) const override;

        /**
         * Declare the parameters this class takes through input files. This
         * class declares the inner and outer boundary compositions.
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
         * Compositions at the inner and outer boundaries.
         *
         * This variable is read from the parameter file through a parameter called 'Inner composition'.
         */
        std::vector<double> inner_composition;

        /**
         * Compositions at the inner and outer boundaries.
         *
         * This variable is read from the parameter file through a parameter called 'Outer composition'.
         */
        std::vector<double> outer_composition;

        /**
         * Boundary indicators for the inner and outer boundaries.
         */
        types::boundary_id inner_boundary_indicator;
        types::boundary_id outer_boundary_indicator;
    };
  }
}


#endif
