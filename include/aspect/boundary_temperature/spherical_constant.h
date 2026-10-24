// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_boundary_temperature_spherical_constant_h
#define _aspect_boundary_temperature_spherical_constant_h

#include <aspect/boundary_temperature/interface.h>

#include <aspect/simulator_access.h>

namespace aspect
{
  namespace BoundaryTemperature
  {
    /**
     * A class that implements a temperature boundary condition for a
     * spherical shell geometry in which the temperature at the inner and
     * outer surfaces (i.e. at the core-mantle and the
     * mantle-lithosphere/atmosphere boundaries) are constant.
     *
     * @ingroup BoundaryTemperatures
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
         * @copydoc aspect::BoundaryTemperature::Interface::boundary_temperature()
         */
        double boundary_temperature (const types::boundary_id boundary_indicator,
                                     const Point<dim> &position) const override;

        /**
         * Return the minimal the temperature on that part of the boundary on
         * which Dirichlet conditions are posed.
         *
         * This value is used in computing dimensionless numbers such as the
         * Nusselt number indicating heat flux.
         */
        double minimal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids) const override;

        /**
         * Return the maximal the temperature on that part of the boundary on
         * which Dirichlet conditions are posed.
         *
         * This value is used in computing dimensionless numbers such as the
         * Nusselt number indicating heat flux.
         */
        double maximal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids) const override;

        /**
         * Declare the parameters this class takes through input files. This
         * class declares the inner and outer boundary temperatures.
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
         * Temperatures at the inner and outer boundaries.
         *
         * This variable is read from the parameter file through a parameter called 'Inner temperature'.
         */
        double inner_temperature;

        /**
         * Temperatures at the inner and outer boundaries.
         *
         * This variable is read from the parameter file through a parameter called 'Outer temperature'.
         */
        double outer_temperature;

        /**
         * The boundary indicator for the inner and outer boundaries.
         */
        types::boundary_id inner_boundary_indicator;
        types::boundary_id outer_boundary_indicator;
    };
  }
}


#endif
