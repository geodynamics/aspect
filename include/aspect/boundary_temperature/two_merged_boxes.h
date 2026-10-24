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


#ifndef _aspect_boundary_temperature_two_merged_boxes_h
#define _aspect_boundary_temperature_two_merged_boxes_h

#include <aspect/boundary_temperature/interface.h>

#include <aspect/simulator_access.h>


namespace aspect
{
  namespace BoundaryTemperature
  {
    /**
     * A class that implements a temperature boundary condition for a
     * box geometry with additional boundary indicators for the
     * lithospheric part of left and right (and front and back in 3D) boundaries.
     *
     * @ingroup BoundaryTemperatures
     */
    template <int dim>
    class TwoMergedBoxes : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * This function returns constant compositions at the boundaries.
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
         * The values of the various temperature variables on each of the
         * 2*dim+2*(dim-1) boundaries of the box.
         *
         * This variable is read from the parameter file through parameters called 'Left composition',
         * 'Right composition', 'Bottom composition', 'Top composition',
         * 'Front composition', and 'Back composition', as well as 'Left composition lithosphere',
         * 'Right composition lithosphere', 'Front composition lithosphere', and 'Back composition lithosphere'.
         */
        double temperature_values[2*dim+2*(dim-1)];
    };
  }
}


#endif
