// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2019 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_postprocess_boundary_pressures_h
#define _aspect_postprocess_boundary_pressures_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that computes the laterally averaged pressure
     * at the top and bottom boundaries of the solution.  This is useful
     * for calculating the dynamic topography at those surfaces.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class BoundaryPressures : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Evaluate the solution for the laterally averaged pressure at
         * the top and bottom of the domain.
         */
        std::pair<std::string,std::string>
        execute (TableHandler &statistics) override;

        /**
         * Query the pressure at the top surface
         */
        double pressure_at_top() const;

        /**
         * Query the pressure at the bottom surface
         */
        double pressure_at_bottom() const;

      private:

        /**
         * Pressure at the top of the domain.
         * Filled when execute() is called.
         */
        double top_pressure;

        /**
         * Pressure at the bottom of the domain.
         * Filled when execute() is called.
         */
        double bottom_pressure;
    };
  }
}


#endif
