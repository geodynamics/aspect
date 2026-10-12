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


#ifndef _aspect_postprocess_boundary_densities_h
#define _aspect_postprocess_boundary_densities_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that computes the laterally averaged density
     * at the top and bottom boundaries of the solution.  This is useful
     * for calculating the geoid at those surfaces.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class BoundaryDensities : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Evaluate the solution for laterally averaged density at the top
         * and bottom of the domain.
         */
        std::pair<std::string,std::string>
        execute (TableHandler &statistics) override;

        /**
         * Query the density at the top surface
         */
        double density_at_top() const;

        /**
         * Query the density at the bottom surface
         */
        double density_at_bottom() const;

      private:

        /**
         * Density at the top of the domain.
         * Filled when execute() is called.
         */
        double top_density;

        /**
         * Density at the bottom of the domain.
         * Filled when execute() is called.
         */
        double bottom_density;
    };
  }
}


#endif
