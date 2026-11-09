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


#ifndef _aspect_postprocess_basic_statistics_h
#define _aspect_postprocess_basic_statistics_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that outputs some simplified statistics
     * like the Rayleigh number and other quantities that only
     * make sense in certain model setups. The output is written
     * after completing initial adaptive refinement steps.
     * The postprocessor assumes a point at the surface at the adiabatic
     * surface temperature and pressure is a reasonable reference condition
     * for computing these properties. Furthermore, the Rayleigh
     * number is computed using the model depth (i.e. not the
     * radius of the Earth), as we need a definition that is geometry
     * independent. Take care when comparing these values to published
     * studies and make sure they use the same definitions.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class BasicStatistics : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Evaluate the solution for some general statistics.
         */
        std::pair<std::string,std::string>
        execute (TableHandler &statistics) override;
    };
  }
}


#endif
