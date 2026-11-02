// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2016-2020 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_postprocess_load_balance_statistics_h
#define _aspect_postprocess_load_balance_statistics_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that computes statistics about
     * the distribution of cells, and if present particles
     * across subdomains.
     *
     * In particular, it computes maximal, average and
     * minimal number of cells across all ranks.
     * If there are particles it also computes the
     * maximal, average, and minimum number of particles across
     * all ranks, and maximal, average, and minimal ratio
     * between local number of particles and local number
     * of cells across all processes. All of these numbers
     * can be useful to assess the load balance between
     * different MPI ranks, as the difference between the
     * minimal and maximal load should be as small as
     * possible.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class LoadBalanceStatistics : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Compute load balance statistics.
         */
        std::pair<std::string,std::string>
        execute (TableHandler &statistics) override;
    };
  }
}


#endif
