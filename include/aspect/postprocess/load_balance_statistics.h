/*
  Copyright (C) 2016 - 2019 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


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
