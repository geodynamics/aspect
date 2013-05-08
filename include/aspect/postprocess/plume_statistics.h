/*
  Copyright (C) 2011, 2012 by the authors of the ASPECT code.

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
/*  $Id: temperature_statistics.h 1143 2012-09-24 18:58:15Z heister $  */


#ifndef __aspect__postprocess_plume_statistics_h
#define __aspect__postprocess_plume_statistics_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator.h>

namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that computes some statistics about a rising mantle plume,
     * such as volume, average temperature etc.
     * Every point that has a non-adiabatic temperature higher than 100 K and is
     * in a depth smaller than 300 km is considered to belong to a plume
     * approaching the surface. For this plume, the maximal temperature,
     * average temperature and volume are calculated.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class PlumeStatistics : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Evaluate the solution for some plume statistics.
         **/
        virtual
        std::pair<std::string,std::string>
        execute (TableHandler &statistics);
    };
  }
}


#endif
