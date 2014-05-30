/*
  Copyright (C) 2014 by the authors of the ASPECT code.

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


#ifndef __aspect__postprocess_topography_h
#define __aspect__postprocess_topography_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator.h>

namespace aspect
{
  namespace Postprocess
  {

    /**
     * Loops over the cells vertices and determines the maximum and minimum
     * topography of the domain, relative to the initial box height for a box
     * geometry model, and the initial radius for spheres or spherical shells.
     * Intended for use with a free surface.
     * @ingroup Postprocessing
     */

    template <int dim>
    class Topography : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Calculate the maximum and minimum topography [m]
         */
        virtual
        std::pair<std::string,std::string>
        execute (TableHandler &statistics);
    };
  }
}


#endif
