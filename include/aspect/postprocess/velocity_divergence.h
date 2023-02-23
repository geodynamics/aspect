/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


#ifndef _aspect_postprocess_velocity_divergence_h
#define _aspect_postprocess_velocity_divergence_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>
#include <cmath>
#include <limits>
#include <deal.II/fe/fe_values.h>
#include <aspect/geometry_model/sphere.h>

namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that computes some statistics about the velocity.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class VelocityDivergence : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:

        /**
         * Evaluate the solution for some velocity statistics.
         */
        std::pair<std::string,std::string>
        execute (TableHandler &statistics) override;      
        // /**
        //  * Evaluate the solution for some velocity statistics.
        //  */
        // std::pair<std::string,std::string>
        // execute (TableHandler &statistics) override;

        // /**
        //  * Declare the parameters this class takes through input files.
        //  */
        // static
        // void
        // declare_parameters (ParameterHandler &prm);

        // /**
        //  * Read the parameters this class declares from the parameter file.
        //  */
        // void
        // parse_parameters (ParameterHandler &prm) override;

      private:
        /**
         * Parameters to define the depth until which the divergence should be calculated and the thresholds to track subductions and rifts
        //  */
        // double max_depth;
        // double subduction_threshold;
        // double rift_threshold;

      //  /**
      //    * Extent of the whole model domain in x-, y-, and z-direction (in 3d).
      //    */
      //   Point<dim> extents;    

    };
  }
}


#endif