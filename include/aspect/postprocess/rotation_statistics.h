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


#ifndef _aspect_postprocess_rotation_statistics_h
#define _aspect_postprocess_rotation_statistics_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that computes some statistics about the rotational
     * velocity of the model (i.e. integrated net rotation and angular momentum).
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class RotationStatistics : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Evaluate the solution for some velocity statistics.
         */
        std::pair<std::string,std::string>
        execute (TableHandler &statistics) override;

        /**
         * Declare the parameters this class takes through input files.
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
         * Whether to use a constant density of one for the computation of the
         * angular momentum and moment of inertia. This is an approximation
         * that assumes that the 'volumetric' rotation is equal to the 'mass'
         * rotation. If this parameter is true this postprocessor computes
         * 'net rotation' instead of 'angular momentum'.
         */
        bool use_constant_density;

        /**
         * Whether to write the full moment of inertia tensor into the
         * statistics output instead of its norm for the current rotation
         * axis. This is a second-order symmetric tensor with
         * 6 components in 3D. In 2D this option has no effect, because
         * the rotation axis is fixed and thus it is always a scalar.
         */
        bool output_full_tensor;
    };
  }
}


#endif
