/*
  Copyright (C) 2014 - 2021 by the authors of the ASPECT code.

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



#ifndef _aspect_postprocess_velocity_points_residual_h
#define _aspect_postprocess_velocity_points_residual_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>

namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that computes some point value statistics about the velocity
     * residual in the model. The velocity residual is the rms residual between the
     * modeled and the observed velocities at individual data points specified in the input
     * ascii file.
     * The postprocessor outputs the velocity residual at individual points into text files
     * for each non linear iteration, and the output format is similar to the input data
     * format such that these generated files can be used as an input in the next iteration.
     * The postprocessor also outputs the rms residual statictics.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class VelocityPointsResidual : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor
         */
        VelocityPointsResidual ();

        /**
         * This function reads the specified input velocity ascii data files.
         */
        void initialize () override;

        /**
         * This function returns the input data point and the velocity at that point
         * from the ascii data file. The ascii data is structured in 1 dimension such
         * that the first column represents the point ids.
         * This function is called from execute() function.
         */
        std::pair <Point<dim>, Tensor<1,dim>>
        get_observed_data (const unsigned int p) const;

        /**
         * Evaluate the solution for some velocity residual statistics.
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
         * Determines if the input ascii data file has velocity components in
         * spherical, i.e., (r, phi, theta) or in cartesian, i.e., (x, y, z)
         * coordinate system.
         */
        bool use_spherical_unit_vectors;

        /**
         * Pointer to the structured data that stores the evaluation points.
         */
        std::unique_ptr<Utilities::StructuredDataLookup<1>> data_lookup;

        /**
         * Directory in which the input data files are present.
         */
        std::string data_directory;

        /**
         * Filename of the input ascii data file.
         */
        std::string data_file_name;

        /**
          * Consecutively counted number indicating the how-manyth time we will
          * create output the next time we get to it.
          */
        unsigned int output_file_number;
    };
  }
}


#endif
