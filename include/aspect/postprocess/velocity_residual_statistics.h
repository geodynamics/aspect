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



#ifndef _aspect_postprocess_velocity_residual_statistics_h
#define _aspect_postprocess_velocity_residual_statistics_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that computes some statistics about the velocity residual in the model.
     * The velocity residual is calculated as the difference between the modeled velocities and
     * input ascii data velocities.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class VelocityResidualStatistics : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * This function reads the specified input velocity ascii data files.
         */
        void initialize () override;

        /**
         * This function returns the input data velocity at a point.
         * This function is called from execute() function.
         */
        Tensor<1,dim>
        get_data_velocity (const Point<dim> &p) const;

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
         * Pointer to the structured data
         */
        std::unique_ptr<Utilities::StructuredDataLookup<dim>> data_lookup;

        /**
         * Directory in which the input data files are present.
         */
        std::string data_directory;

        /**
         * Filename of the input ascii data file.
         */
        std::string data_file_name;

        /**
         * Scale the input data by a scalar factor. Can be used to transform
         * the unit of the data (if they are not specified in SI units (m/s or
         * m/yr depending on the "Use years in output instead of seconds"
         * parameter).
         */
        double scale_factor;
    };
  }
}


#endif
