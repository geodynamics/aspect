/*
  Copyright (C) 2014 - 2020 by the authors of the ASPECT code.

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



#ifndef _aspect_postprocess_boundary_velocity_residual_statistics_h
#define _aspect_postprocess_boundary_velocity_residual_statistics_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/boundary_velocity/gplates.h>


namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that computes some statistics about the velocity residual at the top surface.
     * The velocity residual is calculated as the difference between the modeled velocities and
     * input data velocities (GPlates model or ascii data).
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class BoundaryVelocityResidualStatistics : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * This function reads the specified input velocity data files, i.e., either an ascii data file or
         * a file from the GPlates model.
         */
        void initialize () override;

        /**
         * This function returns the input data velocity, (GPlates model or ascii data)
         * value at a point.  This function is called from execute() function.
         */
        Tensor<1,dim>
        get_data_velocity (const Point<dim> &p) const;

        /**
         * Evaluate the solution statistics for some velocity residual at the top boundary.
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
        * Determines if the model velocities are compared against ascii data
        * files, or a gplates model.
        */
        bool use_ascii_data;

        /**
        * Pointer to the gplates boundary velocity model
        */
        std::unique_ptr<BoundaryVelocity::internal::GPlatesLookup<dim>> gplates_lookup;

        /**
         * Pointer to the ascii data
         */
        std::unique_ptr<Utilities::AsciiDataLookup<dim>> ascii_data_lookup;

        /**
         * Directory in which the input data files, i.e., GPlates model or ascii data
         * are present.
         */
        std::string data_directory;

        /**
         * Filename of the input Gplates model or ascii data file. For GPlates, the file names
         * can contain the specifiers %s and/or %c (in this order), meaning the name of the
         * boundary and the number of the data file time step.
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
