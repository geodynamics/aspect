/*
  Copyright (C) 2014 - 2022 by the authors of the ASPECT code.

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



#ifndef _aspect_postprocess_boundary_strain_rate_residual_statistics_h
#define _aspect_postprocess_boundary_strain_rate_residual_statistics_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>

namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that computes some statistics about the residual of the second invariant of the strain
     * rate at the top surface. The residual is calculated as the difference between the second invariant of
     * the modeled strain rate and a data file containing the second invariant of strain rate observations.
     *
     * @ingroup Postprocessing
     */

    template <int dim>
    class BoundaryStrainRateResidualStatistics : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * This function reads the specified input surface strain rate data files as an ascii data file.
         */
        void initialize () override;

        /**
         * This function returns the reference surface strain rate read from the data file at the given point @p p.
         */
        double
        get_data_surface_strain_rate (const Point<dim> &p) const;

        /**
         * Evaluate the solution to compute statistics about the residual of the second invariant of strain rate
         * residual at the top boundary.
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
         * Pointer to the structured data containing the surface strain rate.
         */
        std::unique_ptr<Utilities::StructuredDataLookup<dim>> strain_rate_data_lookup;

        /**
         * Directory in which the input data files are present.
         */
        std::string data_directory;

        /**
         * Filename of the input ascii data file containing the surface strain rate components.
         */
        std::string data_file_name;

        /**
         * Scale the input data by a scalar factor. Can be used to transform
         * the unit of the data (if they are not specified in SI units (/s or
         * /yr depending on the "Use years in output instead of seconds"
         * parameter).
         */
        double scale_factor;
    };
  }
}


#endif
