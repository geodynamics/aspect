/*
  Copyright (C) 2018 by the authors of the ASPECT code.

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


#ifndef _aspect_postprocess_gravity_point_values_h
#define _aspect_postprocess_gravity_point_values_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that computes gravity, gravity anomalies, gravity potential and
     * gravity gradients for a set of points (e.g. satellites) in or above the model
     * surface for a user-defined range of latitudes, longitudes and radius, or a list
     * of point coordinates.

     * @ingroup Postprocessing
     */
    template <int dim>
    class GravityPointValues : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Specify the creation of output_gravity.txt.
         */
        virtual
        std::pair<std::string,std::string> execute (TableHandler &);

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);
        /**
         * @}
         */

      private:
        /**
         * Quadrature degree increase over the velocity element degree may be required when
         * gravity is calculated near the surface or inside the model. An increase in the
         * quadrature element adds accuracy to the gravity solution from noise due to the
         * model grid.
         */
        double quadrature_degree_increase;

        /**
         * Gravity may be calculated for a sets of points along the radius (e.g. depth
         * profile) between a minimum and maximum radius. Number of points along the radius
         * is specified with n_points_radius.
         */
        double n_points_radius;

        /**
         * Gravity may be calculated for a sets of points along the longitude (e.g. satellite
         * mapping) between a minimum and maximum longitude. Number of points along the
         * longitude is specified with n_points_longitude.
         */
        double n_points_longitude;

        /**
         * Gravity may be calculated for a sets of points along the latitude (e.g. satellite
         * mapping) between a minimum and maximum latitude. Number of points along the
         * latitude is specified with n_points_latitude.
         */
        double n_points_latitude;

        /**
         * Minimum radius for the depth range in case of the map sampling scheme. Presribe
         * a minimum radius for a sampling coverage at a specific height. May be defined
         * in or outside the model.
         */
        double minimum_radius;

        /**
         * Maximum radius for depth-profile in case of the map sampling scheme. May be
         * defined in or outside the model. No need to specify maximum_radius if
         * n_points_radius is 1.
         */
        double maximum_radius;

        /**
         * Minimum longitude for longitude range in case of the map sampling scheme.
         */
        double minimum_colongitude;

        /**
         * Maximum longitude for the longitude range in case of the map sampling scheme.
         * No need to specify maximum_longitude if n_points_longitude is 1.
         */
        double maximum_colongitude;

        /**
         * Minimum latitude for the latitude range in case of the map sampling scheme.
         */
        double minimum_colatitude;

        /**
         * Maximum latitude for the latitude range in case of the map sampling scheme.
         * No need to specify maximum_latitude if n_points_latitude is 1.
         */
        double maximum_colatitude;

        /**
         * A reference density is required to obtained relative density to calculate
         * gravity anomalies;
         */
        double reference_density;

        /**
         * Specify the sampling scheme determining if gravity calculation is performed
         * for a map of points or a list of points.
         */
        enum SamplingScheme
        {
          map,
          list
        } sampling_scheme;

        /**
         * List of radius coordinates for the list sampling scheme. Must be in order
         * with the lists of longitude and latitude.
         */
        std::vector<double> radius_list;

        /**
         * List of longitude coordinates for the list sampling scheme. Must be in order
         * with the lists of radius and latitude.
         */
        std::vector<double> longitude_list;

        /**
         * List of latitude coordinates for the list sampling scheme. Must be in order
         * with the lists of radius and longitude.
         */
        std::vector<double> latitude_list;

    };
  }
}


#endif
