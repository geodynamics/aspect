/*
  Copyright (C) 2018 - 2020 by the authors of the ASPECT code.

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
     * of point coordinates. Spherical coordinates in the output file are radius,
     * colatitude and colongitude. Gravity is here based on the density distribution
     * from the material model (and non adiabatic). This means that the density may
     * come directly from an ascii file. This postprocessor also computes theoretical
     * gravity and its derivatives, which corresponds to the analytical solution of gravity
     * in the same geometry but filled with a reference density. The reference density
     * is also used to determine the density difference for computing gravity anomalies.
     * Thus one must carefully evaluate the meaning of the gravity anomaly output,
     * because the solution may not reflect the actual gravity anomaly (due to
     * differences in the assumed reference density). One way to guarantee correct
     * gravity anomalies is to subtract the gravity of a certain point from the average
     * gravity on the map. Another way is to directly use density anomalies for this
     * postprocessor.
     * The average- minimum- and maximum gravity acceleration and potential are written
     * into the statistics file.

     * @ingroup Postprocessing
     */
    template <int dim>
    class GravityPointValues : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        GravityPointValues ();

        /**
         * Specify the creation of output_gravity.txt.
         */
        std::pair<std::string,std::string> execute (TableHandler &) override;

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

        /**
         * Serialize the contents of this class as far as they are not read
         * from input parameter files.
         */
        template <class Archive>
        void serialize (Archive &ar, const unsigned int version);

        /**
         * Save the state of this object.
         */
        void save (std::map<std::string, std::string> &status_strings) const override;

        /**
         * Restore the state of the object.
         */
        void load (const std::map<std::string, std::string> &status_strings) override;

      private:
        /**
         * Interval between the generation of gravity output.
         */
        double output_interval;

        /**
         * A time (in seconds) at which the last graphical output was supposed
         * to be produced. Used to check for the next necessary output time.
         */
        double last_output_time;

        /**
         * Maximum number of steps between the generation of gravity output.
         */
        unsigned int maximum_timesteps_between_outputs;

        /**
         * Timestep at which the last gravity output was produced
         * Used to check for the next necessary output time.
         */
        unsigned int last_output_timestep;

        /**
         * Consecutively counted number indicating the how-manyth time we will
         * create output the next time we get to it.
         */
        unsigned int output_file_number;

        /**
         * end_time is taken from parameter file. It is used to tell the
         * postprocessor to write gravity output at the end time.
         */
        double end_time;

        /**
         * Set the time output was supposed to be written. In the simplest
         * case, this is the previous last output time plus the interval, but
         * in general we'd like to ensure that it is the largest supposed
         * output time, which is smaller than the current time, to avoid
         * falling behind with last_output_time and having to catch up once
         * the time step becomes larger. This is done after every output.
         */
        void set_last_output_time (const double current_time);

        /**
         * Set the precision of the gravity acceleration, potential and gradients
         * in the gravity output and statistics file.
         */
        unsigned int precision;

        /**
         * Quadrature degree increase over the velocity element degree may be required when
         * gravity is calculated near the surface or inside the model. An increase in the
         * quadrature element adds accuracy to the gravity solution from noise due to the
         * model grid.
         */
        unsigned int quadrature_degree_increase;

        /**
         * Parameter for the fibonacci spiral sampling scheme:
         */
        unsigned int n_points_spiral;

        /**
         * Parameter for the map and fibonacci spiral sampling scheme:
         * Gravity may be calculated for a sets of points along the radius (e.g. depth
         * profile) between a minimum and maximum radius. Number of points along the radius
         * is specified with n_points_radius.
         */
        unsigned int n_points_radius;

        /**
         * Parameter for the map sampling scheme:
         * Gravity may be calculated for a sets of points along the longitude (e.g. satellite
         * mapping) between a minimum and maximum longitude. Number of points along the
         * longitude is specified with n_points_longitude.
         */
        unsigned int n_points_longitude;

        /**
         * Parameter for the map sampling scheme:
         * Gravity may be calculated for a sets of points along the latitude (e.g. satellite
         * mapping) between a minimum and maximum latitude. Number of points along the
         * latitude is specified with n_points_latitude.
         */
        unsigned int n_points_latitude;

        /**
         * Parameter for the map and fibonacci spiral sampling scheme:
         * Prescribe a minimum radius for a sampling coverage at a specific height.
         * May be set in- or outside the model domain.
         */
        double minimum_radius;

        /**
         * Parameter for the map and fibonacci spiral sampling scheme:
         * Maximum radius for the radius range.
         * May be set in- or outside the model domain.
         * No need to specify maximum_radius if n_points_radius is 1.
         */
        double maximum_radius;

        /**
         * Parameter for the map sampling scheme:
         * Minimum longitude for longitude range.
         */
        double minimum_colongitude;

        /**
         * Parameter for the map sampling scheme:
         * Maximum longitude for the longitude range.
         * No need to specify maximum_longitude if n_points_longitude is 1.
         */
        double maximum_colongitude;

        /**
         * Parameter for the map sampling scheme:
         * Minimum latitude for the latitude range.
         */
        double minimum_colatitude;

        /**
         * Parameter for the map sampling scheme:
         * Maximum latitude for the latitude range.
         * No need to specify maximum_latitude if n_points_latitude is 1.
         */
        double maximum_colatitude;

        /**
         * A reference density is required for two purposes:
         * 1) benchmark the gravity postprocessor with computing analytically gravity
         * acceleration and gradients in the domain filled with a constant density.
         * 2) calculate the density difference of the density given by the material
         * model and a constant density, in order to compute gravity anomalies.
         */
        double reference_density;

        /**
         * Specify the sampling scheme determining if gravity calculation is performed.
         */
        enum SamplingScheme
        {
          map,
          fibonacci_spiral,
          list_of_points
        } sampling_scheme;

        /**
         * Parameter for the list of points sampling scheme:
         * List of radius coordinates for the list of points sampling scheme.
         * Must follow the same order as the lists of longitude and latitude.
         */
        std::vector<double> radius_list;

        /**
         * Parameter for the list of points sampling scheme:
         * List of longitude coordinates for the list of points sampling scheme.
         * Must follow the same order as the lists of radius and latitude.
         */
        std::vector<double> longitude_list;

        /**
         * Parameter for the list of points sampling scheme:
         * List of latitude coordinates for the list of points sampling scheme.
         * Must follow the same order as the lists of longitude and longitude.
         */
        std::vector<double> latitude_list;

    };
  }
}


#endif
