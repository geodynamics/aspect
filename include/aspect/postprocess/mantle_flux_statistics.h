/*
  Copyright (C) 2011 - 2026 by the authors of the ASPECT code.

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

#ifndef _aspect_postprocess_mantle_flux_statistics_h
#define _aspect_postprocess_mantle_flux_statistics_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

#include <limits>

namespace aspect
{
  namespace Postprocess
  {
    /**
     * A postprocessor that measures hot upwellings and cold downwellings
     * across layers at selected depths in a spherical model.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class MantleFluxStatistics : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Values used to identify and visualize a plume or slab at one point.
         */
        struct PointDiagnostics
        {
          double outward_velocity = 0.0;
          double temperature_anomaly = 0.0;
          double temperature_flux_density = 0.0;
          double buoyancy_mass_flux_density = 0.0;
          double buoyancy_force_rate_density = 0.0;
          int structure = 0;
        };

        /**
         * Compute fluxes and measurements at all requested depths.
         */
        std::pair<std::string,std::string>
        execute (TableHandler &statistics) override;

        /**
         * Return the local quantities used by this postprocessor. A structure
         * value of 1 marks a plume, -1 marks a slab, and 0 marks background
         * mantle.
         */
        PointDiagnostics
        evaluate_point (const Point<dim> &position,
                        const Tensor<1,dim> &velocity,
                        const double temperature,
                        const double density,
                        const double thermal_expansivity,
                        const Tensor<1,dim> &gravity) const;

        /**
         * Declare the parameters this postprocessor reads.
         */
        static void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters from the input file.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        /**
         * Join nearby points into separate mantle structures.
         */
        std::vector<std::vector<Point<dim>>>
        find_clusters (const std::vector<Point<dim>> &points) const;

        /**
         * Return the largest surface distance between two points in a cluster.
         */
        double
        cluster_length (const std::vector<Point<dim>> &points) const;

        /**
         * Return the mean surface distance from the cluster center.
         */
        double
        cluster_radius (const std::vector<Point<dim>> &points) const;

        /**
         * Surface distance between two points in a spherical model.
         */
        double
        surface_distance (const Point<dim> &first,
                          const Point<dim> &second) const;

        /**
         * Write the centers of cells assigned to each plume and slab.
         */
        void
        write_cluster_file (const std::string &content);

        std::vector<double> measurement_depths;
        double measurement_half_thickness = 20000.0;
        double cold_temperature_threshold = -200.0;
        double hot_temperature_threshold = 200.0;
        double minimum_slab_length = 0.0;
        double maximum_point_spacing = 80000.0;
        unsigned int minimum_points_per_structure = 5;
        bool write_cluster_files = false;
        double cluster_file_interval = 0.0;
        double last_output_time = std::numeric_limits<double>::quiet_NaN();
    };
  }
}

#endif
