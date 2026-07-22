/*
  Copyright (C) 2026 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE. If not see
  <http://www.gnu.org/licenses/>.
*/


#ifndef _aspect_postprocess_continental_fragmentation_statistics_h
#define _aspect_postprocess_continental_fragmentation_statistics_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace Postprocess
  {
    /**
     * A postprocessor that computes area-, perimeter-, connectivity-, and
     * velocity-based statistics for continental material on the top surface.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class ContinentalFragmentationStatistics :
      public Interface<dim>,
      public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Evaluate the solution for continental fragmentation statistics.
         */
        std::pair<std::string,std::string>
        execute(TableHandler &statistics) override;

        /**
         * Parse the parameters of this postprocessor.
         */
        void
        parse_parameters(ParameterHandler &prm) override;

        /**
         * Declare the parameters of this postprocessor.
         */
        static void
        declare_parameters(ParameterHandler &prm);

      private:
        /**
         * Properties measured for one top-surface mesh face.
         */
        struct FaceRecord
        {
          /** Center of the surface face. */
          Point<dim> center;

          /** Vertices used to reconstruct surface connectivity and edges. */
          std::vector<Point<dim>> vertices;

          /** Surface measure of the face. */
          double area = 0.0;

          /** Mean sum of the selected continental compositional fields. */
          double continent_value = 0.0;

          /** Mean velocity magnitude on the face in meters per second. */
          double speed = 0.0;

          /** Whether the mean continent value exceeds the continental threshold. */
          bool is_continent = false;

          /**
           * ID of the retained continent block containing this face, or -1 if
           * the face is non-continental or belongs to a filtered block.
           */
          int block_id = -1;
        };

        /**
         * Integrated properties of one retained connected continent block.
         */
        struct BlockSummary
        {
          /** Block identifier for the current timestep. */
          unsigned int id = numbers::invalid_unsigned_int;

          /** Total surface area of the block. */
          double area = 0.0;

          /** Integral of surface speed over the block area. */
          double speed_area_integral = 0.0;

          /** Number of top-surface mesh faces in the block. */
          unsigned int n_faces = 0;

          /** Numerator of the area-weighted Cartesian centroid. */
          Tensor<1,dim> centroid_numerator;
        };

        /**
         * Convert the selected compositional-field names to internal indices.
         */
        void
        resolve_field_indices();

        /**
         * Return whether the specified cell face lies on the model's top boundary.
         */
        bool
        is_top_boundary_face(const typename DoFHandler<dim>::active_cell_iterator &cell,
                             const unsigned int face_no) const;

        /**
         * Write the measured properties and block classification of every
         * top-surface face.
         */
        void
        write_surface_file(const std::vector<FaceRecord> &all_faces) const;

        /**
         * Write the integrated properties and centroid of every retained
         * continent block.
         */
        void
        write_block_file(const std::vector<BlockSummary> &blocks,
                         const std::vector<Point<dim>> &block_centroids) const;

        /**
         * Names of the compositional fields interpreted as continental material.
         */
        std::vector<std::string> continent_field_names;

        /**
         * Internal compositional-field indices corresponding to the selected
         * continental field names.
         */
        std::vector<unsigned int> continent_field_indices;

        /**
         * A surface face is continental when its mean summed continent
         * composition exceeds this value.
         */
        double continent_threshold = 0.5;

        /**
         * Minimum area in square meters for a connected component to be
         * retained as a continent block.
         */
        double minimum_block_area = 0.0;

        /**
         * Whether to include a compact diagnostic summary in the screen output.
         */
        bool output_verbose_screen_line = false;

        /**
         * Whether to write one diagnostic record per top-surface face.
         */
        bool write_surface_map = false;

        /**
         * Whether to write one diagnostic record per retained continent block.
         */
        bool write_block_summary = false;

        /**
         * Prefix used for the surface-map and block-summary file names.
         */
        std::string output_file_prefix = "continental_fragmentation_statistics";
    };
  }
}


#endif
