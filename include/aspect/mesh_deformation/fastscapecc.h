/*
  Copyright (C) 2023 by the authors of the ASPECT code.
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

#ifndef _aspect_mesh_deformation_fastscapecc_h
#define _aspect_mesh_deformation_fastscapecc_h

#include <aspect/global.h>

#include <aspect/mesh_deformation/interface.h>
#include <aspect/simulator_access.h>

#include <xtensor/xarray.hpp>
#include <xtensor/xmath.hpp>
#include <xtensor/xview.hpp>
#include <xtensor/xadapt.hpp>

#include <fastscapelib/flow/flow_graph.hpp>
#include <fastscapelib/flow/sink_resolver.hpp>
#include <fastscapelib/flow/flow_router.hpp>
#include <fastscapelib/grid/raster_grid.hpp>
#include <fastscapelib/eroders/diffusion_adi.hpp>
#include <fastscapelib/eroders/spl.hpp>

#include <fastscapelib/grid/healpix_grid.hpp>

namespace aspect
{
  using namespace dealii;

  namespace MeshDeformation
  {

    /**
     * A plugin that utilizes the landscape evolution code FastScape
     * to deform the ASPECT boundary through advection, uplift,
     * hillslope diffusion, sediment deposition, marine diffusion,
     * and the stream power law, which describes river incision.
     *
     */
    template<int dim>
    class FastScapecc : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        enum class GeometryType
        {
          Box, SphericalShell, Undefined
        };
        /**
         * Initialize variables for FastScape.
         */
        virtual void initialize ();

        /**
          * A function that creates constraints for the velocity of certain mesh
          * vertices (e.g. the surface vertices) for a specific boundary.
          * The calling class will respect
          * these constraints when computing the new vertex positions.
          */
        virtual
        void
        compute_velocity_constraints_on_boundary(const DoFHandler<dim> &mesh_deformation_dof_handler,
                                                 AffineConstraints<double> &mesh_velocity_constraints,
                                                 const std::set<types::boundary_id> &boundary_id) const;


        /**
         * Returns whether or not the plugin requires surface stabilization
         */
        bool needs_surface_stabilization () const override;

        /**
         * Declare parameters for the FastScape plugin.
         */
        static
        void declare_parameters (ParameterHandler &prm);

        /**
         * Parse parameters for the FastScape plugin.
         */
        void parse_parameters (ParameterHandler &prm);

      private:
        GeometryType geometry_type;

        // // Add unique pointers for FastScape components
        // std::unique_ptr<fastscapelib::healpix_grid<>> grid; // Healpix grid pointer

        // // Unique pointer for fastscapelib flow_graph with healpix_grid as template parameter.
        // std::unique_ptr<fastscapelib::flow_graph<fastscapelib::healpix_grid<>>> flow_graph;

        // // Unique pointer for fastscapelib spl_eroder with flow_graph as template parameter.
        // std::unique_ptr<fastscapelib::spl_eroder<fastscapelib::flow_graph<fastscapelib::healpix_grid<>>>> spl_eroder;



        // Unique pointer for fastscapelib raster_boundary_status.
        std::unique_ptr<fastscapelib::raster_boundary_status> bs;

        // Unique pointer for fastscapelib raster_grid.
        std::unique_ptr<fastscapelib::raster_grid<>> grid_box;

        // Unique pointer for fastscapelib flow_graph with raster_grid as template parameter.
        std::unique_ptr<fastscapelib::flow_graph<fastscapelib::raster_grid<>>> flow_graph_box;

        // Unique pointer for fastscapelib spl_eroder with flow_graph<raster_grid> as template parameter.
        std::unique_ptr<fastscapelib::spl_eroder<fastscapelib::flow_graph<fastscapelib::raster_grid<>>>> spl_eroder_box;

        // Unique pointer for fastscapelib diffusion_adi_eroder with raster_grid as template parameter.
        std::unique_ptr<fastscapelib::diffusion_adi_eroder<fastscapelib::raster_grid<>>> diffusion_eroder_box;

        /**
         * Fill velocity data table to be interpolated back onto the ASPECT mesh.
         */
        Table<dim,double> fill_data_table(std::vector<double> &values,
                                          TableIndices<dim> &size_idx,
                                          const int &array_size) const;


        /**
         * Suggestion for the number of FastScape steps to run for every ASPECT timestep,
         * where the FastScape timestep is determined by ASPECT_timestep_length divided by
         * this parameter.
         */
        unsigned int fastscape_steps_per_aspect_step;

        /**
         * Maximum timestep allowed for FastScape, if the suggested timestep exceeds this
         * limit it is repeatedly divided by 2 until the final timestep is smaller than this parameter.
         */
        double maximum_fastscape_timestep;

        /**
         * Maximum timestep allowed for FastScape, if the suggested timestep exceeds this
         * limit it is repeatedly divided by 2 until the final timestep is smaller than this parameter.
         */

        /**
         * Check whether FastScape needs to be restarted.
         */
        mutable bool restart;

        /**
         * ASPECT end time to check if we should destroy FastScape.
         */
        double end_time;

        /**
         * FastScape cell size in X, dx should always equal dy.
         */
        double dx;

        /**
         * FastScape cell size in Y, dy should always equal dx.
         */
        double dy;

        /**
         * FastScape X extent (ASPECT X extent + 2*dx for ghost nodes).
         */
        double x_extent;

        /**
         * Fastscape Y extent (ASPECT Y extent + 2*dy for ghost nodes).
         */
        double y_extent;

        /**
         * User set FastScape Y extent for a 2D ASPECT model.
         */
        double y_extent_2d;

        /**
         * Number of x points in FastScape array.
         */
        int nx;

        /**
         * Number of y points in FastScape array.
         */
        int ny;

        /**
         * Size of the FastScape array (nx*ny).
         */
        unsigned int array_size;

        /**
         * Number of faces for the healpix grid.
         */
        int nsides;

        /**
         * How many levels FastScape should be refined above the maximum ASPECT surface resolution.
         */
        unsigned int additional_refinement_levels;

        /**
         * Maximum expected refinement level at ASPECT's surface.
         * This and resolution_difference are required to properly transfer node data from
         * ASPECT to FastScape.
         */
        int maximum_surface_refinement_level;

        /**
         * Difference in refinement levels expected at the ASPECT surface,
         * where this would be set to 2 if 3 refinement leves are set at the surface.
         * This and surface_resolution are required to properly transfer node data from
         * ASPECT to FastScape.
         *
         * TODO: Should this be kept this way, or make it so the input is the expected levels
         * of refinement at the surface, and we can subtract one within the code? Also,
         * it would be good to find a way to check these are correct, because they are a
         * common source of errors.
         */
        int surface_refinement_difference;

        /**
         * If set to false, the FastScape surface is averaged along Y and returned
         * to ASPECT. If set to true, the center slice of the FastScape model is
         * returned to ASPECT.
         *
         * TODO: Do we average the ghost nodes or remove them? Need to double check.
         */
        bool center_slice;

        /**
         * Seed number for initial topography noise in FastScape.
         */
        int fs_seed;

        /**
         * Variable to hold ASPECT domain extents.
         */
        std::array<std::pair<double,double>,dim> grid_extent;

        /**
         * Table for interpolating FastScape surface velocities back to ASPECT.
         */
        std::array< unsigned int, dim > table_intervals;


        /**
         * Magnitude (m) of the initial noise applied to FastScape.
         * Applied as either a + or - value to the topography
         * such that the total difference can be up to 2*noise_h.
         */
        double noise_h;


        // FastScape boundary condition variables //
        /**
         * FastScape bottom boundary condition that determines topography at the FastScape bottom boundary.
         * Where 1 represents a fixed height boundary (though this can still be uplifted through uplift velocities), and 0 a
         * reflective boundary. When two opposing boundaries are reflective (e.g., top and bottom are both zero), then the boundaries
         * become cyclic.
         */
        unsigned int bottom;

        /**
         * FastScape top boundary condition that determines topography at the FastScape top boundary.
         * Where 1 represents a fixed height boundary (though this can still be uplifted through uplift velocities), and 0 a
         * reflective boundary. When two opposing boundaries are reflective (e.g., top and bottom are both zero), then the boundaries
         * become cyclic.
         */
        unsigned int top;

        /**
         * FastScape right boundary condition that determines topography at the FastScape right boundary.
         * Where 1 represents a fixed height boundary (though this can still be uplifted through uplift velocities), and 0 a
         * reflective boundary. When two opposing boundaries are reflective (e.g., left and right are both zero), then the boundaries
         * become cyclic.
         */
        unsigned int right;

        /**
         * FastScape left boundary condition that determines topography at the FastScape left boundary.
         * Where 1 represents a fixed height boundary (though this can still be uplifted through uplift velocities), and 0 a
         * reflective boundary. When two opposing boundaries are reflective (e.g., left and right are both zero), then the boundaries
         * become cyclic.
         */
        unsigned int left;

        /**
         * Integer that holds the full boundary conditions sent to FastScape (e.g., 1111).
         */
        int bc;

        // FastScape erosional parameters //
        /**
         * Drainage area exponent for the stream power law (m variable in FastScape surface equation).
         */
        double m;

        /**
         * Slope exponent for the steam power law (n variable in FastScape surface equation).
         */
        double n;

        /**
         * Slope exponent for multi-direction flow, where 0 is uniform, and 10 is steepest descent. (-1 varies with slope)
         * (p variable in FastScape surface equation).
         */
        double p;


        /**
         * Bedrock river incision rate for the stream power law
         * (meters^(1-2m)/yr, kf variable in FastScape surface equation).
         */
        double kff;

        /**
         * Sediment river incision rate for the stream power law (meters^(1-2m)/yr).
         * When set to -1 this is identical to the bedrock value.
         * (kf variable in FastScape surface equation applied to sediment).
         */
        double kfsed;

        /**
         * Bedrock transport coefficient for hillslope diffusion (m^2/yr, kd in FastScape surface equation).
         */
        double kdd;

        /**
         * Bedrock transport coefficient for hillslope diffusion (m^2/yr). When set to -1 this is
         * identical to the bedrock value.
         * (kd in FastScape surface equation applied to sediment).
         */
        double kdsed;

        /**
         * Precision value for how close a ASPECT node must be to the FastScape node
         * for the value to be transferred. This is only necessary if use_v is set to 0
         * and the free surface is used to advect the surface with a normal projection.
         */
        double node_tolerance;

        /**
         * The number of cells in each coordinate direction.
         */
        std::array<unsigned int, dim> repetitions;


        /**
         * FastScape X extent (ASPECT X extent + 2*dx for ghost nodes).
         */
        double precision;

        double inner_radius;
        double outer_radius;
        double opening_angle;
    };
  }
}

#endiff
