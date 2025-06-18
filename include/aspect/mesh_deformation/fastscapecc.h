/*
  Copyright (C) 2023 - 2025 by the authors of the ASPECT code.

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

#ifdef ASPECT_WITH_FASTSCAPELIB

#include <type_traits>

#include <aspect/global.h>
#include <aspect/geometry_model/box.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/mesh_deformation/interface.h>
#include <aspect/mesh_deformation/fastscapecc_adapter.h>
#include <aspect/simulator_access.h>

#include <fastscapelib/flow/flow_graph.hpp>
#include <fastscapelib/flow/flow_router.hpp>
#include <fastscapelib/flow/sink_resolver.hpp>
#include <fastscapelib/eroders/diffusion_adi.hpp>
#include <fastscapelib/eroders/spl.hpp>

#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>


#include <deal.II/lac/affine_constraints.h>

namespace aspect
{
  using namespace dealii;

  namespace MeshDeformation
  {
    /**
     * A comparator for dealii::Point<dim> to use in std::map
     * ------------------------------------------------------
     * std::map requires a way to compare keys (i.e., define a "less than" relation)
     * dealii::Point<dim> does not provide operator< by default,
     * so we must define one ourselves.
     *
     * This comparator compares two Point<dim> objects lexicographically:
     * First compares x, then y (and z if dim == 3).
     *
     * This is required to allow using std::map<Point<dim>, T>
     * to uniquely index surface vertices on unstructured spherical meshes.
     */
    template <int dim>
    struct PointComparator
    {
      bool operator()(const dealii::Point<dim> &a, const dealii::Point<dim> &b) const
      {
        for (unsigned int d = 0; d < dim; ++d)
          {
            if (a[d] < b[d]) return true;
            if (a[d] > b[d]) return false;
          }
        return false;
      }
    };

    /**
     * A plugin that utilizes the landscape evolution code Fastscapelib (C++)
     * to deform the ASPECT boundary through erosion.
     */
    template<int dim>
    class FastScapecc : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        FastScapecc();

        /**
         * Initialize variables for FastScape.
         */
        void initialize () override;

        struct ValuesAtSurfaceVertex
        {
          double topography;
          double surface_uplift_rate;
        };

        /**
          * A function that creates constraints for the velocity of certain mesh
          * vertices (namely, the surface vertices) for a specific boundary.
          * The calling class will respect
          * these constraints when computing the new vertex positions.
          */
        void
        compute_velocity_constraints_on_boundary(const DoFHandler<dim> &mesh_deformation_dof_handler,
                                                 AffineConstraints<double> &mesh_velocity_constraints,
                                                 const std::set<types::boundary_id> &boundary_id) const override;


        double interpolate_surface_velocity(const Point<dim> &p,
                                            const std::vector<double> &V) const;
        /**
         * Returns whether or not the plugin requires surface stabilization.
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
        void parse_parameters (ParameterHandler &prm) override;


      private:

        /**
         * For a given location, typically the location of a vertex of ASPECT's
         * volume mesh, return the DoF index that corresponds to this vertex
         * location within the surface_dof_handler. This index can then also
         * be used to index into vectors defined on surface_dof_handler.
         */
        unsigned int vertex_index(const Point<dim> &p) const;

        /*
         * A map, used only for the spherical geometry, that maps the location
         * of a vertex to an index between zero and the number of vertices
         * in the surface mesh.
         *
         * TODO: This map is only used for the spherical shell geometry. We
         *   should prevent this object from being used altogether in all
         *   other cases, and this could be facilitated by wrapping it in
         *   std_cxx17::optional<>. Only when we encounter the spherical
         *   geometry do we actually set anything in the optional.
         */
        std::map<dealii::Point<dim>, unsigned int, PointComparator<dim>> spherical_vertex_index_map;

        /**
         * The surface mesh is a distributed Triangulation to support large-scale
         * parallel simulations (e.g., global models).
         */
        using SurfaceMeshType = parallel::distributed::Triangulation<dim - 1, dim>;

        SurfaceMeshType surface_mesh;
        DoFHandler<dim - 1, dim> surface_mesh_dof_handler;
        mutable AffineConstraints<double> surface_constraints;

        FE_Q<dim - 1, dim> surface_fe;

        /**
         * Initialize the surface mesh based on the given geometry model.
         */
        void init_surface_mesh(const GeometryModel::Interface<dim> &);

        /**
         * Pointers to Fastscapelib objects
         */
        using GridAdapterType = typename fastscapelib::dealii_grid<SurfaceMeshType>;
        using FlowGraphType = typename fastscapelib::flow_graph<GridAdapterType>;
        std::unique_ptr<GridAdapterType> grid;
        std::unique_ptr<FlowGraphType> flow_graph;
        std::unique_ptr<fastscapelib::spl_eroder<FlowGraphType>> spl_eroder;
        // mutable std::shared_ptr<dealii::GridTools::Cache<dim-1, dim-1>> surface_cache;


        void project_surface_solution(const std::set<types::boundary_id> &boundary_ids,
                                      dealii::LinearAlgebra::distributed::Vector<double> &surface_vertical_velocity,
                                      dealii::LinearAlgebra::distributed::Vector<double> &surface_elevation) const;


        /**
         * Project the Stokes velocity solution onto the
         * free surface. Called by make_constraints()
         */
        // void project_velocity_onto_boundary (const DoFHandler<dim> &free_surface_dof_handler,
        //                                      const IndexSet &mesh_locally_owned,
        //                                      const IndexSet &mesh_locally_relevant,
        //                                      LinearAlgebra::Vector &output) const;


        /**
         * Fill velocity data table to be interpolated back onto the ASPECT mesh.
         */
        // Table<dim,double> fill_data_table(std::vector<double> &values,
        //                                   TableIndices<dim> &size_idx,
        //                                   const int &array_size) const;

        /**
         * Run-time parameters
         * @{
         */
        /**
         * Suggestion for the number of FastScape steps to run for every ASPECT timestep,
         * where the FastScape time step is determined by ASPECT's time step length divided by
         * this parameter.
         */
        unsigned int fastscape_steps_per_aspect_step;

        /**
         * Maximum time step allowed for FastScape. If the suggested time step exceeds this
         * limit it is repeatedly divided by 2 until the final time step is smaller than this parameter.
         */
        double maximum_fastscape_timestep;

        /**
         * Total number of Fastscapelib grid nodes.
         */
        unsigned int n_grid_nodes;

        /**
         * How many levels FastScape should be refined above the maximum ASPECT surface resolution.
         */
        unsigned int additional_refinement_levels;

        /**
         * Maximum expected refinement level at ASPECT's surface.
         * This and resolution_difference are required to properly transfer node data from
         * ASPECT to FastScape.
         */
        unsigned int surface_refinement_level;

        /**
         * Difference in refinement levels expected at the ASPECT surface,
         * where this would be set to 2 if 3 refinement levels are set at the surface.
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
         * Magnitude (m) of the initial noise applied to FastScape.
         * Applied as either a + or - value to the topography
         * such that the total difference can be up to 2*noise_h.
         */
        double noise_h;

        // FastScape erosional parameters //
        /**
         * Drainage area exponent for the stream power law (m variable in FastScape surface equation).
         */
        double area_exp;

        /**
         * Slope exponent for the steam power law (n variable in FastScape surface equation).
         */
        double slope_exp;

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
         * @}
         */


        // Grid extent in each direction [min, max]
        // std::array<std::pair<double, double>, 2> grid_extent_surface;

        std::array<std::pair<double,double>,dim> box_grid_extent;


        // Grid resolution (spacing between points)
        double dx;
        double dy;

        // Number of grid cells in x and y directions
        unsigned int nx;
        unsigned int ny;

        mutable bool elevation_initialized = false;
        mutable xt::xarray<double> elevation;

        /**
         * The number of cells in each coordinate direction.
         */
        std::array<unsigned int, dim> repetitions;
        std::array<std::pair<double, double>, dim - 1> grid_extent_surface;

        // A way to ouput internal fastscape iteration for visualisation
        bool output_internal_fastscape_steps;


    };
  }
}

#endif // #ifdef ASPECT_WITH_FASTSCAPELIB

#endif
