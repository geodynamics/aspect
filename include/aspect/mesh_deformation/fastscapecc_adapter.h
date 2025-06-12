/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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

#ifndef _aspect_mesh_deformation_fastscapecc_adapter_h
#define _aspect_mesh_deformation_fastscapecc_adapter_h

// #ifndef FASTSCAPELIB_DEALII_GRID_H_
// #define FASTSCAPELIB_DEALII_GRID_H_

#include <unordered_map>
#include <vector>

#include <xtensor/xbroadcast.hpp>
#include <xtensor/xmath.hpp>

#include "fastscapelib/grid/base.hpp"
#include "fastscapelib/utils/xtensor_utils.hpp"


namespace fastscapelib
{
  template <class T, class S>
  class dealii_grid;

  /**
   * Dealii grid specialized types
   */
  template <class T, class S>
  struct grid_inner_types<dealii_grid<T, S>>
  {
    static constexpr bool is_structured = false;
    static constexpr bool is_uniform = false;

    using grid_data_type = double;

    using xt_selector = S;
    static constexpr std::size_t xt_ndims = 1;

    static constexpr uint8_t n_neighbors_max = 8;
    using neighbors_cache_type = neighbors_no_cache<n_neighbors_max>;
  };

  /**
   * @brief Grid adapter class for using a DEAL.II surface or boundary mesh with
   * Fastscapelib
   *
   * This adapter should work with various kinds of DEAL.II meshes (e.g., box,
   * sphere, etc.). The only requirement is that the mesh must be uniform (i.e.,
   * all cells have an equal area) at a given refinement level.
   *
   * Note: a ``dealii_grid`` node cooresponds to a DEAL.II mesh cell.
   *
   * @tparam T The DEAL.II mesh type.
   * @tparam S The container selector for data array members.
   */
  template <class T, class S = xt_selector>
  class dealii_grid : public grid<dealii_grid<T, S>>
  {
    public:
      using self_type = dealii_grid<T, S>;
      using base_type = grid<self_type>;
      using inner_types = grid_inner_types<self_type>;

      using grid_data_type = typename base_type::grid_data_type;

      using xt_selector = typename base_type::xt_selector;
      using container_type = xt_tensor_t<xt_selector, grid_data_type, inner_types::xt_ndims>;
      using neighbors_type = typename base_type::neighbors_type;
      using neighbors_indices_type = typename base_type::neighbors_indices_type;
      using neighbors_distances_type = typename base_type::neighbors_distances_type;

      using nodes_status_type = typename base_type::nodes_status_type;
      using nodes_status_array_type = xt_tensor_t<xt_selector, node_status, 1>;

      using size_type = typename base_type::size_type;
      using shape_type = typename base_type::shape_type;

      dealii_grid(T &triangulation, double cell_area);

      void set_nodes_status(const nodes_status_array_type &nodes_status);

    protected:
      T &m_triangulation;

      shape_type m_shape;
      size_type m_size;
      double m_node_area;

      nodes_status_type m_nodes_status;

      using neighbors_distances_impl_type = typename base_type::neighbors_distances_impl_type;
      using neighbors_indices_impl_type = typename base_type::neighbors_indices_impl_type;

      std::vector<size_type> m_neighbors_count;
      std::vector<neighbors_indices_impl_type> m_neighbors_indices;
      std::vector<neighbors_distances_impl_type> m_neighbors_distances;

      void compute_connectivity();

      inline container_type nodes_areas_impl() const;
      inline grid_data_type nodes_areas_impl(const size_type &idx) const noexcept;

      inline size_type neighbors_count_impl(const size_type &idx) const;

      void neighbors_indices_impl(neighbors_indices_impl_type &neighbors,
                                  const size_type &idx) const;

      inline const neighbors_distances_impl_type &neighbors_distances_impl(
        const size_type &idx) const;

      static constexpr std::size_t dimension_impl() noexcept;

      friend class grid<self_type>;
  };


  /**
   * Creates a new adapter from an existing DEAL.II surface mesh.
   *
   * @param triangulation The DEAL.II triangulation object.
   * @param cell_area the uniform area of the mesh cells.
   */
  template <class T, class S>
  dealii_grid<T, S>::dealii_grid(T &triangulation, double cell_area)
    : base_type(0)
    , m_triangulation(triangulation)
    , m_node_area(cell_area)
  {

    m_size = triangulation.n_global_active_cells();
    m_shape = { static_cast<typename shape_type::value_type>(m_size) };

    // pre-compute explicit grid connectivity
    compute_connectivity();

    // Set all grid nodes as "active"
    // - Fastscapelib doesn't support yet labelling grid nodes as ghost nodes
    // - we assume that all the erosion processes that we apply here depend on flow routing
    //   (i.e., they rely on a `fastscapelib::flow_graph` object built on top of this grid,
    //   such as `fastscapelib::spl_eroder`)
    // - It is better (more explicit) to set base level nodes via the flow graph
    //   object (i.e., `fastscapelib::flow_graph::set_base_levels`)
    // - Masks (e.g., oceans) can be defined via the flow graph too
    //   (i.e., `fastscapelib::flow_graph::set_mask`)
    // - TODO: periodic boundary conditions (if any) should be handled in `compute_connectivity`
    //   (there is no need to label nodes as periodic boundaries since the grid connectivity
    //   is computed explicitly here)
    // m_nodes_status = node_status::core;
    m_nodes_status = xt::xarray<node_status>::from_shape({m_size});
    std::fill(m_nodes_status.begin(), m_nodes_status.end(), node_status::core);


  }

  template <class T, class S>
  void dealii_grid<T, S>::set_nodes_status(const nodes_status_array_type &nodes_status)
  {
    if (!xt::same_shape(nodes_status.shape(), m_shape))
      {
        throw std::invalid_argument(
          "invalid shape for nodes_status array (expects shape [N] where N is the total number of nodes)");
      }
    m_nodes_status = nodes_status;
  }

//   template <class T, class S>
//   void dealii_grid<T, S>::compute_connectivity()
//   {
//     m_neighbors_count.resize(m_size);
//     m_neighbors_indices.resize(m_size);
//     m_neighbors_distances.resize(m_size);
//     unsigned int counter = 0;
//     std::unordered_map<unsigned int, unsigned int> global_to_local_index;
//     unsigned int local_index = 0;
//     // for (const auto &cell : m_triangulation.get_dof_handler().active_cell_iterators())
//     for (const auto &cell : m_triangulation.active_cell_iterators())
//       {
//         m_neighbors_count[counter] = cell->n_faces();
//         auto center = cell.center();
//         for (auto face_i = 0; face_i < m_neighbors_count[counter]; ++face_i)
//           {
//             unsigned int global_neighbor_index = 0;
//             for (const auto &neighbor_cell : cell.neighbor(face_i))
//               {
//                 global_neighbor_index = neighbor_cell.global_active_cell_index();
//                 m_neighbors_distances[counter][face_i] = center.distance(neighbor_cell.center());
//                 // TODO: I assume we only have one neighbor per face, since we start with a uniform grid, we will only have one, but cwe could generatlize this
//               }
//             // TODO: add assert to check if neighbor_index is correclty set
//             global_to_local_index.emplace(global_neighbor_index, counter);
//           }

//         counter++;
//       }
//     counter = 0;
//     // for (const auto &cell : this->get_dof_handler().active_cell_iterators())
//     for (const auto &cell : m_triangulation.active_cell_iterators())
//       {
//         for (auto k = 0; k <= m_neighbors_count[counter]; k++)
//           {
//             auto global_nb_index = m_neighbors_indices[counter][k];
//             m_neighbors_indices[counter][k] = global_to_local_index.at(global_nb_index);
//           }
//         counter++;
//       }
//   }


        template <class T, class S>
        void fastscapelib::dealii_grid<T, S>::compute_connectivity()
        {
        // Resize internal storage for neighbor counts, indices, and distances
        m_neighbors_count.resize(m_size);
        m_neighbors_indices.resize(m_size);
        m_neighbors_distances.resize(m_size);

        // Step 1: Map global active cell indices to local indices
        // This is necessary to translate global indices into local array indices
        std::unordered_map<unsigned int, unsigned int> global_to_local_index;
        unsigned int local_index = 0;
        for (const auto &cell : m_triangulation.active_cell_iterators())
        {
            global_to_local_index[cell->global_active_cell_index()] = local_index++;
        }

        // Step 2: Loop through each active cell and determine valid neighbors
        local_index = 0;
        for (const auto &cell : m_triangulation.active_cell_iterators())
        {
            const auto &center = cell->center();  // Get cell center for distance computation
            unsigned int valid_neighbors = 0;

            for (unsigned int face_i = 0; face_i < cell->n_faces(); ++face_i)
            {
            // Only consider internal faces with valid neighbor
            // if (!cell->at_boundary(face_i) && cell->neighbor_index_is_valid(face_i))
            if (!cell->at_boundary(face_i) && cell->neighbor(face_i).state() == dealii::IteratorState::valid)
            {
                const auto &neighbor = cell->neighbor(face_i);
                const unsigned int global_nb_index = neighbor->global_active_cell_index();

                auto it = global_to_local_index.find(global_nb_index);
                    AssertThrow(it != global_to_local_index.end(),
                                aspect::ExcMessage("Neighbor index not found in global_to_local_index map."));


                // Fill neighbor index and distance arrays
                m_neighbors_indices[local_index][valid_neighbors] = it->second;
                m_neighbors_distances[local_index][valid_neighbors] = center.distance(neighbor->center());

                ++valid_neighbors;
            }
            }

            // Save how many neighbors this node has
            m_neighbors_count[local_index] = valid_neighbors;

            ++local_index;
        }
        }


  template <class T, class S>
  inline auto dealii_grid<T, S>::nodes_areas_impl() const -> container_type
  {
    return xt::broadcast(m_node_area, m_shape);
  }

  template <class T, class S>
  inline auto dealii_grid<T, S>::nodes_areas_impl(const size_type & /*idx*/) const noexcept
  -> grid_data_type
  {
    return m_node_area;
  }

  template <class T, class S>
  inline auto dealii_grid<T, S>::neighbors_count_impl(const size_type &idx) const -> size_type
  {
    return m_neighbors_count[idx];
  }

  template <class T, class S>
  void dealii_grid<T, S>::neighbors_indices_impl(neighbors_indices_impl_type &neighbors,
                                                 const size_type &idx) const
  {
    const auto &size = m_neighbors_count[idx];

    for (size_type i = 0; i < size; i++)
      {
        neighbors[i] = m_neighbors_indices[idx][i];
      }
  }

  template <class T, class S>
  auto dealii_grid<T, S>::neighbors_distances_impl(const size_type &idx) const
  -> const neighbors_distances_impl_type &
  {
    return m_neighbors_distances[idx];
  }
}

#endif  // _aspect_mesh_deformation_fastscapecc_adapter_h
