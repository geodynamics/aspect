#ifndef FASTSCAPELIB_DEALII_GRID_H_
#define FASTSCAPELIB_DEALII_GRID_H_

#include <unordered_map>

#include <xtensor/xbroadcast.hpp>
#include <xtensor/xmath.hpp>

#include "fastscapelib/grid/base.hpp"
#include "fastscapelib/utils/xtensor_utils.hpp"

#include <tuple>
#include <utility>
#include <vector>


namespace fastscapelib
{
  template <class S, class T>
  class dealii_grid;

  /**
   * Dealii grid specialized types
   */
  template <class S, class T>
  struct grid_inner_types<dealii_grid<S, T>>
  {
    static constexpr bool is_structured = false;
    static constexpr bool is_uniform = false;

    using grid_data_type = double;

    using container_selector = S;
    static constexpr std::size_t container_ndims = 1;

    static constexpr uint8_t n_neighbors_max = 8;
    using neighbors_cache_type = neighbors_no_cache<n_neighbors_max>;
  };

  /**
   * @brief 2-dimensional grid on the sphere (Dealii).
   *
   * Fastscapelib grid adapter for a Dealii (Hierarchical Equal Area
   * isoLatitude Pixelation of a sphere) grid.
   *
   * @tparam S The container selector for data array members.
   * @tparam T The integer type used to store the Dealii grid node indices.
   */
  template <class S = xt_selector, class T = int>
  class dealii_grid : public grid<dealii_grid<S, T>>
  {
    public:
      using self_type = dealii_grid<S, T>;
      using base_type = grid<self_type>;
      using inner_types = grid_inner_types<self_type>;

      using grid_data_type = typename base_type::grid_data_type;

      using container_selector = typename base_type::container_selector;
      using container_type = xt_tensor_t<container_selector,
            grid_data_type,
            inner_types::container_ndims>;

      using neighbors_type = typename base_type::neighbors_type;
      using neighbors_indices_type = typename base_type::neighbors_indices_type;
      using neighbors_distances_type = typename base_type::neighbors_distances_type;

      using nodes_status_type = typename base_type::nodes_status_type;
      using nodes_status_array_type = xt_tensor_t<container_selector, node_status, 1>;

      using size_type = typename base_type::size_type;
      using shape_type = typename base_type::shape_type;

      dealii_grid(T &triangulation,
                  const nodes_status_array_type &nodes_status,
                  double radius = 6378137.0);

      void set_nodes_status(const nodes_status_array_type &nodes_status);

      std::pair<container_type, container_type> nodes_lonlat() const;
      std::pair<double, double> nodes_lonlat(const size_type &idx) const;
      std::tuple<container_type, container_type, container_type> nodes_xyz() const;
      std::tuple<double, double, double> nodes_xyz(const size_type &idx) const;

      T nside() const;
      double radius() const;

    protected:
      T &m_triangulation;

      shape_type m_shape;
      size_type m_size;
      double m_radius;
      double m_node_area;

      nodes_status_type m_nodes_status;

      using neighbors_distances_impl_type = typename base_type::neighbors_distances_impl_type;
      using neighbors_indices_impl_type = typename base_type::neighbors_indices_impl_type;

      std::vector<size_type> m_neighbors_count;
      std::vector<neighbors_indices_impl_type> m_neighbors_indices;
      std::vector<neighbors_distances_impl_type> m_neighbors_distances;

      void set_neighbors();

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
   * @name Constructors
   */
  /**
   * Creates a new Dealii grid adapter for fastscapelib
   *
   * @param triangulation The Dealii triangulation object.
   * @param radius The radius of the sphere (default: Earth radius in meters).
   */
  template <class S, class T>
  dealii_grid<S, T>::dealii_grid(T &triangulation,
                                 const nodes_status_array_type &nodes_status,
                                 double radius)
    : base_type(0)
    , m_triangulation(triangulation)
    , m_radius(radius)
  {

    m_size = triangulation.n_global_active_cells();
    m_shape = { static_cast<typename shape_type::value_type>(m_size) };
    m_node_area = 4.0 * M_PI * m_radius * m_radius / m_size;

    // will also compute grid node neighbors
    set_nodes_status(nodes_status);
  }
  //@}

  template <class S, class T>
  void dealii_grid<S, T>::set_nodes_status(const nodes_status_array_type &nodes_status)
  {
    if (!xt::same_shape(nodes_status.shape(), m_shape))
      {
        throw std::invalid_argument(
          "invalid shape for nodes_status array (expects shape [N] where N is the total number of nodes)");
      }
    m_nodes_status = nodes_status;

    // maybe invalidates the grid node neighbors so it must be (re)computed
    set_neighbors();
  }

  template <class S, class T>
  void dealii_grid<S, T>::set_neighbors()
  {
    m_neighbors_count.resize(m_size);
    m_neighbors_indices.resize(m_size);
    m_neighbors_distances.resize(m_size);
    unsigned int counter = 0;
    std::unordered_map<unsigned int, unsigned int> global_to_local_index;
    for (const auto &cell : m_triangulation->get_dof_handler().active_cell_iterators())
      {
        m_neighbors_count[counter] = cell->n_faces();
        auto center = cell.center();
        for (auto face_i = 0; face_i < m_neighbors_count[counter]; ++face_i)
          {
            unsigned int global_neighbor_index = 0;
            for (const auto &neighbor_cell :cell.neighbor(face_i))
              {
                global_neighbor_index = neighbor_cell.global_active_cell_index();
                m_neighbors_distances[counter][face_i] = center.distance(neighbor_cell.center());
                // TODO: I assume we only have one neighbor per face, since we start with a uniform grid, we will only have one, but cwe could generatlize this
              }
            // TODO: add assert to check if neighbor_index is correclty set
            global_to_local_index.emplace(global_neighbor_index, counter);
          }

        counter++;
      }
    counter = 0;
    for (const auto &cell : this->get_dof_handler().active_cell_iterators())
      {
        for (auto k = 0; k <= m_neighbors_count[counter]; k++)
          {
            auto global_nb_index = m_neighbors_indices[counter][k];
            m_neighbors_indices[counter][k] = global_to_local_index.at(global_nb_index);
          }
        counter++;
      }
  }

  template <class S, class T>
  inline auto dealii_grid<S, T>::nodes_areas_impl() const -> container_type
  {
    return xt::broadcast(m_node_area, m_shape);
  }

  template <class S, class T>
  inline auto dealii_grid<S, T>::nodes_areas_impl(const size_type & /*idx*/) const noexcept
  -> grid_data_type
  {
    return m_node_area;
  }

  template <class S, class T>
  inline auto dealii_grid<S, T>::neighbors_count_impl(const size_type &idx) const -> size_type
  {
    return m_neighbors_count[idx];
  }

  template <class S, class T>
  void dealii_grid<S, T>::neighbors_indices_impl(neighbors_indices_impl_type &neighbors,
                                                 const size_type &idx) const
  {
    const auto &size = m_neighbors_count[idx];

    for (size_type i = 0; i < size; i++)
      {
        neighbors[i] = m_neighbors_indices[idx][i];
      }
  }

  template <class S, class T>
  auto dealii_grid<S, T>::neighbors_distances_impl(const size_type &idx) const
  -> const neighbors_distances_impl_type &
  {
    return m_neighbors_distances[idx];
  }

  template <class S, class T>
  constexpr std::size_t dealii_grid<S, T>::dimension_impl() noexcept
  {
    return 2;
  }
}

#endif  // FASTSCAPELIB_DEALII_GRID_H__
