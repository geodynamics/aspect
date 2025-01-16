#ifndef FASTSCAPELIB_DEALII_GRID_H_
#define FASTSCAPELIB_DEALII_GRID_H_

#include "dealii_cxx/dealii_base.h"
#include "dealii_cxx/dealii_tables.h"
#include "dealii_cxx/vec3.h"
#include <unordered_map>

// conflict between dealii xcomplex macro and xtl xcomplex
#undef xcomplex
#include <xtensor/xbroadcast.hpp>
#include <xtensor/xmath.hpp>

#include "fastscapelib/grid/base.hpp"
#include "fastscapelib/utils/consts.hpp"
#include "fastscapelib/utils/xtensor_containers.hpp"

#include <math.h>
#include <memory>
#include <tuple>
#include <utility>
#include <vector>


namespace fastscapelib
{
  namespace detail
  {
    inline double vec3_distance(vec3 a, vec3 b)
    {
      return std::sqrt(std::pow(a.x - b.x, 2) + std::pow(a.y - b.y, 2)
                       + std::pow(a.z - b.z, 2));
    }

  }

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
      using container_type = fixed_shape_container_t<container_selector,
            grid_data_type,
            inner_types::container_ndims>;

      using neighbors_type = typename base_type::neighbors_type;
      using neighbors_indices_type = typename base_type::neighbors_indices_type;
      using neighbors_distances_type = typename base_type::neighbors_distances_type;

      using nodes_status_type = typename base_type::nodes_status_type;
      using nodes_status_array_type = fixed_shape_container_t<container_selector, node_status, 1>;

      using size_type = typename base_type::size_type;
      using shape_type = typename base_type::shape_type;

      dealii_grid(T nside,
                  const nodes_status_array_type &nodes_status,
                  double radius = numeric_constants<>::EARTH_RADIUS_METERS);

      // TODO: factory calculating nside from a given approx. cell area.

      void set_nodes_status(const nodes_status_array_type &nodes_status);

      std::pair<container_type, container_type> nodes_lonlat() const;
      std::pair<double, double> nodes_lonlat(const size_type &idx) const;
      std::tuple<container_type, container_type, container_type> nodes_xyz() const;
      std::tuple<double, double, double> nodes_xyz(const size_type &idx) const;

      T nside() const;
      double radius() const;

    protected:
      using dealii_type = T_Dealii_Base<T>;
      std::unique_ptr<dealii_type> m_dealii_obj_ptr;

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
   * Creates a new Dealii grid
   *
   * @param nside The number of divisions along the side of a base-resolution Dealii pixel.
   * @param radius The radius of the sphere (default: Earth radius in meters).
   */
  template <class S, class T>
  dealii_grid<S, T>::dealii_grid(T &triangulation,
                                 const nodes_status_array_type &nodes_status,
                                 double radius)
    : base_type(0).
    m_triangulation(triangulation),
    , m_radius(radius)
  {

    m_size = triangulation.n_global_active_cells()
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

  /**
   * @name Dealii specific methods
   */
  /**
   * Returns the longitude and latitude of a given grid node (Dealii cell centroid), in radians.
   *
   * @param idx The grid node indice.
   */
  template <class S, class T>
  std::pair<double, double> dealii_grid<S, T>::nodes_lonlat(const size_type &idx) const
  {
    auto ang = m_dealii_obj_ptr->pix2ang(static_cast<int>(idx));

    return std::make_pair<double, double>(double(ang.phi),
                                          xt::numeric_constants<double>::PI_2 - ang.theta);
  }

  /**
   * Returns the longitude and latitude of all grid nodes (Dealii cell centroids), in radians.
   */
  template <class S, class T>
  auto dealii_grid<S, T>::nodes_lonlat() const -> std::pair<container_type, container_type>
  {
    container_type lon(m_shape);
    container_type lat(m_shape);

    for (size_type idx = 0; idx < m_size; ++idx)
      {
        auto ang = m_dealii_obj_ptr->pix2ang(static_cast<int>(idx));
        lon(idx) = ang.phi;
        lat(idx) = xt::numeric_constants<double>::PI_2 - ang.theta;
      }

    return std::make_pair<container_type, container_type>(std::move(lon), std::move(lat));
  }

  /**
   * Returns the x, y, z cartesian coordinates of a given grid node (Dealii cell centroid).
   *
   * @param idx The grid node indice.
   */
  template <class S, class T>
  std::tuple<double, double, double> dealii_grid<S, T>::nodes_xyz(const size_type &idx) const
  {
    auto vec = m_dealii_obj_ptr->pix2vec(static_cast<int>(idx));

    return std::make_tuple<double, double, double>(
             vec.x * m_radius, vec.y * m_radius, vec.z * m_radius);
  }

  /**
   * Returns the x, y, z cartesian coordinates of all grid nodes (Dealii cell centroids).
   */
  template <class S, class T>
  auto dealii_grid<S, T>::nodes_xyz() const
  -> std::tuple<container_type, container_type, container_type>
  {
    container_type x(m_shape);
    container_type y(m_shape);
    container_type z(m_shape);

    for (size_type idx = 0; idx < m_size; ++idx)
      {
        auto vec = m_dealii_obj_ptr->pix2vec(static_cast<int>(idx));
        x(idx) = vec.x * m_radius;
        y(idx) = vec.y * m_radius;
        z(idx) = vec.z * m_radius;
      }

    return std::make_tuple<container_type, container_type, container_type>(
      std::move(x), std::move(y), std::move(z));
  }
  //@}

  template <class S, class T>
  void dealii_grid<S, T>::set_neighbors()
  {
    m_neighbors_count.resize(m_size);
    m_neighbors_indices.resize(m_size);
    m_neighbors_distances.resize(m_size);
    vector<Point<3> centers;
    centers.resize(m_size);
    unsigned int counter = 0;
    std::unordered_map<unsigned int, unsigned int> m_global_to_local_index;
    for (const auto &cell : this->get_dof_handler().active_cell_iterators())
      {
        // to fill:
        // m_global_to_local_index
        // m_neighbors_count
        // m_neighbors_indices
        // m_neighbors_distances
        // m_node_area

        m_neighbors_count[counter] = cell->n_faces();
        Point<3> center = cell.center();
        for (face_i = 0; face_i < m_neighbors_count[counter]; ++face_i)
          {
            unsigned int global_neighbor_index = 0;
            for (const auto &neighbor_cell :cell.neighbor(face_i))
              {
                global_neighbor_index = neighbor_cell.global_active_cell_index();
                m_neighbors_distances[counter][face_i] = center.distance(neighbor_cell.center());
                // TODO: I assume we only have one neighbor per face, since we start with a uniform grid, we will only have one, but cwe could generatlize this
              }
            // TODO: add assert to check if neighbor_index is correclty set
            m_global_to_local_index.emplace(std::make_pair(global_neighbor_index,counter));

          }

        counter++;
      }
    counter = 0;
    for (const auto &cell : this->get_dof_handler().active_cell_iterators())
      {
        for (const auto &neighbor_cell :cell.neighbor(face_i))
          {
            m_neighbors_indices[counter][face_i] = global_to_local.at(neighbor_cell.global_active_cell_index());
          }
        counter++;
      }
  }

  /**
   * @name Dealii specific properties
   */

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
