#ifndef FASTSCAPELIB_GRID_TRIMESH_H
#define FASTSCAPELIB_GRID_TRIMESH_H

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "xtensor/xstrided_view.hpp"
#include "xtensor/xhistogram.hpp"

#include "fastscapelib/grid/base.hpp"
#include "fastscapelib/utils/xtensor_utils.hpp"


namespace fastscapelib
{

    namespace detail
    {

        /**
         * Used to extract unique edges in the mesh (or count their occurence).
         *
         * This hash function yields the same output for simple permutations,
         * which is what we want since the pairs of node indices (2, 5) and
         * (5, 2) both refer to the same edge.
         */
        template <class T>
        struct tri_edge_hash
        {
            std::size_t operator()(const std::pair<T, T>& p) const
            {
                auto h1 = std::hash<T>()(p.first);
                auto h2 = std::hash<T>()(p.second);

                return h1 ^ h2;
            }
        };

        /**
         * Used to extract unique edges in the mesh (or count their occurence).
         *
         * This comparator returns true for simple permuations. It is needed in
         * addition to ``edge_hash`` in order to build a map of unique edges
         * ignoring permutations.
         */
        template <class T>
        struct tri_edge_equal
        {
            using pair_type = std::pair<T, T>;

            bool operator()(const pair_type& p1, const pair_type& p2) const
            {
                if (p1.first == p2.second && p1.second == p2.first)
                {
                    return true;
                }
                else
                {
                    return p1.first == p2.first && p1.second == p2.second;
                }
            }
        };

        template <class T>
        using tri_edge_type = std::pair<T, T>;

        template <class T>
        using tri_edge_map = std::
            unordered_map<tri_edge_type<T>, std::size_t, tri_edge_hash<T>, tri_edge_equal<T>>;
    }


    template <class S, unsigned int N>
    class trimesh_xt;

    /**
     * 2-d triangular mesh specialized types.
     */
    template <class S, unsigned int N>
    struct grid_inner_types<trimesh_xt<S, N>>
    {
        static constexpr bool is_structured = false;
        static constexpr bool is_uniform = false;

        using grid_data_type = double;

        using xt_selector = S;
        static constexpr std::size_t xt_ndims = 1;

        static constexpr uint8_t n_neighbors_max = N;
        using neighbors_cache_type = neighbors_no_cache<0>;
    };

    /**
     * @brief 2-dimensional triangular (unstructured) mesh.
     *
     * Fastscapelib grid adapter for a 2-d triangular mesh. This class
     * requires an input mesh (it doesn't provide any meshing capability).
     *
     * @tparam S The xtensor container selector for data array members.
     * @tparam N The maximum number of node neighbors.
     */
    template <class S, unsigned int N = 20>
    class trimesh_xt : public grid<trimesh_xt<S, N>>
    {
    public:
        using self_type = trimesh_xt<S, N>;
        using base_type = grid<self_type>;

        using grid_data_type = typename base_type::grid_data_type;

        using xt_selector = typename base_type::xt_selector;
        using size_type = typename base_type::size_type;
        using shape_type = typename base_type::shape_type;

        using points_type = xt_tensor_t<xt_selector, grid_data_type, 2>;
        using triangles_type = xt_tensor_t<xt_selector, size_type, 2>;
        using indices_type = xt_tensor_t<xt_selector, size_type, 1>;
        using areas_type = xt_tensor_t<xt_selector, grid_data_type, 1>;

        using neighbors_type = typename base_type::neighbors_type;
        using neighbors_indices_type = typename base_type::neighbors_indices_type;
        using neighbors_distances_type = typename base_type::neighbors_distances_type;

        using nodes_status_type = typename base_type::nodes_status_type;
        using nodes_status_map_type = typename std::map<size_type, node_status>;
        using nodes_status_array_type = xt_tensor_t<xt_selector, node_status, 1>;

        trimesh_xt(const points_type& points,
                   const triangles_type& triangles,
                   const nodes_status_map_type& nodes_status = {});

        trimesh_xt(const points_type& points,
                   const triangles_type& triangles,
                   const nodes_status_array_type& nodes_status);

    protected:
        using neighbors_distances_impl_type = typename base_type::neighbors_distances_impl_type;
        using neighbors_indices_impl_type = typename base_type::neighbors_indices_impl_type;

        shape_type m_shape;
        size_type m_size;

        points_type m_nodes_points;
        areas_type m_nodes_areas;
        std::unordered_set<size_type> m_boundary_nodes;
        nodes_status_type m_nodes_status;

        std::vector<neighbors_indices_impl_type> m_neighbors_indices;
        std::vector<neighbors_distances_impl_type> m_neighbors_distances;

        void set_size_shape(const points_type& points, const triangles_type& triangles);
        void set_neighbors(const points_type& points, const triangles_type& triangles);
        void set_nodes_status(const nodes_status_map_type& nodes_status);
        void set_nodes_status(const nodes_status_array_type& nodes_status);
        void set_nodes_areas(const points_type& points, const triangles_type& triangles);

        inline areas_type nodes_areas_impl() const;
        inline grid_data_type nodes_areas_impl(const size_type& idx) const noexcept;

        inline size_type neighbors_count_impl(const size_type& idx) const;

        void neighbors_indices_impl(neighbors_indices_impl_type& neighbors,
                                    const size_type& idx) const;

        inline const neighbors_distances_impl_type& neighbors_distances_impl(
            const size_type& idx) const;

        friend class grid<self_type>;
    };


    /**
     * @name Constructors
     */
    //@{
    /**
     * Creates a new triangular mesh.
     *
     * @param points The mesh node x,y coordinates (array of shape [N, 2]).
     * @param triangles The node indices of the triangles (array of shape [K, 3]).
     * @param nodes_status Manually define the status at any node on the mesh.
     *
     * If ``nodes_status`` is empty, a "fixed value" status is set for all
     * boundary nodes (i.e., the end-points of all the edges that are not shared
     * by more than one triangle).
     */
    template <class S, unsigned int N>
    trimesh_xt<S, N>::trimesh_xt(const points_type& points,
                                 const triangles_type& triangles,
                                 const nodes_status_map_type& nodes_status)
        : base_type(0)
        , m_nodes_points(points)
    {
        set_size_shape(points, triangles);
        set_neighbors(points, triangles);
        set_nodes_status(nodes_status);
        set_nodes_areas(points, triangles);
    }

    /**
     * Creates a new triangular mesh.
     *
     * @param points The mesh node x,y coordinates (array of shape [N, 2]).
     * @param triangles The node indices of the triangles (array of shape [K, 3]).
     * @param nodes_status The status of all nodes on the mesh (array of shape [N]).
     */
    template <class S, unsigned int N>
    trimesh_xt<S, N>::trimesh_xt(const points_type& points,
                                 const triangles_type& triangles,
                                 const nodes_status_array_type& nodes_status)
        : base_type(0)
        , m_nodes_points(points)
    {
        set_size_shape(points, triangles);
        set_neighbors(points, triangles);
        set_nodes_status(nodes_status);
        set_nodes_areas(points, triangles);
    }
    //@}

    template <class S, unsigned int N>
    void trimesh_xt<S, N>::set_size_shape(const points_type& points,
                                          const triangles_type& triangles)
    {
        if (points.shape()[1] != 2)
        {
            throw std::invalid_argument("invalid shape for points array (expects shape [N, 2])");
        }
        if (triangles.shape()[1] != 3)
        {
            throw std::invalid_argument("invalid shape for triangles array (expects shape [K, 2])");
        }

        m_size = points.shape()[0];
        m_shape = { static_cast<typename shape_type::value_type>(m_size) };
    }

    template <class S, unsigned int N>
    void trimesh_xt<S, N>::set_neighbors(const points_type& points, const triangles_type& triangles)
    {
        // extract and count triangle edges

        using edge_type = detail::tri_edge_type<size_type>;
        using edge_map = detail::tri_edge_map<size_type>;

        edge_map edges_count;
        const std::array<std::array<size_type, 2>, 3> tri_local_indices{
            { { 1, 2 }, { 2, 0 }, { 0, 1 } }
        };

        size_type n_triangles = triangles.shape()[0];

        for (size_type i = 0; i < n_triangles; i++)
        {
            for (const auto& edge_idx : tri_local_indices)
            {
                const edge_type key(triangles(i, edge_idx[0]), triangles(i, edge_idx[1]));

                auto result = edges_count.insert({ key, 1 });
                // increment edge count if already inserted
                if (!result.second)
                {
                    result.first->second += 1;
                }
            }
        }

        // fill node neighbor data and find boundary nodes

        m_boundary_nodes.clear();
        m_neighbors_indices.resize(m_size);
        m_neighbors_distances.resize(m_size);

        for (const auto& edge : edges_count)
        {
            const edge_type& edge_points = edge.first;
            size_type count = edge.second;

            if (count == 1)
            {
                m_boundary_nodes.insert(edge_points.first);
                m_boundary_nodes.insert(edge_points.second);
            }

            m_neighbors_indices[edge_points.first].push_back(edge_points.second);
            m_neighbors_indices[edge_points.second].push_back(edge_points.first);

            const auto x1 = points(edge_points.first, 0);
            const auto y1 = points(edge_points.first, 1);
            const auto x2 = points(edge_points.second, 0);
            const auto y2 = points(edge_points.second, 1);
            auto distance = std::sqrt(((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)));
            m_neighbors_distances[edge_points.first].push_back(distance);
            m_neighbors_distances[edge_points.second].push_back(distance);
        }
    }

    template <class S, unsigned int N>
    void trimesh_xt<S, N>::set_nodes_areas(const points_type& points,
                                           const triangles_type& triangles)
    {
        size_type n_points = points.shape()[0];
        size_type n_triangles = triangles.shape()[0];
        double just_above_zero = std::numeric_limits<double>::min();

        std::array<std::array<size_type, 2>, 3> local_idx{ { { 1, 2 }, { 2, 0 }, { 0, 1 } } };

        std::array<size_type, 3> coords_shape{ { 3, n_triangles, 2 } };
        xt::xtensor<double, 3> half_edge_coords = xt::empty<double>(coords_shape);

        for (size_type t = 0; t < n_triangles; t++)
        {
            for (size_type i = 0; i < 3; i++)
            {
                auto v1 = local_idx[i][0];
                auto v2 = local_idx[i][1];

                for (size_type j = 0; j < 2; j++)
                {
                    auto p1 = triangles(t, v1);
                    auto p2 = triangles(t, v2);
                    half_edge_coords(i, t, j) = points(p2, j) - points(p1, j);
                }
            }
        }

        xt::xtensor<double, 2> ei_dot_ei = xt::sum(half_edge_coords * half_edge_coords, 2);
        xt::xtensor<double, 2> ei_dot_ej = ei_dot_ei - xt::sum(ei_dot_ei, 0) / 2.0;

        xt::xtensor<double, 1> triangles_areas;
        triangles_areas.resize({ n_triangles });

        for (size_type t = 0; t < n_triangles; t++)
        {
            double area_square
                = 0.25
                  * (ei_dot_ej(2, t) * ei_dot_ej(0, t) + ei_dot_ej(0, t) * ei_dot_ej(1, t)
                     + ei_dot_ej(1, t) * ei_dot_ej(2, t));

            // prevent negative values due to round-off errors
            triangles_areas(t) = std::sqrt(std::max(area_square, just_above_zero));
        }

        auto ce_ratios = -ei_dot_ej * 0.25 / triangles_areas;
        xt::xtensor<double, 2> tri_partitions = ei_dot_ei / 2 * ce_ratios / (3 - 1);

        xt::xtensor<double, 1> weights;
        weights.resize({ n_triangles * 3 });

        for (size_type t = 0; t < n_triangles; t++)
        {
            weights(t) = tri_partitions(1, t) + tri_partitions(2, t);
            weights(n_triangles + t) = tri_partitions(2, t) + tri_partitions(0, t);
            weights(n_triangles * 2 + t) = tri_partitions(0, t) + tri_partitions(1, t);
        }

        auto triangles_t_flat = xt::flatten(xt::transpose(triangles));
        m_nodes_areas = xt::bincount(triangles_t_flat, weights, n_points);

        // avoid area = 0 for isolated nodes (not nice for logarithm)
        for (size_type i = 0; i < n_points; i++)
        {
            if (m_nodes_areas(i) == 0 && neighbors_count_impl(i) == 0)
            {
                m_nodes_areas(i) = just_above_zero;
            }
        }
    }

    template <class S, unsigned int N>
    void trimesh_xt<S, N>::set_nodes_status(const nodes_status_map_type& nodes_status)
    {
        nodes_status_type temp_nodes_status(m_shape, node_status::core);

        if (nodes_status.size() > 0)
        {
            for (const auto& [idx, status] : nodes_status)
            {
                if (status == node_status::looped)
                {
                    throw std::invalid_argument("node_status::looped is not allowed in "
                                                "triangular meshes");
                }

                temp_nodes_status.at(idx) = status;
            }
        }
        else
        {
            // if no status at node is given, set fixed value boundaries for all boundary nodes
            for (const size_type& idx : m_boundary_nodes)
            {
                temp_nodes_status[idx] = node_status::fixed_value;
            }
        }

        m_nodes_status = temp_nodes_status;
    }

    template <class S, unsigned int N>
    void trimesh_xt<S, N>::set_nodes_status(const nodes_status_array_type& nodes_status)
    {
        if (!xt::same_shape(nodes_status.shape(), m_shape))
        {
            throw std::invalid_argument(
                "invalid shape for nodes_status array (expects shape [N] where N is the total number of nodes)");
        }
        m_nodes_status = nodes_status;
    }

    template <class S, unsigned int N>
    inline auto trimesh_xt<S, N>::neighbors_count_impl(const size_type& idx) const -> size_type
    {
        return m_neighbors_indices[idx].size();
    }

    template <class S, unsigned int N>
    inline auto trimesh_xt<S, N>::nodes_areas_impl() const -> areas_type
    {
        return m_nodes_areas;
    }

    template <class S, unsigned int N>
    inline auto trimesh_xt<S, N>::nodes_areas_impl(const size_type& idx) const noexcept
        -> grid_data_type
    {
        return m_nodes_areas(idx);
    }

    template <class S, unsigned int N>
    void trimesh_xt<S, N>::neighbors_indices_impl(neighbors_indices_impl_type& neighbors,
                                                  const size_type& idx) const
    {
        const auto& size = m_neighbors_indices[idx].size();
        neighbors.resize(size);

        for (size_type i = 0; i < size; i++)
        {
            neighbors[i] = m_neighbors_indices[idx][i];
        }
    }

    template <class S, unsigned int N>
    auto trimesh_xt<S, N>::neighbors_distances_impl(const size_type& idx) const
        -> const neighbors_distances_impl_type&
    {
        return m_neighbors_distances[idx];
    }


    /**
     * @typedef trimesh
     *
     * \rst
     * Alias template on ``trimesh_xt`` with :cpp:type:`xt::xtensor`
     * used as array container type for data members.
     *
     * This is mainly for convenience when using in C++ applications.
     * \endrst
     */
    using trimesh = trimesh_xt<xt_selector>;
}

#endif  // FASTSCAPELIB_GRID_TRIMESH_H
