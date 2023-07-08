/**
 * A profile grid is a one dimensional grid.
 */
#ifndef FASTSCAPELIB_GRID_PROFILE_GRID_H
#define FASTSCAPELIB_GRID_PROFILE_GRID_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <map>
#include <vector>

#include <xtensor/xbroadcast.hpp>

#include "fastscapelib/grid/base.hpp"
#include "fastscapelib/grid/structured_grid.hpp"


namespace fastscapelib
{

    /**
     * Status at grid boundary nodes.
     */
    class profile_boundary_status : public boundary_status
    {
    public:
        node_status left = node_status::core;  /**< Status at left edge/border node(s) */
        node_status right = node_status::core; /**< Status at right edge/border node(s) */

        profile_boundary_status(const node_status status);
        profile_boundary_status(const node_status left_status, const node_status right_status);
        profile_boundary_status(const std::array<node_status, 2>& status);

        bool is_horizontal_looped() const;

    protected:
        void check_looped_symmetrical();
    };

    /**
     * @name Constructors
     */
    //@{
    /**
     * Set the same status for all boundary nodes.
     */
    inline profile_boundary_status::profile_boundary_status(node_status status)
        : left(status)
        , right(status)
    {
        check_looped_symmetrical();
    }

    /**
     * Set status at the left and right edge nodes of a profile grid.
     */
    inline profile_boundary_status::profile_boundary_status(node_status left_status,
                                                            node_status right_status)
        : left(left_status)
        , right(right_status)
    {
        check_looped_symmetrical();
    }

    /**
     * Set status at the left and right edge nodes of a profile grid.
     */
    inline profile_boundary_status::profile_boundary_status(
        const std::array<node_status, 2>& status)
        : left(status[0])
        , right(status[1])
    {
        check_looped_symmetrical();
    }
    //@}

    /**
     * Return true if periodic (looped) conditions are set for ``left`` and ``right``.
     */
    inline bool profile_boundary_status::is_horizontal_looped() const
    {
        return is_looped(left) && is_looped(right);
    }

    inline void profile_boundary_status::check_looped_symmetrical()
    {
        if (is_looped(left) ^ is_looped(right))
        {
            throw std::invalid_argument("looped boundaries are not symmetrical");
        }
    }

    //*******************
    //* Profile grid (1D)
    //*******************

    template <class S, class C>
    class profile_grid_xt;

    /**
     * Profile grid specialized types.
     */
    template <class S, class C>
    struct grid_inner_types<profile_grid_xt<S, C>>
    {
        static constexpr bool is_structured = true;
        static constexpr bool is_uniform = true;

        using grid_data_type = double;

        using xt_selector = S;
        static constexpr std::size_t xt_ndims = 1;

        static constexpr uint8_t n_neighbors_max = 2u;
        using neighbors_cache_type = C;

        // structured_grid types
        using length_type = double;
        using spacing_type = double;
    };

    /**
     * 1-dimensional uniform grid.
     *
     * Useful for modeling single channel or hillslope profiles.
     *
     * @tparam S The xtensor container selector for data array members.
     * @tparam C The grid neighbor nodes cache type.
     */
    template <class S, class C = neighbors_cache<2>>
    class profile_grid_xt : public structured_grid<profile_grid_xt<S, C>>
    {
    public:
        using self_type = profile_grid_xt<S, C>;
        using base_type = structured_grid<self_type>;
        using inner_types = grid_inner_types<self_type>;

        using grid_data_type = typename base_type::grid_data_type;

        using xt_selector = typename base_type::xt_selector;
        using xt_type = xt_tensor_t<xt_selector, grid_data_type, inner_types::xt_ndims>;

        using size_type = typename base_type::size_type;
        using shape_type = typename base_type::shape_type;

        using length_type = typename base_type::length_type;
        using spacing_type = typename base_type::spacing_type;

        using code_type = std::uint8_t;

        using neighbors_type = typename base_type::neighbors_type;
        using neighbors_indices_type = typename base_type::neighbors_indices_type;
        using neighbors_distances_type = typename base_type::neighbors_distances_type;

        using node_status_type = typename base_type::node_status_type;

        profile_grid_xt(size_type size,
                        spacing_type spacing,
                        const profile_boundary_status& bounds_status,
                        const std::vector<node>& nodes_status = {});

        static profile_grid_xt from_length(size_type size,
                                           length_type length,
                                           const profile_boundary_status& bounds_status,
                                           const std::vector<node>& nodes_status = {});

    protected:
        using neighbors_distances_impl_type = typename base_type::neighbors_distances_impl_type;
        using neighbors_indices_impl_type = typename base_type::neighbors_indices_impl_type;
        using coded_ndistances_type = std::array<neighbors_distances_impl_type, 3>;

        shape_type m_shape;
        size_type m_size;
        spacing_type m_spacing;
        length_type m_length;
        grid_data_type m_node_area;

        xt::xtensor<code_type, 1> m_gcode_idx;

        node_status_type m_nodes_status;
        profile_boundary_status m_bounds_status;

        void set_nodes_status(const std::vector<node>& nodes_status);

        static constexpr std::array<std::ptrdiff_t, 3> offsets{ { 0, -1, 1 } };
        std::vector<neighbors_type> m_all_neighbors;

        coded_ndistances_type m_neighbors_distances;
        void build_neighbors_distances();

        void build_gcode();
        code_type gcode(const size_type& idx) const;

        std::array<size_type, 3> m_neighbors_count;
        void build_neighbors_count();

        inline xt_type nodes_areas_impl() const;
        inline grid_data_type nodes_areas_impl(const size_type& idx) const noexcept;

        inline size_type neighbors_count_impl(const size_type& idx) const noexcept;

        void neighbors_indices_impl(neighbors_indices_impl_type& neighbors,
                                    const size_type& idx) const;

        const neighbors_distances_impl_type& neighbors_distances_impl(const size_type& idx) const;

        friend class structured_grid<self_type>;
        friend class grid<self_type>;
    };

    template <class S, class C>
    constexpr std::array<std::ptrdiff_t, 3> profile_grid_xt<S, C>::offsets;

    /**
     * @name Constructors
     */
    //@{
    /**
     * Creates a new profile grid.
     *
     * @param size Total number of grid nodes.
     * @param spacing Distance between two adjacent grid nodes.
     * @param bounds_status Status at boundary nodes (left/right grid edges).
     * @param nodes_status Manually define the status at any node on the grid.
     */
    template <class S, class C>
    profile_grid_xt<S, C>::profile_grid_xt(size_type size,
                                           spacing_type spacing,
                                           const profile_boundary_status& bounds_status,
                                           const std::vector<node>& nodes_status)
        : base_type(size)
        , m_size(size)
        , m_spacing(spacing)
        , m_bounds_status(bounds_status)
    {
        m_shape = { static_cast<typename shape_type::value_type>(m_size) };
        m_length = static_cast<spacing_type>(size - 1) * spacing;
        m_node_area = spacing;
        set_nodes_status(nodes_status);
        build_gcode();
        build_neighbors_distances();
        build_neighbors_count();
    }
    //@}

    /**
     * @name Factories
     */
    //@{
    /**
     * Creates a new profile grid.
     *
     * @param size Total number of grid nodes.
     * @param length Total physical length of the grid.
     * @param bounds_status Status at boundary nodes (left & right grid edges).
     * @param nodes_status Manually define the status at any node on the grid.
     */
    template <class S, class C>
    profile_grid_xt<S, C> profile_grid_xt<S, C>::from_length(
        size_type size,
        length_type length,
        const profile_boundary_status& bounds_status,
        const std::vector<node>& nodes_status)
    {
        spacing_type spacing = length / static_cast<length_type>(size - 1);
        return profile_grid_xt<S, C>(size, spacing, bounds_status, nodes_status);
    }
    //@}

    template <class S, class C>
    void profile_grid_xt<S, C>::build_gcode()
    {
        m_gcode_idx.resize({ m_size });

        m_gcode_idx.fill(1);
        m_gcode_idx[0] = 0;
        m_gcode_idx[m_size - 1] = 2;
    }

    template <class S, class C>
    auto profile_grid_xt<S, C>::gcode(const size_type& idx) const -> code_type
    {
        return m_gcode_idx[idx];
    }

    template <class S, class C>
    void profile_grid_xt<S, C>::build_neighbors_distances()
    {
        m_neighbors_distances.fill({ m_spacing, m_spacing });
    }

    template <class S, class C>
    void profile_grid_xt<S, C>::build_neighbors_count()
    {
        if (m_bounds_status.is_horizontal_looped())
        {
            m_neighbors_count = std::array<size_type, 3>({ 2, 2, 2 });
        }
        else
        {
            m_neighbors_count = std::array<size_type, 3>({ 1, 2, 1 });
        }
    }

    template <class S, class C>
    void profile_grid_xt<S, C>::set_nodes_status(const std::vector<node>& nodes_status)
    {
        node_status_type temp_nodes_status(m_shape, node_status::core);

        temp_nodes_status(0) = m_bounds_status.left;
        temp_nodes_status(m_size - 1) = m_bounds_status.right;

        for (const node& n : nodes_status)
        {
            if (n.status == node_status::looped)
            {
                throw std::invalid_argument("node_status::looped is not allowed in "
                                            "'nodes_status' "
                                            "(use 'bounds_status' instead)");
            }
            else if (temp_nodes_status.at(n.idx) == node_status::looped)
            {
                throw std::invalid_argument("cannot overwrite the status of a "
                                            "looped boundary node");
            }

            temp_nodes_status.at(n.idx) = n.status;
        }

        m_nodes_status = temp_nodes_status;
    }

    template <class S, class C>
    inline auto profile_grid_xt<S, C>::nodes_areas_impl() const -> xt_type
    {
        return xt::broadcast(m_node_area, m_shape);
    }

    template <class S, class C>
    inline auto profile_grid_xt<S, C>::nodes_areas_impl(const size_type& /*idx*/) const noexcept
        -> grid_data_type
    {
        return m_node_area;
    }

    template <class S, class C>
    inline auto profile_grid_xt<S, C>::neighbors_count_impl(const size_type& idx) const noexcept
        -> size_type
    {
        return m_neighbors_count[gcode(idx)];
    }

    template <class S, class C>
    auto profile_grid_xt<S, C>::neighbors_distances_impl(const size_type& /*idx*/) const
        -> const neighbors_distances_impl_type&
    {
        return m_neighbors_distances[0];
    }

    template <class S, class C>
    inline auto profile_grid_xt<S, C>::neighbors_indices_impl(
        neighbors_indices_impl_type& neighbors, const size_type& idx) const -> void
    {
        if (idx == 0)
        {
            if (m_bounds_status.is_horizontal_looped())
            {
                neighbors[0] = m_size - 1;
                neighbors[1] = 1;
            }
            else
            {
                neighbors[0] = 1;
            }
        }
        else if (idx == m_size - 1)
        {
            neighbors[0] = m_size - 2;
            if (m_bounds_status.is_horizontal_looped())
            {
                neighbors[1] = 0;
            }
        }
        else
        {
            for (std::size_t k = 1; k < 3; ++k)
            {
                std::size_t nb_idx = detail::add_offset(idx, offsets[k]);
                neighbors[k - 1] = nb_idx;
            }
        }
    }

    /**
     * @typedef profile_grid
     *
     * \rst
     * Alias template on ``profile_grid_xt`` with :cpp:type:`xt::xtensor`
     * used as array container type for data members.
     *
     * This is mainly for convenience when using in C++ applications.
     * \endrst
     */
    using profile_grid = profile_grid_xt<xt_selector>;
}

#endif
