/**
 * A raster grid represents a 2D uniform grid.
 */
#ifndef FASTSCAPELIB_GRID_RASTER_GRID_H
#define FASTSCAPELIB_GRID_RASTER_GRID_H

#include <map>
#include <vector>
#include <utility>

#include <xtensor/xbroadcast.hpp>
#include <xtensor/xtensor.hpp>

#include "fastscapelib/grid/base.hpp"
#include "fastscapelib/grid/structured_grid.hpp"
#include "fastscapelib/grid/profile_grid.hpp"


namespace fastscapelib
{
    /*****************
     * Grid elements *
     *****************/

    /**
     * Represents a raster grid node.
     */
    struct raster_node
    {
        std::size_t row;    /**< Row index */
        std::size_t col;    /**< Column index */
        node_status status; /**< Node status */
    };

    /**
     * Represents a raster grid node neighbor.
     */
    struct raster_neighbor
    {
        std::size_t flatten_idx; /**< Flat index of the neighbor node */
        std::size_t row;         /**< Row index of the neighbor node */
        std::size_t col;         /**< Column index of the neighbor node */
        double distance;         /**< Distance to the neighbor node */
        node_status status;      /**< Status at the neighbor node */

        bool operator==(const raster_neighbor& rhs) const
        {
            return (flatten_idx == rhs.flatten_idx) && (row == rhs.row) && (col == rhs.col)
                   && (distance == rhs.distance) && (status == rhs.status);
        }
    };

    /**
     * Status at grid boundary nodes.
     *
     * To disambiguate the cases at each of the four raster grid corners, node
     * status is set from one of their two overlaping borders according to the
     * following precedance order:
     *
     * fixed value > fixed gradient > looped > core
     */
    class raster_boundary_status : public boundary_status
    {
    public:
        node_status left = node_status::core;   /**< Status at left border nodes */
        node_status right = node_status::core;  /**< Status at right border nodes */
        node_status top = node_status::core;    /**< Status at top border nodes */
        node_status bottom = node_status::core; /**< Status at bottom border nodes */

        raster_boundary_status(node_status status);
        raster_boundary_status(const std::array<node_status, 4>& status);

        bool is_vertical_looped() const;
        bool is_horizontal_looped() const;

    protected:
        void check_looped_symmetrical() const;
    };

    /**
     * @name Constructors
     */
    //@{
    /**
     * Set the same status for all boundary nodes.
     */
    inline raster_boundary_status::raster_boundary_status(node_status status)
        : left(status)
        , right(status)
        , top(status)
        , bottom(status)
    {
        check_looped_symmetrical();
    }

    /**
     * Set status at the left, right, top and bottom border nodes of a raster grid.
     */
    inline raster_boundary_status::raster_boundary_status(const std::array<node_status, 4>& status)
        : left(status[0])
        , right(status[1])
        , top(status[2])
        , bottom(status[3])
    {
        check_looped_symmetrical();
    }
    //@}

    /**
     * Return true if periodic (looped) conditions are set for ``left`` and ``right``.
     */
    inline bool raster_boundary_status::is_horizontal_looped() const
    {
        return is_looped(left) && is_looped(right);
    }

    /**
     * Return true if periodic (looped) conditions are set for ``top`` and ``bottom``.
     */
    inline bool raster_boundary_status::is_vertical_looped() const
    {
        return is_looped(top) && is_looped(bottom);
    }

    /**
     * Throw an error if the boundary status is not symmetrical.
     */
    inline void raster_boundary_status::check_looped_symmetrical() const
    {
        if (is_looped(left) ^ is_looped(right) || is_looped(top) ^ is_looped(bottom))
        {
            throw std::invalid_argument("looped boundaries are not symmetrical");
        }
    }

    /*******************
     * Raster grid (2D) *
     ********************/

    /**
     * Raster grid node connectivity.
     */
    enum class raster_connect
    {
        rook = 0,  /**< 4-connectivity (vertical or horizontal) */
        queen = 1, /**< 8-connectivity (including diagonals) */
        bishop = 2 /**< 4-connectivity (only diagonals) */
    };

    /*************************
     * raster_neighbors_base *
     *************************/

    struct raster_neighbors_base
    {
        using neighbors_offsets_type = std::vector<std::array<std::ptrdiff_t, 2>>;
        using size_type = std::size_t;
    };

    /********************
     * raster_neighbors *
     ********************/

    template <raster_connect RC>
    struct raster_neighbors
    {
    };


    template <>
    struct raster_neighbors<raster_connect::queen> : public raster_neighbors_base
    {
        static constexpr std::uint8_t _n_neighbors_max = 8u;

        neighbors_offsets_type node_neighbors_offsets(std::ptrdiff_t up,
                                                      std::ptrdiff_t down,
                                                      std::ptrdiff_t left,
                                                      std::ptrdiff_t right) const;

        std::array<size_type, 9> build_neighbors_count(
            const raster_boundary_status& bounds_status) const;
    };

    /**
     * Given the possible moves (up/down/left/right) depending on grid periodicity,
     * return a list of offsets corresponding to "queen" neighbors relative locations.
     * Offsets are returned in row-major layout when accessible (else dropped).
     *
     * Examples:
     *
     *         All          w/o left       w/o top
     *      0   1   2        0   1
     *        \ | /          | /
     *      3 - N - 4        N - 2        0 - N - 1
     *        / | \          | \            / | \
     *      5   6   7        3   4        2   3   4
     */
    inline auto raster_neighbors<raster_connect::queen>::node_neighbors_offsets(
        std::ptrdiff_t up, std::ptrdiff_t down, std::ptrdiff_t left, std::ptrdiff_t right) const
        -> neighbors_offsets_type
    {
        std::array<bool, 8> mask{
            (up != 0 && left != 0),   up != 0,   (up != 0 && right != 0),  left != 0, right != 0,
            (down != 0 && left != 0), down != 0, (down != 0 && right != 0)
        };

        std::array<neighbors_offsets_type::value_type, 8> offsets{ { { up, left },
                                                                     { up, 0 },
                                                                     { up, right },
                                                                     { 0, left },
                                                                     { 0, right },
                                                                     { down, left },
                                                                     { down, 0 },
                                                                     { down, right } } };


        neighbors_offsets_type selected_offsets{};
        selected_offsets.reserve(8);
        for (std::size_t i = 0; i < 8; ++i)
            if (mask[i])
                selected_offsets.push_back(offsets[i]);

        return selected_offsets;
    }

    /**
     * Build the neighbors count of each node code (e.g. position on the grid)
     * depending on the periodicity of the grid.
     */
    inline auto raster_neighbors<raster_connect::queen>::build_neighbors_count(
        const raster_boundary_status& bounds_status) const -> std::array<size_type, 9>
    {
        std::array<size_type, 9> neighbors_count;

        if (bounds_status.is_vertical_looped() && bounds_status.is_horizontal_looped())
        {
            neighbors_count.fill(8);
        }
        else if (bounds_status.is_vertical_looped())
        {
            neighbors_count = std::array<size_type, 9>({ 5, 8, 5, 5, 8, 5, 5, 8, 5 });
        }
        else if (bounds_status.is_horizontal_looped())
        {
            neighbors_count = std::array<size_type, 9>({ 5, 5, 5, 8, 8, 8, 5, 5, 5 });
        }
        else
        {
            neighbors_count = std::array<size_type, 9>({ 3, 5, 3, 5, 8, 5, 3, 5, 3 });
        }

        return neighbors_count;
    }


    template <>
    struct raster_neighbors<raster_connect::rook> : public raster_neighbors_base
    {
        static constexpr std::uint8_t _n_neighbors_max = 4u;

        neighbors_offsets_type node_neighbors_offsets(std::ptrdiff_t up,
                                                      std::ptrdiff_t down,
                                                      std::ptrdiff_t left,
                                                      std::ptrdiff_t right) const;

        std::array<size_type, 9> build_neighbors_count(
            const raster_boundary_status& bounds_status) const;
    };

    /**
     * Given the possible moves (up/down/left/right) depending on grid periodicity,
     * return a list of offsets corresponding to "rook" neighbors relative locations.
     * Offsets are returned in row-major layout when accessible (else dropped).
     *
     * Examples:
     *
     *       All         w/o left        w/o top         etc.
     *        0            0
     *        |            |
     *   1 -- N -- 2       N -- 1      0 -- N -- 1
     *        |            |                |
     *        3            2                2
     */
    inline auto raster_neighbors<raster_connect::rook>::node_neighbors_offsets(
        std::ptrdiff_t up, std::ptrdiff_t down, std::ptrdiff_t left, std::ptrdiff_t right) const
        -> neighbors_offsets_type
    {
        std::array<bool, 4> mask{ up != 0, left != 0, right != 0, down != 0 };

        std::array<neighbors_offsets_type::value_type, 4> offsets{
            { { up, 0 }, { 0, left }, { 0, right }, { down, 0 } }
        };

        neighbors_offsets_type selected_offsets{};
        selected_offsets.reserve(4);
        for (std::size_t i = 0; i < 4; ++i)
            if (mask[i])
                selected_offsets.push_back(offsets[i]);

        return selected_offsets;
    }

    /**
     * Build the neighbors count of each node code (e.g. position on the grid)
     * depending on the periodicity of the grid.
     */
    inline auto raster_neighbors<raster_connect::rook>::build_neighbors_count(
        const raster_boundary_status& bounds_status) const -> std::array<size_type, 9>
    {
        std::array<size_type, 9> neighbors_count;

        if (bounds_status.is_vertical_looped() && bounds_status.is_horizontal_looped())
        {
            neighbors_count.fill(4);
        }
        else if (bounds_status.is_vertical_looped())
        {
            neighbors_count = std::array<size_type, 9>({ 3, 4, 3, 3, 4, 3, 3, 4, 3 });
        }
        else if (bounds_status.is_horizontal_looped())
        {
            neighbors_count = std::array<size_type, 9>({ 3, 3, 3, 4, 4, 4, 3, 3, 3 });
        }
        else
        {
            neighbors_count = std::array<size_type, 9>({ 2, 3, 2, 3, 4, 3, 2, 3, 2 });
        }

        return neighbors_count;
    }


    template <>
    struct raster_neighbors<raster_connect::bishop> : public raster_neighbors_base
    {
        static constexpr std::uint8_t _n_neighbors_max = 4u;

        neighbors_offsets_type node_neighbors_offsets(std::ptrdiff_t up,
                                                      std::ptrdiff_t down,
                                                      std::ptrdiff_t left,
                                                      std::ptrdiff_t right) const;

        std::array<size_type, 9> build_neighbors_count(
            const raster_boundary_status& bounds_status) const;
    };

    /**
     * Given the possible moves (up/down/left/right) depending on grid periodicity,
     * return a list of offsets corresponding to "bishop" neighbors relative locations.
     * Offsets are returned in row-major layout when accessible (else dropped).
     *
     * Examples:
     *
     *         All          w/o left       w/o top
     *      0       1            0
     *        \   /            /
     *          N            N                N
     *        /   \            \            /   \
     *      2       3            1        0       1
     */
    inline auto raster_neighbors<raster_connect::bishop>::node_neighbors_offsets(
        std::ptrdiff_t up, std::ptrdiff_t down, std::ptrdiff_t left, std::ptrdiff_t right) const
        -> neighbors_offsets_type
    {
        std::array<bool, 4> mask{ (up != 0 && left != 0),
                                  (up != 0 && right != 0),
                                  (down != 0 && left != 0),
                                  (down != 0 && right != 0) };

        std::array<neighbors_offsets_type::value_type, 4> offsets{
            { { up, left }, { up, right }, { down, left }, { down, right } }
        };

        neighbors_offsets_type selected_offsets{};
        selected_offsets.reserve(4);
        for (std::size_t i = 0; i < 4; ++i)
            if (mask[i])
                selected_offsets.push_back(offsets[i]);

        return selected_offsets;
    }

    /**
     * Build the neighbors count of each node code (e.g. position on the grid)
     * depending on the periodicity of the grid.
     */
    inline auto raster_neighbors<raster_connect::bishop>::build_neighbors_count(
        const raster_boundary_status& bounds_status) const -> std::array<size_type, 9>
    {
        std::array<size_type, 9> neighbors_count;

        if (bounds_status.is_vertical_looped() && bounds_status.is_horizontal_looped())
        {
            neighbors_count.fill(4);
        }
        else if (bounds_status.is_vertical_looped())
        {
            neighbors_count = std::array<size_type, 9>({ 2, 4, 2, 2, 4, 2, 2, 4, 2 });
        }
        else if (bounds_status.is_horizontal_looped())
        {
            neighbors_count = std::array<size_type, 9>({ 2, 2, 2, 4, 4, 4, 2, 2, 2 });
        }
        else
        {
            neighbors_count = std::array<size_type, 9>({ 1, 2, 1, 2, 4, 2, 1, 2, 1 });
        }

        return neighbors_count;
    }

    /********************
     * grid_inner_types *
     ********************/

    template <class S, raster_connect RC, class C>
    class raster_grid;

    /**
     * Raster grid specialized types.
     */
    template <class S, raster_connect RC, class C>
    struct grid_inner_types<raster_grid<S, RC, C>>
    {
        static constexpr bool is_structured = true;
        static constexpr bool is_uniform = true;

        using grid_data_type = double;

        using container_selector = S;
        static constexpr std::size_t container_ndims = 2;

        static constexpr uint8_t n_neighbors_max = raster_neighbors<RC>::_n_neighbors_max;
        using neighbors_cache_type = C;

        // structured_grid types
        using length_type = std::array<grid_data_type, 2>;
        using spacing_type = std::array<grid_data_type, 2>;
    };

    /**
     * 2-dimensional uniform (raster) grid.
     *
     * @tparam S The container selector for data array members.
     * @tparam RC The kind of raster node connectivity.
     * @tparam C The grid neighbor nodes cache type.
     */
    template <class S = xt_selector,
              raster_connect RC = raster_connect::queen,
              class C = neighbors_cache<raster_neighbors<RC>::_n_neighbors_max>>
    class raster_grid
        : public structured_grid<raster_grid<S, RC, C>>
        , public raster_neighbors<RC>
    {
    public:
        using self_type = raster_grid<S, RC, C>;
        using base_type = structured_grid<self_type>;
        using inner_types = grid_inner_types<self_type>;

        using grid_data_type = typename base_type::grid_data_type;

        using container_selector = typename base_type::container_selector;
        using container_type = fixed_shape_container_t<container_selector,
                                                       grid_data_type,
                                                       inner_types::container_ndims>;

        using size_type = typename base_type::size_type;
        using shape_type = typename base_type::shape_type;

        using length_type = typename base_type::length_type;
        using spacing_type = typename base_type::spacing_type;

        using code_type = std::uint8_t;

        // row, col index pair
        using raster_idx_type = std::pair<size_type, size_type>;

        using neighbors_type = typename base_type::neighbors_type;
        using neighbors_indices_type = typename base_type::neighbors_indices_type;
        using neighbors_distances_type = typename base_type::neighbors_distances_type;
        using neighbors_indices_raster_type = std::vector<raster_idx_type>;
        using neighbors_raster_type = std::vector<raster_neighbor>;

        using nodes_status_type = typename base_type::nodes_status_type;
        using nodes_status_map_type = typename std::map<raster_idx_type, node_status>;

        raster_grid(const shape_type& shape,
                    const spacing_type& spacing,
                    const raster_boundary_status& bounds_status,
                    const nodes_status_map_type& nodes_status = {});

        static raster_grid from_length(const shape_type& shape,
                                       const length_type& length,
                                       const raster_boundary_status& bounds_status,
                                       const nodes_status_map_type& nodes_status = {});

        shape_type shape() const noexcept;
        raster_boundary_status bounds_status() const noexcept;

        using base_type::neighbors_indices;

        inline neighbors_indices_raster_type& neighbors_indices(
            const size_type& row,
            const size_type& col,
            neighbors_indices_raster_type& neighbors_indices);

        inline neighbors_indices_raster_type neighbors_indices(const size_type& row,
                                                               const size_type& col);

        using base_type::neighbors_distances;

        using base_type::neighbors;
        inline neighbors_raster_type& neighbors(const size_type& row,
                                                const size_type& col,
                                                neighbors_raster_type& neighbors);
        inline neighbors_raster_type neighbors(const size_type& row, const size_type& col);

        code_type nodes_codes(const size_type& row, const size_type& col) const noexcept;
        code_type nodes_codes(const size_type& idx) const noexcept;

    private:
        using neighbors_distances_impl_type = typename base_type::neighbors_distances_impl_type;
        using neighbors_indices_impl_type = typename base_type::neighbors_indices_impl_type;
        using neighbors_offsets_type = typename raster_neighbors<RC>::neighbors_offsets_type;

        using coded_ncount_type = std::array<size_type, 9>;
        using coded_noffsets_type = std::array<neighbors_offsets_type, 9>;
        using coded_ndistances_type = std::array<neighbors_distances_impl_type, 9>;
        using nodes_codes_type = std::vector<code_type>;

        shape_type m_shape;
        size_type m_size;
        spacing_type m_spacing;
        length_type m_length;
        grid_data_type m_node_area;

        nodes_status_type m_nodes_status;
        raster_boundary_status m_bounds_status;

        struct corner_node
        {
            size_type row;
            size_type col;
            node_status row_border;
            node_status col_border;
        };

        nodes_codes_type m_nodes_codes;

        coded_ncount_type m_neighbors_count;
        coded_noffsets_type m_neighbor_offsets;
        coded_ndistances_type m_neighbor_distances;

        inline size_type ravel_idx(const size_type& row, const size_type& col) const noexcept;
        inline raster_idx_type unravel_idx(const size_type& idx) const noexcept;

        void set_nodes_status(const nodes_status_map_type& nodes_status);

        void build_nodes_codes();
        coded_noffsets_type build_coded_neighbors_offsets();
        coded_ndistances_type build_coded_neighbors_distances();

        inline const neighbors_offsets_type& neighbor_offsets(code_type code) const noexcept;

        inline container_type nodes_areas_impl() const;
        inline grid_data_type nodes_areas_impl(const size_type& idx) const noexcept;

        inline size_type neighbors_count_impl(const size_type& idx) const noexcept;

        inline const neighbors_distances_impl_type& neighbors_distances_impl(
            const size_type& idx) const noexcept;

        void neighbors_indices_impl(neighbors_indices_impl_type& neighbors,
                                    const size_type& idx) const;

        friend class structured_grid<self_type>;
        friend class grid<self_type>;
    };

    /**
     * @name Constructors
     */
    //@{
    /**
     * Creates a new raster grid.
     *
     * @param shape Shape of the grid (number of rows and cols).
     * @param spacing Distance between two adjacent rows / cols.
     * @param bounds_status Status at boundary nodes (left/right/top/bottom borders).
     * @param nodes_status Manually define the status at any node on the grid.
     */
    template <class S, raster_connect RC, class C>
    raster_grid<S, RC, C>::raster_grid(const shape_type& shape,
                                       const spacing_type& spacing,
                                       const raster_boundary_status& bounds_status,
                                       const nodes_status_map_type& nodes_status)
        : base_type(shape[0] * shape[1])
        , m_shape(shape)
        , m_spacing(spacing)
        , m_bounds_status(bounds_status)
    {
        m_size = shape[0] * shape[1];
        m_length = { (static_cast<double>(shape[0]) - 1) * spacing[0],
                     (static_cast<double>(shape[1]) - 1) * spacing[1] };
        m_node_area = spacing[0] * spacing[1];

        build_nodes_codes();
        m_neighbors_count = this->build_neighbors_count(bounds_status);

        set_nodes_status(nodes_status);

        m_neighbor_offsets = build_coded_neighbors_offsets();
        m_neighbor_distances = build_coded_neighbors_distances();
    }
    //@}

    /**
     * @name Factories
     */
    //@{
    /**
     * Creates a new raster grid.
     *
     * @param shape Shape of the grid (number of rows and cols).
     * @param length Total physical lengths of the grid.
     * @param bounds_status Status at boundary nodes (left/right/top/bottom borders).
     * @param nodes_status Manually define the status at any node on the grid.
     */
    template <class S, raster_connect RC, class C>
    raster_grid<S, RC, C> raster_grid<S, RC, C>::from_length(
        const shape_type& shape,
        const length_type& length,
        const raster_boundary_status& bounds_status,
        const nodes_status_map_type& nodes_status)
    {
        spacing_type spacing = { length[0] / (static_cast<double>(shape[0]) - 1),
                                 length[1] / (static_cast<double>(shape[1]) - 1) };
        return raster_grid<S, RC, C>(shape, spacing, bounds_status, nodes_status);
    }
    //@}

    /**
     * @name Grid properties
     */
    //@{
    /**
     * Returns the shape of the grid node arrays.
     */
    template <class S, raster_connect RC, class C>
    auto raster_grid<S, RC, C>::shape() const noexcept -> shape_type
    {
        return m_shape;
    }

    /**
     * Returns the grid node status at grid left / right / top / bottom borders.
     */
    template <class S, raster_connect RC, class C>
    auto raster_grid<S, RC, C>::bounds_status() const noexcept -> raster_boundary_status
    {
        return m_bounds_status;
    }
    //@}

    template <class S, raster_connect RC, class C>
    inline auto raster_grid<S, RC, C>::ravel_idx(const size_type& row,
                                                 const size_type& col) const noexcept -> size_type
    {
        // TODO: assumes row-major layout -> support col-major?
        return row * m_shape[1] + col;
    }

    template <class S, raster_connect RC, class C>
    inline auto raster_grid<S, RC, C>::unravel_idx(const size_type& idx) const noexcept
        -> raster_idx_type
    {
        // TODO: assumes row-major layout -> support col-major?
        auto ncols = m_shape[1];
        size_type row = idx / ncols;
        size_type col = idx - row * ncols;

        return std::make_pair(row, col);
    }

    template <class S, raster_connect RC, class C>
    void raster_grid<S, RC, C>::set_nodes_status(const nodes_status_map_type& nodes_status)
    {
        nodes_status_type temp_nodes_status
            = container_impl<nodes_status_type>::init(m_shape, node_status::core);
        const auto nrows = static_cast<size_type>(m_shape[0]);
        const auto ncols = static_cast<size_type>(m_shape[1]);

        // set border nodes
        container_impl<container_type>::get_left_view(temp_nodes_status) = m_bounds_status.left;
        container_impl<container_type>::get_right_view(temp_nodes_status) = m_bounds_status.right;
        container_impl<container_type>::get_top_view(temp_nodes_status) = m_bounds_status.top;
        container_impl<container_type>::get_bottom_view(temp_nodes_status) = m_bounds_status.bottom;

        // set corner nodes
        std::vector<corner_node> corners
            = { { 0, 0, m_bounds_status.top, m_bounds_status.left },
                { 0, ncols - 1, m_bounds_status.top, m_bounds_status.right },
                { nrows - 1, 0, m_bounds_status.bottom, m_bounds_status.left },
                { nrows - 1, ncols - 1, m_bounds_status.bottom, m_bounds_status.right } };

        for (const auto& c : corners)
        {
            node_status cs = std::max(c.row_border, c.col_border, detail::node_status_cmp);
            temp_nodes_status(c.row, c.col) = cs;
        }

        // set user-defined nodes
        for (const auto& [idx, status] : nodes_status)
        {
            container_impl<container_type>::check_size(temp_nodes_status, idx.first, idx.second);

            if (status == node_status::looped)
            {
                throw std::invalid_argument("node_status::looped is not allowed in "
                                            "'nodes_status' "
                                            "(use 'bounds_status' instead)");
            }
            else if (temp_nodes_status(idx.first, idx.second) == node_status::looped)
            {
                throw std::invalid_argument("cannot overwrite the status of a "
                                            "looped boundary node");
            }

            temp_nodes_status(idx.first, idx.second) = status;
        }

        m_nodes_status = temp_nodes_status;
    }

    /**
     * Pre-store for each (row, col) dimension 1-d arrays
     * that will be used to get the characteristic location of a node
     * on the grid.
     */
    template <class S, raster_connect RC, class C>
    void raster_grid<S, RC, C>::build_nodes_codes()
    {
        std::array<std::vector<code_type>, 2> gcode_rc;

        for (std::uint8_t dim = 0; dim < 2; ++dim)
        {
            auto fill_value = static_cast<std::uint8_t>(3 - dim * 2);
            std::vector<std::uint8_t> gcode_component(m_shape[dim], fill_value);

            gcode_component[0] = 0;
            gcode_component[m_shape[dim] - 1] = static_cast<std::uint8_t>(fill_value * 2);

            gcode_rc[dim] = gcode_component;
        }

        m_nodes_codes.resize({ m_size });
        for (std::size_t r = 0; r < m_shape[0]; ++r)
        {
            for (std::size_t c = 0; c < m_shape[1]; ++c)
            {
                m_nodes_codes[ravel_idx(r, c)]
                    = static_cast<std::uint8_t>(gcode_rc[0][r] + gcode_rc[1][c]);
            }
        }
    }

    /**
     * Pre-store the (row, col) index offsets of grid node neighbors for
     * each of the 9-characteristic locations on the grid.
     *
     * Those offsets take into account looped boundary conditions (if any).
     * The order of the returned offsets corresponds to the row-major layout.
     */
    template <class S, raster_connect RC, class C>
    auto raster_grid<S, RC, C>::build_coded_neighbors_offsets() -> coded_noffsets_type
    {
        auto dr = static_cast<std::ptrdiff_t>(m_shape[0] - 1);
        auto dc = static_cast<std::ptrdiff_t>(m_shape[1] - 1);

        if (!m_bounds_status.is_vertical_looped())
        {
            dr = 0;
        }
        if (!m_bounds_status.is_horizontal_looped())
        {
            dc = 0;
        }

        return {
            this->node_neighbors_offsets(dr, 1, dc, 1),      // top-left corner
            this->node_neighbors_offsets(dr, 1, -1, 1),      // top border
            this->node_neighbors_offsets(dr, 1, -1, -dc),    // top-right corner
            this->node_neighbors_offsets(-1, 1, dc, 1),      // left border
            this->node_neighbors_offsets(-1, 1, -1, 1),      // inner
            this->node_neighbors_offsets(-1, 1, -1, -dc),    // right border
            this->node_neighbors_offsets(-1, -dr, dc, 1),    // bottom-left corner
            this->node_neighbors_offsets(-1, -dr, -1, 1),    // bottom border
            this->node_neighbors_offsets(-1, -dr, -1, -dc),  // bottom-right corner
        };
    }

    template <class S, raster_connect RC, class C>
    inline auto raster_grid<S, RC, C>::neighbor_offsets(code_type code) const noexcept
        -> const neighbors_offsets_type&
    {
        return m_neighbor_offsets[code];
    }

    template <class S, raster_connect RC, class C>
    auto raster_grid<S, RC, C>::build_coded_neighbors_distances() -> coded_ndistances_type
    {
        coded_ndistances_type nb_distances;
        auto xspacing = m_spacing;

        auto to_dist = [&xspacing](auto&& offset) -> double
        { return container_impl<container_type>::compute_distance(offset, xspacing); };

        for (std::uint8_t k = 0; k < 9; ++k)
        {
            auto offsets = neighbor_offsets(k);
            auto distances = neighbors_distances_impl_type();

            std::transform(offsets.cbegin(), offsets.cend(), distances.begin(), to_dist);

            nb_distances[k] = distances;
        }

        return nb_distances;
    }

    /**
     * @name Node methods
     */
    /**
     * Given row and col indices, return a code in the range [0,8], which
     * corresponds to one of the following characteristic locations on the
     * grid (i.e., inner/border/corner). Use a row-major layout:
     *
     *   0 -- 1 -- 2
     *   |         |
     *   3    4    5
     *   |         |
     *   6 -- 7 -- 8
     */
    template <class S, raster_connect RC, class C>
    inline auto raster_grid<S, RC, C>::nodes_codes(const size_type& row,
                                                   const size_type& col) const noexcept -> code_type
    {
        return m_nodes_codes[ravel_idx(row, col)];
    }

    /**
     * Given a flat index, return a code in the range [0,8], which
     * corresponds to one of the following characteristic locations on the
     * grid (i.e., inner/border/corner). Use a row-major layout:
     *
     *   0 -- 1 -- 2
     *   |         |
     *   3    4    5
     *   |         |
     *   6 -- 7 -- 8
     */
    template <class S, raster_connect RC, class C>
    inline auto raster_grid<S, RC, C>::nodes_codes(const size_type& idx) const noexcept -> code_type
    {
        return m_nodes_codes[idx];
    }
    //@}

    /**
     * @name Neighbor methods
     */
    /**
     * Returns an array of the indices of the neighbors of a given grid node.
     *
     * Follows looped boundary conditions, if any.
     *
     * @param row The grid node row index.
     * @param col The grid node column index.
     */
    template <class S, raster_connect RC, class C>
    inline auto raster_grid<S, RC, C>::neighbors_indices(const size_type& row, const size_type& col)
        -> neighbors_indices_raster_type
    {
        neighbors_indices_raster_type indices;
        neighbors_indices(row, col, indices);

        return indices;
    }

    /**
     * Resize and fills an array with the neighbors indices of a given grid node.
     *
     * Follows looped boundary conditions, if any.
     *
     * This method prevents allocating a new container for better performance.
     *
     * @param row The grid node row index.
     * @param col The grid node column index.
     * @param neighbors_indices Reference to the container to be updated with the neighbors indices.
     */
    template <class S, raster_connect RC, class C>
    inline auto raster_grid<S, RC, C>::neighbors_indices(
        const size_type& row,
        const size_type& col,
        neighbors_indices_raster_type& neighbors_indices) -> neighbors_indices_raster_type&
    {
        const size_type flat_idx = ravel_idx(row, col);
        const auto& n_count = neighbors_count_impl(flat_idx);
        const auto& n_indices = this->get_nb_indices_from_cache(flat_idx);

        if (neighbors_indices.size() != n_count)
        {
            neighbors_indices.resize({ n_count });
        }

        for (size_type i = 0; i < n_count; ++i)
        {
            neighbors_indices[i] = unravel_idx(n_indices[i]);
        }

        return neighbors_indices;
    }

    /**
     * Returns a vector of the neighbors of a given grid node.
     *
     * Follows looped boundary conditions, if any.
     *
     * @param row The grid node row index.
     * @param col The grid node column index.
     * @return A vector of neighbor node objects.
     */
    template <class S, raster_connect RC, class C>
    inline auto raster_grid<S, RC, C>::neighbors(const size_type& row, const size_type& col)
        -> neighbors_raster_type
    {
        neighbors_raster_type nb;
        neighbors(row, col, nb);

        return nb;
    }

    /**
     * Resize and fills a vactor with the neighbors of a given grid node.
     *
     * Follows looped boundary conditions, if any.
     *
     * This method prevents allocating a new container for better performance.
     *
     * @param row The grid node row index.
     * @param col The grid node column index.
     * @param neighbors Reference to the vector to be updated with the neighbor objects.
     */
    template <class S, raster_connect RC, class C>
    inline auto raster_grid<S, RC, C>::neighbors(const size_type& row,
                                                 const size_type& col,
                                                 neighbors_raster_type& neighbors)
        -> neighbors_raster_type&
    {
        size_type n_flat_idx;
        raster_idx_type n_raster_idx;

        const size_type flat_idx = ravel_idx(row, col);
        const auto& n_count = neighbors_count_impl(flat_idx);
        const auto& n_indices = this->get_nb_indices_from_cache(flat_idx);
        const auto& n_distances = neighbors_distances_impl(flat_idx);

        if (neighbors.size() != n_count)
        {
            neighbors.resize({ n_count });
        }

        for (size_type i = 0; i < n_count; ++i)
        {
            n_flat_idx = n_indices[i];
            n_raster_idx = unravel_idx(n_flat_idx);
            neighbors[i] = raster_neighbor({ n_flat_idx,
                                             n_raster_idx.first,
                                             n_raster_idx.second,
                                             n_distances[i],
                                             this->nodes_status()(n_flat_idx) });
        }

        return neighbors;
    }
    //@}

    template <class S, raster_connect RC, class C>
    inline auto raster_grid<S, RC, C>::nodes_areas_impl() const -> container_type
    {
        return xt::broadcast(m_node_area, m_shape);
    }

    template <class S, raster_connect RC, class C>
    inline auto raster_grid<S, RC, C>::nodes_areas_impl(const size_type& /*idx*/) const noexcept
        -> grid_data_type
    {
        return m_node_area;
    }

    template <class S, raster_connect RC, class C>
    inline auto raster_grid<S, RC, C>::neighbors_count_impl(const size_type& idx) const noexcept
        -> size_type
    {
        return m_neighbors_count[m_nodes_codes[idx]];
    }

    template <class S, raster_connect RC, class C>
    inline auto raster_grid<S, RC, C>::neighbors_distances_impl(const size_type& idx) const noexcept
        -> const neighbors_distances_impl_type&
    {
        return m_neighbor_distances[nodes_codes(idx)];
    }

    template <class S, raster_connect RC, class C>
    inline auto raster_grid<S, RC, C>::neighbors_indices_impl(
        neighbors_indices_impl_type& neighbors, const size_type& idx) const -> void
    {
        const auto& offsets = neighbor_offsets(nodes_codes(idx));

        for (size_type i = 0; i < offsets.size(); ++i)
        {
            const auto offset = offsets[i];
            neighbors.at(i) = static_cast<size_type>((offset)[0]) * m_shape[1]
                              + static_cast<size_type>((offset)[1]) + idx;
        }
    }
}

#endif
