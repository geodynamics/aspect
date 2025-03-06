/**
 * structured grid abstract class.
 */
#ifndef FASTSCAPELIB_GRID_STRUCTURED_GRID_H
#define FASTSCAPELIB_GRID_STRUCTURED_GRID_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <stdexcept>
#include <type_traits>
#include <map>
#include <vector>

#include "xtensor/xbuilder.hpp"
#include "xtensor/xarray.hpp"
#include "xtensor/xfixed.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xview.hpp"
#include "xtensor/xindex_view.hpp"
#include "xtensor/xnoalias.hpp"

#include "xtl/xiterator_base.hpp"

#include "fastscapelib/grid/base.hpp"
#include "fastscapelib/utils/utils.hpp"
#include "fastscapelib/utils/containers.hpp"


namespace fastscapelib
{

    /**
     * Base class for setting the status at the border nodes of structured grids.
     */
    class boundary_status
    {
    protected:
        bool is_looped(node_status status) const;
    };

    inline bool boundary_status::is_looped(const node_status status) const
    {
        return status == node_status::looped;
    }

    /**
     * Extends the common grid interface for all structured grid types.
     *
     * @tparam G The derived grid type.
     */
    template <class G>
    class structured_grid : public grid<G>
    {
    public:
        using base_type = grid<G>;
        using inner_types = grid_inner_types<G>;

        using shape_type = typename base_type::shape_type;
        using length_type = typename inner_types::length_type;
        using spacing_type = typename inner_types::spacing_type;

        using spacing_t = std::conditional_t<std::is_arithmetic<spacing_type>::value,
                                             spacing_type,
                                             const spacing_type&>;

        spacing_t spacing() const noexcept;
        length_type length() const noexcept;
        shape_type shape() const noexcept;

    protected:
        using grid<G>::grid;
        ~structured_grid() = default;
    };

    /**
     * @name Grid properties
     */
    //@{
    /**
     * Returns the (uniform) spacing between two adjacent grid nodes.
     *
     * Depending on the dimensions of the grid, returns either a single value
     * or an array (constant reference).
     */
    template <class G>
    inline auto structured_grid<G>::spacing() const noexcept -> spacing_t
    {
        return this->derived_grid().m_spacing;
    }

    /**
     * Returns the length of the grid for all its dimensions.
     */
    template <class G>
    inline auto structured_grid<G>::length() const noexcept -> length_type
    {
        return this->derived_grid().m_length;
    }

    /**
     * Returns the shape of the grid node arrays.
     */
    template <class G>
    inline auto structured_grid<G>::shape() const noexcept -> shape_type
    {
        return this->derived_grid().m_shape;
    }
    //@}
}

#endif
