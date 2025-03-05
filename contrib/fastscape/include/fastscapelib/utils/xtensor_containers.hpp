/**
 * xtensor utils (container tags, etc.)
 */
#ifndef FASTSCAPELIB_UTILS_XTENSOR_CONTAINERS_HPP
#define FASTSCAPELIB_UTILS_XTENSOR_CONTAINERS_HPP

#include "xtensor/xarray.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xadapt.hpp"
#include "xtensor/xview.hpp"

#include <stdexcept>


namespace fastscapelib
{
    /**
     * The xtensor container selector, used by default in Fastscapelib C++ API.
     *
     * \rst
     *
     * Fastscapelib template classes specialized with this selector will use
     * :cpp:type:`xt::xtensor` or :cpp:type:`xt::xarray` as container types for
     * their public array members.
     *
     * \endrst
     *
     */
    struct xt_selector
    {
    };

    /**
     * Defines the container types for the `xt_selector` container selector.
     */
    template <class T, std::size_t N>
    struct container_selection<xt_selector, T, N>
    {
        using fixed_shape_type = xt::xtensor<T, N>;
        using dynamic_shape_type = xt::xarray<T>;
    };

    template <class CT>
    struct xtensor_shared_utils
    {
        using container_type = CT;

        template <class O, class X>
        static double compute_distance(O&& offset, X&& xspacing)
        {
            auto drc = xt::where(xt::equal(xt::adapt(offset), 0), 0., 1.) * xt::adapt(xspacing);
            return std::sqrt(xt::sum(xt::square(drc))(0));
        }

        template <class T, class I>
        static auto get_view(T&& data, I&& max_idx)
        {
            return xt::view(xt::adapt(std::forward<T>(data)),
                            xt::range(0, std::forward<I>(max_idx)));
        }

        template <class T>
        static auto get_top_view(T&& data)
        {
            return xt::view(data, 0, xt::all());
        }

        template <class T>
        static auto get_bottom_view(T&& data)
        {
            return xt::view(data, xt::keep(-1), xt::all());
        }

        template <class T>
        static auto get_left_view(T&& data)
        {
            return xt::view(data, xt::all(), 0);
        }

        template <class T>
        static auto get_right_view(T&& data)
        {
            return xt::view(data, xt::all(), xt::keep(-1));
        }

        template <class T, class I>
        static void check_size(T& data, I& row_index, I& col_index)
        {
            if (row_index >= data.shape(0) || col_index >= data.shape(1))
                throw std::out_of_range("Invalid index");
        }

        template <class S, class V>
        static auto init(S&& shape, V&& value)
        {
            return container_type(shape, value);
        }
    };

    template <class T>
    struct container_impl<xt::xarray<T>> : public xtensor_shared_utils<xt::xarray<T>>
    {
        using base_type = xtensor_shared_utils<xt::xarray<T>>;
        using container_type = typename xt::xarray<T>;
        using size_type = typename container_type::size_type;
        using shape_type = typename container_type::shape_type;

        using base_type::compute_distance;
        using base_type::get_view;
        using base_type::init;
    };

    template <class T, std::size_t N>
    struct container_impl<xt::xtensor<T, N>> : public xtensor_shared_utils<xt::xtensor<T, N>>
    {
        using base_type = xtensor_shared_utils<xt::xtensor<T, N>>;
        using container_type = typename xt::xtensor<T, N>;
        using size_type = typename container_type::size_type;
        using shape_type = typename container_type::shape_type;

        using base_type::compute_distance;
        using base_type::get_view;
        using base_type::init;
    };
}  // namespace fastscapelib

#endif  // FASTSCAPELIB_UTILS_XTENSOR_CONTAINERS_HPP
