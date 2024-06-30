/**
 * xtensor utils (container tags, etc.)
 */
#ifndef FASTSCAPELIB_UTILS_XTENSOR_UTILS_H
#define FASTSCAPELIB_UTILS_XTENSOR_UTILS_H


#include <cstddef>

#include "xtensor/xarray.hpp"
#include "xtensor/xtensor.hpp"


namespace fastscapelib
{

    /**
     * The xtensor selector used by default in Fastscapelib C++ API.
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
     * Used to get the actual xtensor container type from a given selector.
     *
     * @tparam S The xtensor selector type.
     * @tparam T The container value type.
     * @tparam N The number of dimensions (only for static dimension containers)
     */
    template <class S, class T, std::size_t N = 0>
    struct xt_container
    {
    };

    template <class T, std::size_t N>
    struct xt_container<xt_selector, T, N>
    {
        using tensor_type = xt::xtensor<T, N>;
        using array_type = xt::xarray<T>;
    };

    /**
     * Alias for the selected (static dimension) xtensor container type.
     *
     * @tparam S The xtensor selector type.
     * @tparam T The container value type.
     * @tparam N The fixed number of dimensions.
     */
    template <class S, class T, std::size_t N>
    using xt_tensor_t = typename xt_container<S, T, N>::tensor_type;

    /**
     * Alias for the selected (dynamic dimension) xtensor container type.
     *
     * @tparam S The xtensor selector type.
     * @tparam T The container value type.
     */
    template <class S, class T>
    using xt_array_t = typename xt_container<S, T>::array_type;

}  // namespace fastscapelib

#endif
