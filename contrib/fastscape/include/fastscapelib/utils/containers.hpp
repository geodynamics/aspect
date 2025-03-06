#ifndef FASTSCAPELIB_UTILS_CONTAINERS_HPP
#define FASTSCAPELIB_UTILS_CONTAINERS_HPP

#include <cstddef>


namespace fastscapelib
{
    /**
     * Used to get the actual container type from a given selector.
     *
     * @tparam S The selector type.
     * @tparam T The container value type.
     * @tparam N The number of dimensions required (only for static dimension containers)
     */
    template <class S, class T, std::size_t N = 0>
    struct container_selection
    {
    };

    /**
     * Used to get the details (types and implementations) of a container type.
     *
     * @tparam C The container type.
     */
    template <class C>
    struct container_impl
    {
    };

    /**
     * Alias for the selected static dimension container type.
     *
     * @tparam S The container selector type.
     * @tparam T The container value type.
     * @tparam N The fixed number of dimensions.
     */
    template <class S, class T, std::size_t N>
    using fixed_shape_container_t = typename container_selection<S, T, N>::fixed_shape_type;

    /**
     * Alias for the selected dynamic dimension container type.
     *
     * @tparam S The container selector type.
     * @tparam T The container value type.
     */
    template <class S, class T>
    using dynamic_shape_container_t = typename container_selection<S, T>::dynamic_shape_type;
}  // namespace fastscapelib

#include "fastscapelib/utils/xtensor_containers.hpp"

#endif  // FASTSCAPELIB_UTILS_CONTAINERS_HPP
