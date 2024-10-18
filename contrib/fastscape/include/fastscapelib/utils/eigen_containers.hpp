/**
 * xtensor utils (container tags, etc.)
 */
#ifndef FASTSCAPELIB_UTILS_EIGEN_CONTAINERS_HPP
#define FASTSCAPELIB_UTILS_EIGEN_CONTAINERS_HPP

#include "fastscapelib/utils/containers.hpp"

#include <Eigen/Dense>

#include <stdexcept>


namespace fastscapelib
{
    struct eigen_selector
    {
    };

    template <class T, std::size_t N>
    struct container_selection<eigen_selector, T, N>
    {
        // static_assert(N <= 2 && N > 0, "Not handled");
        using fixed_shape_type = Eigen::Array<T, -1, -1>;
        using dynamic_shape_type = Eigen::Array<T, -1, -1>;
    };

    template <class T, int R, int C>
    struct container_impl<Eigen::Array<T, R, C>>
    {
        using container_type = typename Eigen::Array<T, R, C>;
        using size_type = std::size_t;
        using shape_type = std::array<std::size_t, 2>;

        template <class S, class V>
        static auto init(S&& shape, V&& value)
        {
            return container_type::Constant(shape[0], shape[1], value);
        }

        template <class O, class X>
        static double compute_distance(O& offset, X& xspacing)
        {
            using namespace Eigen;

            auto drc = (Map<const Array<typename std::remove_reference_t<O>::value_type, 2, 1>>(
                            offset.data())
                        == 0)
                           .select(container_type::Zero(2, 1), container_type::Constant(2, 1, 1))
                       * Map<const Array<typename std::remove_reference_t<X>::value_type, 2, 1>>(
                           xspacing.data());
            return std::sqrt(drc.square().sum());
        }

        template <class D, class I>
        static auto get_view(D&& data, I&& max_idx)
        {
            using namespace Eigen;

            return Map<const Array<typename std::remove_reference_t<D>::value_type, -1, 1>>(
                data.data())(seq(0, max_idx));
        }

        template <class D>
        static auto get_top_view(D&& data)
        {
            return data(0, Eigen::all);
        }

        template <class D>
        static auto get_bottom_view(D&& data)
        {
            return data(Eigen::last, Eigen::all);
        }

        template <class D>
        static auto get_left_view(D&& data)
        {
            return data(Eigen::all, 0);
        }

        template <class D>
        static auto get_right_view(D&& data)
        {
            return data(Eigen::all, Eigen::last);
        }

        template <class D, class I>
        static void check_size(D& data, I& row_index, I& col_index)
        {
            if (Eigen::Index(row_index) >= data.rows() || Eigen::Index(col_index) >= data.cols())
                throw std::out_of_range("Invalid index");
        }
    };
}  // namespace fastscapelib

#endif  // FASTSCAPELIB_UTILS_EIGEN_CONTAINERS_HPP
