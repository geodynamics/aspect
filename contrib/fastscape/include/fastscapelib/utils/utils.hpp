/**
 * @file
 * @brief Some utilities used internally.
 */
#ifndef FASTSCAPELIB_UTILS_UTILS_H
#define FASTSCAPELIB_UTILS_UTILS_H


#include <cstddef>
#include <type_traits>

#include "xtensor.hpp"
#include "xtensor/xcontainer.hpp"


// type used for array indexers
using index_t = std::ptrdiff_t;

// type used when working with both shapes or sizes (unsigned) and indexers (signed)
using sshape_t = std::make_signed_t<std::size_t>;

// type used as an alias to xtensor's xexpression
// TODO: use xexpression shaped when it's ready (see xtensor GH #994)
template <class E>
using xtensor_t = xt::xexpression<E>;


namespace fastscapelib
{
    namespace detail
    {

        /**
         * @brief Return true if a given (row, col) index is in bounds (false
         *        otherwise).
         */
        template <class S>
        bool in_bounds(const S& shape, index_t row, index_t col)
        {
            return (row >= 0 && row < static_cast<index_t>(shape[0]) && col >= 0
                    && col < static_cast<index_t>(shape[1]));
        }

        /**
         * @brief Return a pair of the (row, col) coordinate of a node
         *        in a dem of ncols columns
         */
        template <class Node_T>
        std::pair<Node_T, Node_T> coords(Node_T node, Node_T ncols)
        {
            return { node / ncols, node % ncols };
        }

        /**
         * @brief Return the node index of a node given its (row, col) coordinate
         *        in a dem of ncols columns
         */
        template <class Node_T>
        Node_T index(Node_T row, Node_T col, Node_T ncols)
        {
            return col + row * ncols;
        }

        /**
         * @brief proxy for flattened view waiting it to be implmented in xtensor
         * Unsafe in the genral case, use it as a replacement for maintance help ...
         */

        template <class Array_T>
        class Flattened2D
        {
        public:
            using value_type = typename Array_T::value_type;
            using size_type = typename Array_T::size_type;

            Flattened2D(Array_T& ref)
                : _array_ref{ ref }
            {
                _ncols = ref.shape()[1];
            }

            typename Array_T::value_type operator()(size_type index)
            {
                auto c = coords(index, _ncols);
                return _array_ref(std::get<0>(c), std::get<1>(c));
            }

            typename Array_T::value_type at(size_type index)
            {
                auto c = coords(index, _ncols);
                return _array_ref.at(std::get<0>(c), std::get<1>(c));
            }

        private:
            Array_T& _array_ref;
            size_type _ncols;
        };

        template <class Array_T>
        Flattened2D<Array_T> make_flattened(Array_T& ref)
        {
            return Flattened2D<Array_T>(ref);
        }
    }
}

#endif
