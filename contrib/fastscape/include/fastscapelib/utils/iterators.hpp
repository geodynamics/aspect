#ifndef FASTSCAPELIB_UTILS_ITERATORS_H
#define FASTSCAPELIB_UTILS_ITERATORS_H

#include <functional>

#include "xtl/xiterator_base.hpp"


namespace fastscapelib
{

    namespace detail
    {

        template <class G>
        struct grid_node_index_iterator
            : public xtl::xbidirectional_iterator_base<grid_node_index_iterator<G>,
                                                       typename G::size_type>
        {
        public:
            using self_type = grid_node_index_iterator<G>;
            using base_type = xtl::xbidirectional_iterator_base<self_type, typename G::size_type>;

            using value_type = typename base_type::value_type;
            using reference = typename base_type::reference;
            using pointer = typename base_type::pointer;
            using difference_type = typename base_type::difference_type;

            using filter_func_type = std::function<bool(const G&, typename G::size_type)>;

            grid_node_index_iterator(const G& grid, filter_func_type func, value_type position = 0)
                : m_idx(position)
                , m_grid(grid)
                , m_filter_func(func)
            {
                // iterate to the actual starting position (1st node that pass the filter test)
                while ((!m_filter_func(m_grid, m_idx)) && (m_idx < m_grid.size()))
                {
                    ++m_idx;
                }
            }

            inline self_type& operator++()
            {
                do
                {
                    ++m_idx;
                } while ((!m_filter_func(m_grid, m_idx)) && (m_idx < m_grid.size()));

                return *this;
            }

            inline self_type& operator--()
            {
                do
                {
                    --m_idx;
                } while ((!m_filter_func(m_grid, m_idx)) && (m_idx > 0));

                return *this;
            }

            inline reference operator*() const
            {
                return m_idx;
            }

        private:
            mutable value_type m_idx;

            const G& m_grid;
            filter_func_type m_filter_func;

            template <class _G>
            friend bool operator==(const grid_node_index_iterator<_G>&,
                                   const grid_node_index_iterator<_G>&);
        };

        template <class G>
        inline bool operator==(const grid_node_index_iterator<G>& lhs,
                               const grid_node_index_iterator<G>& rhs)
        {
            return lhs.m_idx == rhs.m_idx;
        }
    }


    /**
     * STL-compatible, immutable, virtual container for iterating through grid node
     * indices.
     *
     * Currently only implements a bidirectional iterator.
     *
     * @tparam G The grid type.
     */
    template <class G>
    class grid_nodes_indices
    {
    public:
        using filter_func_type = std::function<bool(const G&, typename G::size_type)>;
        using iterator = detail::grid_node_index_iterator<G>;

        grid_nodes_indices(const G& grid, filter_func_type func = nullptr)
            : m_grid(grid)
        {
            if (!func)
            {
                m_filter_func = [](const G&, typename G::size_type) { return true; };
            }
            else
            {
                m_filter_func = func;
            }
        }

        inline iterator begin() const
        {
            return iterator(m_grid, m_filter_func, 0);
        }

        inline iterator end() const
        {
            return iterator(m_grid, m_filter_func, m_grid.size());
        }

        inline std::reverse_iterator<iterator> rbegin() const
        {
            return std::reverse_iterator<iterator>(end());
        }

        inline std::reverse_iterator<iterator> rend() const
        {
            return std::reverse_iterator<iterator>(begin());
        }

    private:
        const G& m_grid;
        filter_func_type m_filter_func;
    };


    /**
     * Thin wrapper around iterators of any STL-compatible container.
     *
     * Only const forward and reverse iterators are exposed.
     *
     * @tparam C The container type.
     *
     */
    template <class C>
    class stl_container_iterator_wrapper
    {
    public:
        using const_iterator_type = typename C::const_iterator;
        using reverse_const_iterator_type = std::reverse_iterator<const_iterator_type>;

        stl_container_iterator_wrapper(const C& container)
            : m_container(container)
        {
        }

        inline const_iterator_type begin() const
        {
            return m_container.begin();
        }

        inline const_iterator_type end() const
        {
            return m_container.end();
        }

        inline reverse_const_iterator_type rbegin() const
        {
            return m_container.rbegin();
        }

        inline reverse_const_iterator_type rend() const
        {
            return m_container.rend();
        }

    private:
        const C& m_container;
    };
}
#endif
