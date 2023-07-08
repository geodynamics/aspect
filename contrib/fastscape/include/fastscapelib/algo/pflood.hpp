/**
 * Provides implementation of priority-flood algorithms for
 * depression filling or pit resolving.
 */
#ifndef FASTSCAPELIB_ALGO_PFLOOD_H
#define FASTSCAPELIB_ALGO_PFLOOD_H

#include <cmath>
#include <functional>
#include <algorithm>
#include <queue>
#include <limits>
#include <type_traits>

#include "xtensor/xtensor.hpp"

#include "fastscapelib/utils/utils.hpp"
#include "fastscapelib/utils/consts.hpp"
#include "fastscapelib/grid/structured_grid.hpp"


namespace fastscapelib
{
    namespace detail
    {

        /**
         * Grid node structure.
         *
         * Stores both grid position and elevation at that position.
         * Also defines  operator '>' that compares only on `elevation`.
         *
         * The main purpose of this container is for using with
         * priority-flood algorithms.
         */
        template <class G, class T>
        struct pflood_node
        {
            using size_type = typename G::size_type;

            size_type m_idx;
            T m_elevation;

            pflood_node()
            {
            }
            pflood_node(size_type idx, T elevation)
                : m_idx(idx)
                , m_elevation(elevation)
            {
            }

            bool operator>(const pflood_node<G, T>& other) const
            {
                return m_elevation > other.m_elevation;
            }
        };


        template <class G, class T>
        using pflood_pr_queue = std::priority_queue<pflood_node<G, T>,
                                                    std::vector<pflood_node<G, T>>,
                                                    std::greater<pflood_node<G, T>>>;


        template <class G, class T>
        using pflood_queue = std::queue<pflood_node<G, T>>;


        /**
         * Initialize priority flood algorithms.
         *
         * Add fixed value grid nodes to the priority queue and mark them as
         * resolved.
         */
        template <class G, class E, class elev_t = typename std::decay_t<E>::value_type>
        void init_pflood(G& grid,
                         E&& elevation,
                         xt::xtensor<bool, 1>& closed,
                         pflood_pr_queue<G, elev_t>& open)
        {
            using size_type = typename G::size_type;

            // TODO: assert elevation shape match grid shape

            const auto elevation_flat = xt::flatten(elevation);

            for (size_type idx = 0; idx < grid.size(); ++idx)
            {
                if (grid.nodes_status()[idx] == fastscapelib::node_status::fixed_value)
                {
                    open.emplace(pflood_node<G, elev_t>(idx, elevation_flat(idx)));
                    closed(idx) = true;
                }
            }
        }


        /**
         * fill_sinks_flat implementation.
         */
        template <class G, class E>
        void fill_sinks_flat_impl(G& grid, E&& elevation)
        {
            using neighbors_indices_type = typename G::neighbors_indices_type;
            using elev_t = typename std::decay_t<E>::value_type;

            auto elevation_flat = xt::flatten(elevation);
            neighbors_indices_type neighbors_indices;

            pflood_pr_queue<G, elev_t> open;
            xt::xtensor<bool, 1> closed = xt::zeros<bool>({ grid.size() });

            init_pflood(grid, elevation, closed, open);

            while (open.size() > 0)
            {
                pflood_node<G, elev_t> inode = open.top();
                open.pop();

                for (auto n_idx : grid.neighbors_indices(inode.m_idx, neighbors_indices))
                {
                    if (closed(n_idx))
                    {
                        continue;
                    }

                    elevation_flat(n_idx) = std::max(elevation_flat(n_idx), inode.m_elevation);
                    open.emplace(pflood_node<G, elev_t>(n_idx, elevation_flat(n_idx)));
                    closed(n_idx) = true;
                }
            }
        }


        /**
         * fill_sinks_sloped implementation.
         */
        template <class G, class E>
        void fill_sinks_sloped_impl(G& grid, E&& elevation)
        {
            using neighbors_indices_type = typename G::neighbors_indices_type;
            using elev_t = typename std::decay_t<E>::value_type;

            neighbors_indices_type neighbors_indices;

            pflood_pr_queue<G, elev_t> open;
            pflood_queue<G, elev_t> pit;
            xt::xtensor<bool, 1> closed = xt::zeros<bool>({ grid.size() });

            init_pflood(grid, elevation, closed, open);

            while (!open.empty() || !pit.empty())
            {
                pflood_node<G, elev_t> inode, knode;

                if (!pit.empty()
                    && (open.empty() || open.top().m_elevation == pit.front().m_elevation))
                {
                    inode = pit.front();
                    pit.pop();
                }
                else
                {
                    inode = open.top();
                    open.pop();
                }

                elev_t elev_tiny_step
                    = std::nextafter(inode.m_elevation, std::numeric_limits<elev_t>::infinity());

                for (auto n_idx : grid.neighbors_indices(inode.m_idx, neighbors_indices))
                {
                    if (closed(n_idx))
                    {
                        continue;
                    }

                    if (elevation.flat(n_idx) <= elev_tiny_step)
                    {
                        elevation.flat(n_idx) = elev_tiny_step;
                        knode = pflood_node<G, elev_t>(n_idx, elevation.flat(n_idx));
                        pit.emplace(knode);
                    }
                    else
                    {
                        knode = pflood_node<G, elev_t>(n_idx, elevation.flat(n_idx));
                        open.emplace(knode);
                    }

                    closed(n_idx) = true;
                }
            }
        }

    }  // namespace detail


    /**
     * Fill all pits and remove all digital dams from a digital
     * elevation model (rectangular grid).
     *
     * Elevation values may be updated so that no depression (sinks or
     * local minima) remains.
     *
     * The algorithm is based on priority queues and is detailed
     * in Barnes (2014). It fills sinks with flat areas.
     *
     * @param grid Grid object
     * @param elevation Elevation values at grid nodes
     */
    template <class G, class E>
    void fill_sinks_flat(G& grid, xtensor_t<E>& elevation)
    {
        detail::fill_sinks_flat_impl(grid, elevation.derived_cast());
    }


    /**
     * Fill all pits and remove all digital dams from a digital
     * elevation model (rectangular grid).
     *
     * Elevation values may be updated so that no depression (sinks or
     * local minima) remains.
     *
     * The algorithm is based on priority queues and is detailed in Barnes
     * (2014). This variant fills sinks with nearly flat areas
     * (i.e. elevation is increased by small amount) so that there is no
     * drainage singularities.
     *
     * @param grid Grid object
     * @param elevation Elevation values at grid nodes
     *
     * @sa fill_sinks_flat
     */
    template <class G, class E>
    void fill_sinks_sloped(G& grid, xtensor_t<E>& elevation)
    {
        detail::fill_sinks_sloped_impl(grid, elevation.derived_cast());
    }
}  // namespace fastscapelib

#endif
