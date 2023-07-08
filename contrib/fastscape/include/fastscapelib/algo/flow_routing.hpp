/**
 * @brief Functions used to route (water) flow on a topographic
 * surface and compute flow path-related features or structures.
 *
 */
#ifndef FASTSCAPELIB_ALGO_FLOW_ROUTING_H
#define FASTSCAPELIB_ALGO_FLOW_ROUTING_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <type_traits>

#include "xtensor/xtensor.hpp"
#include "xtensor/xview.hpp"
#include "xtensor/xmanipulation.hpp"

#include "fastscapelib/utils/utils.hpp"
#include "fastscapelib/utils/consts.hpp"


namespace fastscapelib
{
    namespace detail
    {


        inline auto get_d8_distances(double dx, double dy) -> std::array<double, 9>
        {
            std::array<double, 9> d8_dists;

            for (std::size_t k = 0; k < 9; ++k)
            {
                d8_dists[k]
                    = std::sqrt(std::pow(dy * fastscapelib::consts::d8_row_offsets[k], 2.0)
                                + std::pow(dx * fastscapelib::consts::d8_col_offsets[k], 2.0));
            }

            return d8_dists;
        }


        template <class S, class N, class D>
        void add2stack(index_t& nstack, S&& stack, N&& ndonors, D&& donors, index_t inode)
        {
            for (index_t k = 0; k < ndonors(inode); ++k)
            {
                const auto idonor = donors(inode, k);
                stack(nstack) = idonor;
                ++nstack;
                add2stack(nstack, stack, ndonors, donors, idonor);
            }
        }


        /**
         * compute_receivers_d8 implementation.
         */
        template <class R, class D, class E, class A>
        void compute_receivers_d8_impl(R&& receivers,
                                       D&& dist2receivers,
                                       E&& elevation,
                                       A&& active_nodes,
                                       double dx,
                                       double dy)
        {
            using elev_t = typename std::decay_t<E>::value_type;
            const auto d8_dists = detail::get_d8_distances(dx, dy);

            const auto elev_shape = elevation.shape();
            const auto nrows = static_cast<index_t>(elev_shape[0]);
            const auto ncols = static_cast<index_t>(elev_shape[1]);

            for (index_t r = 0; r < nrows; ++r)
            {
                for (index_t c = 0; c < ncols; ++c)
                {
                    const index_t inode = r * ncols + c;

                    receivers(inode) = inode;
                    dist2receivers(inode) = 0.;

                    if (!active_nodes(r, c))
                    {
                        continue;
                    }

                    double slope_max = std::numeric_limits<double>::min();

                    for (std::size_t k = 1; k <= 8; ++k)
                    {
                        const index_t kr = r + fastscapelib::consts::d8_row_offsets[k];
                        const index_t kc = c + fastscapelib::consts::d8_col_offsets[k];

                        if (!fastscapelib::detail::in_bounds(elev_shape, kr, kc))
                        {
                            continue;
                        }

                        const index_t ineighbor = kr * ncols + kc;
                        const double slope = (elevation(r, c) - elevation(kr, kc)) / d8_dists[k];

                        if (slope > slope_max)
                        {
                            slope_max = slope;
                            receivers(inode) = ineighbor;
                            dist2receivers(inode) = d8_dists[k];
                        }
                    }
                }
            }
        }


        /**
         * compute_donors implementation.
         */
        template <class N, class D, class R>
        void compute_donors_impl(N&& ndonors, D&& donors, R&& receivers)
        {
            const auto nnodes = static_cast<index_t>(receivers.size());

            std::fill(ndonors.begin(), ndonors.end(), 0);

            for (index_t inode = 0; inode < nnodes; ++inode)
            {
                if (receivers(inode) != inode)
                {
                    index_t irec = receivers(inode);
                    donors(irec, ndonors(irec)) = inode;
                    ++ndonors(irec);
                }
            }
        }


        /**
         * compute_stack implementation.
         */
        template <class S, class N, class D, class R>
        void compute_stack_impl(S&& stack, N&& ndonors, D&& donors, R&& receivers)
        {
            const auto nnodes = static_cast<index_t>(receivers.size());
            index_t nstack = 0;

            for (index_t inode = 0; inode < nnodes; ++inode)
            {
                if (receivers(inode) == inode)
                {
                    stack(nstack) = inode;
                    ++nstack;
                    add2stack(nstack, stack, ndonors, donors, inode);
                }
            }
        }


        /**
         * compute_basins implementation.
         */
        template <class B, class O, class S, class R>
        index_t compute_basins_impl(B&& basins, O&& outlets_or_pits, S&& stack, R&& receivers)
        {
            index_t ibasin = -1;

            for (auto&& istack : stack)
            {
                const auto irec = receivers(istack);

                if (irec == istack)
                {
                    ++ibasin;
                    outlets_or_pits(ibasin) = istack;
                }

                basins(istack) = ibasin;
            }

            index_t nbasins = ibasin + 1;

            return nbasins;
        }


        /**
         * find_pits implementation.
         */
        template <class P, class O, class A>
        index_t find_pits_impl(P&& pits, O&& outlets_or_pits, A&& active_nodes, index_t nbasins)
        {
            index_t ipit = 0;
            const auto active_nodes_flat = xt::flatten(active_nodes);

            for (index_t ibasin = 0; ibasin < nbasins; ++ibasin)
            {
                const index_t inode = outlets_or_pits(ibasin);

                if (active_nodes_flat(inode))
                {
                    pits(ipit) = inode;
                    ++ipit;
                }
            }

            index_t npits = ipit;

            return npits;
        }


        /**
         * compute_drainage_area implementation.
         */
        template <class D, class C, class S, class R>
        void compute_drainage_area_impl(D&& drainage_area, C&& cell_area, S&& stack, R&& receivers)
        {
            // reset drainage area values (must use a view to prevent resizing
            // drainage_area to 0-d when cell_area is 0-d!)
            auto drainage_area_ = xt::view(drainage_area, xt::all(), xt::all());
            drainage_area_ = cell_area;

            // update drainage area values
            auto drainage_area_flat = xt::flatten(drainage_area);

            for (auto inode = stack.crbegin(); inode != stack.crend(); ++inode)
            {
                if (receivers(*inode) != *inode)
                {
                    drainage_area_flat(receivers(*inode)) += drainage_area_flat(*inode);
                }
            }
        }
    }  // namespace detail


    /**
     * Compute flow receivers on a rectangular grid using D8 single flow
     * routing method.
     *
     * Each node on the grid is assigned one receiver among its 8 direct
     * neighboring nodes according to the steepest slope (O’Callaghan and
     * Mark, 1984).
     *
     * When no downslope neighbor exist (i.e., pit nodes or grid
     * boundaries), the assigned receiver is the node itself. When two or
     * more neighbors have the same slope, the chosen neighbor is the
     * first one considered by the algorithm.
     *
     * This function also computes the planimetric distance between the
     * node and its receiver, which equals to the grid spacing in x or y
     * or to the distance between two diagonal neighbor nodes or 0.
     *
     * @param receivers : ``[intent=out, shape=(nnodes)]``
     *     Index of flow receiver at grid node.
     * @param dist2receivers : ``[intent=out, shape=(nnodes)]``
     *     Distance to receiver at grid node.
     * @param elevation : ``[intent=in, shape=(nrows, ncols)]``
     *     Topographic elevation at grid node
     * @param active_nodes : ``[intent=in, shape=(nrows, ncols)]``
     *     Boolean array for boundaries
     * @param dx : ``[intent=in]``
     *     Grid spacing in x
     * @param dy : ``[intent=in]``
     *     Grid spacing in y
     */
    template <class R, class D, class E, class A>
    void compute_receivers_d8(xtensor_t<R>& receivers,
                              xtensor_t<D>& dist2receivers,
                              const xtensor_t<E>& elevation,
                              const xtensor_t<A>& active_nodes,
                              double dx,
                              double dy)
    {
        detail::compute_receivers_d8_impl(receivers.derived_cast(),
                                          dist2receivers.derived_cast(),
                                          elevation.derived_cast(),
                                          active_nodes.derived_cast(),
                                          dx,
                                          dy);
    }


    namespace detail
    {
        template <class R, class D, class E, class G>
        void compute_receivers_impl(R&& receivers, D&& dist2receivers, E&& elevation, G& grid)
        {
            using neighbors_type = typename G::neighbors_type;

            double slope, slope_max;
            neighbors_type neighbors;

            for (std::size_t i = 0; i < grid.size(); ++i)
            {
                receivers(i, 0) = i;
                dist2receivers(i, 0) = 0;
                slope_max = std::numeric_limits<double>::min();

                grid.neighbors(i, neighbors);

                for (auto n = neighbors.begin(); n != neighbors.end(); ++n)
                {
                    slope = (elevation.data()[i] - elevation.data()[n->idx]) / n->distance;

                    if (slope > slope_max)
                    {
                        slope_max = slope;
                        receivers(i, 0) = n->idx;
                        dist2receivers(i, 0) = n->distance;
                    }
                }
            }
        }
    }

    template <class R, class D, class E, class G>
    void compute_receivers(xtensor_t<R>& receivers,
                           xtensor_t<D>& dist2receivers,
                           const xtensor_t<E>& elevation,
                           G& grid)
    {
        detail::compute_receivers_impl(receivers.derived_cast(),
                                       dist2receivers.derived_cast(),
                                       elevation.derived_cast(),
                                       grid);
    }


    /**
     * Compute flow donors for each grid/mesh node.
     *
     * Flow donors are retrieved by simply inverting flow
     * receivers.
     *
     * @param ndonors : ``[intent=out, shape=(nnodes)]``
     *     Number of flow donors at grid node.
     * @param donors : ``[intent=out, shape=(nnodes, :)]``
     *     Indexes of flow donors at grid node.
     * @param receivers : ``[intent=in, shape=(nnodes)]``
     *     Index of flow receiver at grid node.
     */
    template <class N, class D, class R>
    void compute_donors(xtensor_t<N>& ndonors, xtensor_t<D>& donors, const xtensor_t<R>& receivers)
    {
        detail::compute_donors_impl(
            ndonors.derived_cast(), donors.derived_cast(), receivers.derived_cast());
    }


    /**
     * Compute a stack of grid/mesh nodes to be used for flow tree
     * traversal.
     *
     * The stack is calculated recursively from outlets (or sinks) to
     * sources, using Braun and Willet's (2013) algorithm.
     *
     * @param stack : ``[intent=out, shape=(nnodes)]``
     *     Stack position at grid node.
     * @param ndonors : ``[intent=in, shape=(nnodes)]``
     *     Number of flow donors at grid node.
     * @param donors : ``[intent=in, shape=(nnodes, :)]``
     *     Indexes of flow donors at grid node.
     * @param receivers : ``[intent=in, shape=(nnodes)]``
     *     Index of flow receiver at grid node.
     */
    template <class S, class N, class D, class R>
    void compute_stack(xtensor_t<S>& stack,
                       const xtensor_t<N>& ndonors,
                       const xtensor_t<D>& donors,
                       const xtensor_t<R>& receivers)
    {
        detail::compute_stack_impl(stack.derived_cast(),
                                   ndonors.derived_cast(),
                                   donors.derived_cast(),
                                   receivers.derived_cast());
    }


    /**
     * Assign an id (integer) to each node of the grid/mesh that
     * corresponds to the catchment to which it belongs.
     *
     * A catchment (or drainage basin) is defined by an ensemble of
     * adjacent nodes through which all flow converges towards a common,
     * single node (outlet or pit).
     *
     * The algorithm performs a single traversal of the flow tree (in the
     * stack order) and increments the catchment id each time an outlet or
     * a pit is found (i.e., when the index of the flow receiver equals
     * the index of the node itself).
     *
     * This functions also computes the grid/mesh node indexes of
     * catchment outlets (or pits) and returns the total number of
     * catchments found inside the domain.
     *
     * @param basins: ``[intent=out, shape=(nnodes)]``
     *     Basin id at grid node.
     * @param outlets_or_pits : ``[intent=out, shape=(nnodes)]``
     *     Grid node index of the outlet (or pit)
     *     for basin id=0,1,...,nbasins-1.
     * @param stack :``[intent=in, shape=(nnodes)]``
     *     Stack position at grid node.
     * @param receivers : ``[intent=in, shape=(nnodes)]``
     *     Index of flow receiver at grid node.
     *
     * @returns
     *     Total number of drainage basins
     *     (``1 <= nbasins <= nnodes``).
     */
    template <class B, class O, class S, class R>
    index_t compute_basins(xtensor_t<B>& basins,
                           xtensor_t<O>& outlets_or_pits,
                           const xtensor_t<S>& stack,
                           const xtensor_t<R>& receivers)
    {
        return detail::compute_basins_impl(basins.derived_cast(),
                                           outlets_or_pits.derived_cast(),
                                           stack.derived_cast(),
                                           receivers.derived_cast());
    }


    /**
     * Find grid/mesh nodes that are pits.
     *
     * @param pits: ``[intent=out, shape=(nnodes)]``
     *     Grid node index of the pit for pits in [0,1,...,npits-1].
     * @param outlets_or_pits : ``[intent=in, shape=(nnodes)]``
     *     Grid node index of the outlet (or pit)
     *     for basin id=0,1,...,nbasins-1.
     * @param active_nodes : ``[intent=in, shape=(nrows, ncols)]``
     *     Boolean array for boundaries
     * @param nbasins : ``[intent=in]``
     *     Total number of drainage basins (``1 <= nbasins <= nnodes``).
     *
     * @returns
     *     Total number of pits found (``0 <= npits <= nbasins``).
     */
    template <class P, class O, class A>
    index_t find_pits(xtensor_t<P>& pits,
                      const xtensor_t<O>& outlets_or_pits,
                      const xtensor_t<A>& active_nodes,
                      index_t nbasins)
    {
        return detail::find_pits_impl(pits.derived_cast(),
                                      outlets_or_pits.derived_cast(),
                                      active_nodes.derived_cast(),
                                      nbasins);
    }


    /**
     * Compute drainage area for each node on a generic grid/mesh.
     *
     * @param drainage_area : ``[intent=out, shape=(nrows, ncols)||(nnodes)]``
     *     Drainage area at grid node.
     * @param cell_area : ``[intent=in, shape=(nrows, ncols)||(nnodes)||()]``
     *     Grid/mesh cell area at grid node (also accepts a scalar).
     * @param stack :``[intent=in, shape=(nnodes)]``
     *     Stack position at grid node.
     * @param receivers : ``[intent=in, shape=(nnodes)]``
     *     Index of flow receiver at grid node.
     */
    template <class D, class C, class S, class R>
    void compute_drainage_area(xtensor_t<D>& drainage_area,
                               const xtensor_t<C>& cell_area,
                               const xtensor_t<S>& stack,
                               const xtensor_t<R>& receivers)
    {
        detail::compute_drainage_area_impl(drainage_area.derived_cast(),
                                           cell_area.derived_cast(),
                                           stack.derived_cast(),
                                           receivers.derived_cast());
    }


    /**
     * Compute drainage area for each node on a 2-d rectangular grid.
     *
     * @param drainage_area : ``[intent=inout, shape=(nrows, ncols)]``
     *     Drainage area at grid node.
     * @param stack :``[intent=in, shape=(nnodes)]``
     *     Stack position at grid node.
     * @param receivers : ``[intent=in, shape=(nnodes)]``
     *     Index of flow receiver at grid node.
     * @param dx : ``[intent=in]``
     *     Grid spacing in x.
     * @param dy : ``[intent=in]``
     *     Grid spacing in y.
     */
    template <class D, class S, class R>
    void compute_drainage_area(xtensor_t<D>& drainage_area,
                               const xtensor_t<S>& stack,
                               const xtensor_t<R>& receivers,
                               double dx,
                               double dy)
    {
        xt::xtensor<double, 1> cell_area = { dx * dy };
        compute_drainage_area(drainage_area, cell_area, stack, receivers);
    }
}  // namespace fastscapelib

#endif
