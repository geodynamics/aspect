#ifndef FASTSCAPELIB_FLOW_BASIN_GRAPH_H
#define FASTSCAPELIB_FLOW_BASIN_GRAPH_H

#pragma once

#include <algorithm>
#include <array>
#include <iostream>
#include <limits>
#include <numeric>
#include <queue>
#include <stack>
#include <stdexcept>
#include <vector>

#include "fastscapelib/grid/base.hpp"
#include "fastscapelib/flow/flow_graph_impl.hpp"
#include "fastscapelib/utils/utils.hpp"
#include "fastscapelib/utils/consts.hpp"
#include "fastscapelib/utils/union_find.hpp"


namespace fastscapelib
{

    /**
     * Algorithm used to compute the reduced graph (tree) of basins (minimum
     * spanning tree).
     *
     * Kruskal's algorithm is rather simple and has non-linear complexity \f$O(E
     * \log V)\f$ where \f$E\f$ is the number of basin connections and \f$V\f$
     * is the number of basins.
     *
     * Boruvka's algorithm is more advanced and has linear complexity for planar
     * graphs (most flow graphs are planar except, e.g., those built on raster
     * grids with connectivity including the diagonals).
     *
     * The gain in performance of Boruvka's algorithm over Kruskal's one is
     * significant only for (very) large graphs.
     */
    enum class mst_method
    {
        kruskal, /**< Kruskal's algorithm. */
        boruvka  /**< Boruvka's algorithm. */
    };

    inline std::ostream& operator<<(std::ostream& os, const mst_method meth)
    {
        switch (meth)
        {
            case mst_method::kruskal:
                os << "kruskal";
                break;
            case mst_method::boruvka:
                os << "boruvka";
                break;
        }
        return os;
    }

    namespace testing
    {
        // The corresponding test needs to access private members of basin_graph
        class basin_graph_orient_edges_Test;
    }


    /**
     * Represents a graph of adjacent basins (catchments) connected via a two
     * "pass" nodes where the water flow spills from one basin into another.
     *
     * \rst
     * Such graph is built and used by :cpp:class:`~fastscapelib::mst_sink_resolver`
     * to resolve flow routing accross closed depressions in the surface topography.
     * \endrst
     *
     * @tparam FG The flow graph implementation type (only
     * ``flow_graph_fixed_array_tag`` is supported).
     */
    template <class FG>
    class basin_graph
    {
    public:
        using flow_graph_impl_type = FG;

        static_assert(
            std::is_same<typename flow_graph_impl_type::tag, flow_graph_fixed_array_tag>::value,
            "basin graph requires the fixed array flow graph implementation");

        using size_type = typename flow_graph_impl_type::size_type;

        using data_type = typename flow_graph_impl_type::data_type;
        using data_array_type = typename flow_graph_impl_type::data_array_type;

        /**
         * Represents an edge of the graph of basins.
         */
        struct edge
        {
            size_type link[2]; /**< Indices of the two connected basins */
            size_type pass[2]; /**< Indices of the grid nodes forming the pass accross the basins */
            data_type pass_elevation; /**< Maximum elevation among the two pass nodes */
            data_type pass_length;    /**< Distance between the two pass nodes */

            /**
             * Create a new edge with no assigned pass yet.
             */
            static edge make_edge(const size_type& from, const size_type& to)
            {
                return edge{ { from, to },
                             { size_type(-1), size_type(-1) },
                             std::numeric_limits<data_type>::lowest(),
                             0. };
            }

            bool operator==(const edge& other)
            {
                return link[0] == other.link[0] && link[1] == other.link[1]
                       && pass[0] == other.pass[0] && pass[1] == other.pass[1]
                       && pass_elevation == other.pass_elevation
                       && pass_length == other.pass_length;
            }
        };

        /**
         * Create a new graph of basins from a flow graph.
         *
         * @param flow_graph_impl A flow graph implementation instance (fixed array).
         * @param basin_method The algorithm used to compute the reduced tree of basins.
         *
         * \rst
         * .. warning::
         *    A basin graph requires a flow graph with single direction flow paths as
         *    its current state. However, no validation is done here since the state of
         *    the same flow graph may later change from single to multiple direction
         *    after the basin graph has been (re)computed.
         * \endrst
         */
        basin_graph(const flow_graph_impl_type& flow_graph_impl, mst_method basin_method)
            : m_flow_graph_impl(flow_graph_impl)
            , m_mst_method(basin_method)
        {
            m_perf_boruvka = 0;
        }

        /**
         * Returns the total number of basins.
         */
        inline size_type basins_count() const
        {
            return m_flow_graph_impl.outlets().size();
        }

        /**
         * Returns a vector of node indices that represent the outlets of each
         * of the basins in this graph.
         */
        inline const std::vector<size_type>& outlets() const
        {
            return m_flow_graph_impl.outlets();
        }

        /**
         * Returns the edges of the graph of basins.
         */
        const std::vector<edge>& edges() const
        {
            return m_edges;
        }

        /**
         * Returns the reduced graph (tree) of basins after applying the
         * selected minimum spanning tree algorithm.
         */
        const std::vector<size_type>& tree() const
        {
            return m_tree;
        }

        /**
         * Update (re-compute from-scratch) the graph and its reduced tree from
         * the given topographic surface.
         *
         * @param elevation The topographic surface elevation at each grid node.
         */
        void update_routes(const data_array_type& elevation);

        /**
         * Boruvka's algorithm performance diagnostic.
         */
        size_t perf_boruvka() const
        {
            return m_perf_boruvka;
        }

    protected:
        void connect_basins(const data_array_type& elevation);
        void compute_tree_kruskal();
        void compute_tree_boruvka();
        void orient_edges();

    private:
        const flow_graph_impl_type& m_flow_graph_impl;

        mst_method m_mst_method;

        std::vector<size_type> m_outlets;  // bottom nodes of basins
        std::vector<edge> m_edges;
        std::vector<size_type> m_tree;  // indices of edges

        // root is used as a virtual basin graph node to which all outer basins are
        // connected, thus collecting the flow from the whole modeled domain. It
        // is used for the computation of the basin tree (minimum spanning
        // tree). As an implementation detail, this virtual node is actually
        // assigned to one of the outer basins, chosen arbitrarily.
        size_type m_root;

        // optimization for basin connections (arrays storing the positions of
        // already created edges)
        std::vector<size_type> m_edge_positions;
        std::vector<size_type> m_edge_positions_tmp;

        // kruskal
        std::vector<size_type> m_edges_indices;
        detail::union_find<size_type> m_basins_uf;

        // boruvka
        std::vector<std::array<size_type, 2>> m_link_basins;

        struct Connect
        {
            size_type begin;
            size_type size;
        };

        std::vector<Connect> m_adjacency;

        struct EdgeParse
        {
            size_type link_id;
            size_type next;
        };

        std::vector<EdgeParse> m_adjacency_list;
        std::vector<size_type> m_low_degrees;
        std::vector<size_type> m_large_degrees;
        std::vector<size_type> m_edge_bucket;
        std::vector<size_type> m_edge_in_bucket;

        // TODO: make it an option
        // 16 for 8-connectivity raster grid, 8 for plannar graph
        size_type m_max_low_degree = 16;

        template <class T>
        inline void check_capacity(std::vector<T> vec) const
        {
            ((void) (vec));
            assert(vec.size() <= vec.capacity());
        }

        size_t m_perf_boruvka;

        inline void increase_perf_boruvka()
        {
            ++m_perf_boruvka;
        }

        // reorder tree
        std::vector<size_type> m_nodes_connects_size;
        std::vector<size_type> m_nodes_connects_ptr;
        std::vector<size_type> m_nodes_adjacency;
        std::vector<std::tuple<size_type /* node */,
                               size_type /* parent */,
                               data_type /* pass elevation */,
                               data_type /* parent pass elevation */>>
            m_reorder_stack;

        // TODO: make it an option (true for sink route fill sloped)
        bool m_keep_order = false;

        std::vector<size_type> m_pass_stack;
        std::vector<size_type> m_parent_basins;

        friend class testing::basin_graph_orient_edges_Test;
    };


    template <class FG>
    void basin_graph<FG>::update_routes(const data_array_type& elevation)
    {
        connect_basins(elevation);

        if (m_mst_method == mst_method::kruskal)
        {
            compute_tree_kruskal();
        }
        else
        {
            compute_tree_boruvka();
        }

        orient_edges();
    }

    template <class FG>
    void basin_graph<FG>::connect_basins(const data_array_type& elevation)
    {
        using neighbors_type = typename flow_graph_impl_type::grid_type::neighbors_type;

        auto nbasins = basins_count();

        const auto& basins = m_flow_graph_impl.basins();
        const auto& receivers = m_flow_graph_impl.receivers();
        const auto& dfs_indices = m_flow_graph_impl.dfs_indices();

        auto& grid = m_flow_graph_impl.grid();

        // used to (re)initalize container index / position
        size_type init_idx = static_cast<size_type>(-1);

        neighbors_type neighbors;

        size_type ibasin;
        size_type current_basin = init_idx;

        // assume iteration starts at a base level node (outer basin)!
        bool is_inner_basin = false;

        // reset
        m_root = init_idx;
        m_edges.clear();
        m_edges.reserve(4 * nbasins);

        m_edge_positions.resize(nbasins);
        std::fill(m_edge_positions.begin(), m_edge_positions.end(), init_idx);
        m_edge_positions_tmp.reserve(nbasins);
        m_edge_positions_tmp.clear();

        for (const auto idfs : dfs_indices)
        {
            if (m_flow_graph_impl.is_masked(idfs))
            {
                continue;
            }

            const auto irec = receivers(idfs, 0);

            // any new basin visited (inner or outer)
            if (irec == idfs)
            {
                ibasin = basins(idfs);
                is_inner_basin = !m_flow_graph_impl.is_base_level(idfs);

                if (!is_inner_basin)
                {
                    if (m_root == init_idx)
                    {
                        // assign root to an existing outer basin
                        m_root = ibasin;
                    }
                    else
                    {
                        m_edges.push_back(edge::make_edge(m_root, ibasin));
                    }
                }
            }

            // any node in an inner basin
            if (is_inner_basin)
            {
                const data_type ielev = elevation.flat(idfs);

                for (auto n : grid.neighbors(idfs, neighbors))
                {
                    if (m_flow_graph_impl.is_masked(n.idx))
                    {
                        continue;
                    }

                    const size_type nbasin = basins(n.idx);

                    // skip if neighbor node is in the same basin or in an
                    // already connected adjacent basin unless the latter is an
                    // outer basin
                    bool skip = ibasin >= nbasin;
                    bool is_inner_nbasin = !m_flow_graph_impl.is_base_level(outlets()[nbasin]);
                    if (skip && is_inner_nbasin)
                    {
                        continue;
                    }

                    const data_type pass_elevation = std::max(ielev, elevation.flat(n.idx));

                    // just jumped from one basin to another
                    // -> update current basin and reset its visited neighbor basins
                    if (current_basin != ibasin)
                    {
                        for (const auto& ivisited : m_edge_positions_tmp)
                        {
                            m_edge_positions[ivisited] = init_idx;
                        }
                        m_edge_positions_tmp.clear();
                        current_basin = ibasin;
                    }

                    // try getting the position (index) of the edge connecting
                    // with the adjacent basin:
                    // - if it is undefined (-1), add a new edge
                    // - if it is defined, update the edge if a pass of lower elevation is found
                    const size_type edge_idx = m_edge_positions[nbasin];

                    if (edge_idx == init_idx)
                    {
                        m_edge_positions[nbasin] = m_edges.size();
                        m_edge_positions_tmp.push_back(nbasin);

                        m_edges.push_back(
                            { { ibasin, nbasin }, { idfs, n.idx }, pass_elevation, n.distance });
                    }
                    else if (pass_elevation < m_edges[edge_idx].pass_elevation)
                    {
                        m_edges[edge_idx] = edge{
                            { ibasin, nbasin }, { idfs, n.idx }, pass_elevation, n.distance
                        };
                    }
                }
            }
        }
    }

    template <class FG>
    void basin_graph<FG>::compute_tree_kruskal()
    {
        m_tree.reserve(basins_count() - 1);
        m_tree.clear();

        // sort edges indices by edge weight (pass elevation)
        m_edges_indices.resize(m_edges.size());
        std::iota(m_edges_indices.begin(), m_edges_indices.end(), 0);
        std::sort(m_edges_indices.begin(),
                  m_edges_indices.end(),
                  [&m_edges = m_edges](const size_type& i0, const size_type& i1)
                  { return m_edges[i0].pass_elevation < m_edges[i1].pass_elevation; });

        m_basins_uf.resize(basins_count());
        m_basins_uf.clear();

        for (size_type edge_idx : m_edges_indices)
        {
            size_type* link = m_edges[edge_idx].link;

            if (m_basins_uf.find(link[0]) != m_basins_uf.find(link[1]))
            {
                m_tree.push_back(edge_idx);
                m_basins_uf.merge(link[0], link[1]);
            }
        }
    }

    template <class FG>
    void basin_graph<FG>::compute_tree_boruvka()
    {
        const auto nbasins = basins_count();

        // used to (re)initalize container index / position
        size_type init_idx = static_cast<size_type>(-1);

        m_adjacency.clear();
        m_adjacency.resize(nbasins, { 0, 0 });
        m_low_degrees.reserve(nbasins);
        m_large_degrees.reserve(nbasins);

        m_edge_bucket.clear();
        m_edge_bucket.resize(nbasins, init_idx);

        // copy link basins
        m_link_basins.resize(m_edges.size());
        for (size_t i = 0; i < m_link_basins.size(); ++i)
        {
            m_link_basins[i][0] = m_edges[i].link[0];
            m_link_basins[i][1] = m_edges[i].link[1];
        }

        // first pass: create edge vector and compute adjacency size
        for (size_t lid = 0; lid < m_edges.size(); ++lid)
        {
            ++m_adjacency[m_link_basins[lid][0]].size;
            ++m_adjacency[m_link_basins[lid][1]].size;
        }

        // compute adjacency pointers
        m_adjacency[0].begin = 0;
        for (size_t nid = 1; nid < nbasins; ++nid)
        {
            m_adjacency[nid].begin = m_adjacency[nid - 1].begin + m_adjacency[nid - 1].size;
            m_adjacency[nid - 1].size = 0;
        }

        m_adjacency_list.resize(m_adjacency.back().begin + m_adjacency.back().size);
        m_adjacency.back().size = 0;

        for (size_t adj_data_i = 0; adj_data_i < m_adjacency_list.size(); ++adj_data_i)
            m_adjacency_list[adj_data_i].next = adj_data_i + 1;

        // second pass on edges: fill adjacency list
        for (size_t lid = 0; lid < m_edges.size(); ++lid)
        {
            auto& basins = m_link_basins[lid];

            m_adjacency_list[m_adjacency[basins[0]].begin + m_adjacency[basins[0]].size].link_id
                = lid;
            m_adjacency_list[m_adjacency[basins[1]].begin + m_adjacency[basins[1]].size].link_id
                = lid;

            ++m_adjacency[basins[0]].size;
            ++m_adjacency[basins[1]].size;
        }

        for (size_t nid = 0; nid < nbasins; ++nid)
        {
            check_capacity(m_low_degrees);
            check_capacity(m_large_degrees);
            // if degree is low enough
            if (m_adjacency[nid].size <= m_max_low_degree)
                m_low_degrees.push_back(nid);
            else
                m_large_degrees.push_back(nid);
        }

        m_perf_boruvka = 0;

        // compute the minimum spanning tree
        m_tree.reserve(nbasins - 1);
        m_tree.clear();

        while (m_low_degrees.size())
        {
            for (size_type nid : m_low_degrees)
            {
                // the node may have large degree after collapse
                if (m_adjacency[nid].size > m_max_low_degree)
                {
                    check_capacity(m_large_degrees);
                    m_large_degrees.push_back(nid);
                    continue;
                }

                // get the minimal weight edge that leaves that node
                size_type found_edge = init_idx;
                size_type node_B_id = init_idx;
                data_type found_edge_weight = std::numeric_limits<data_type>::max();

                size_type adjacency_data_ptr = m_adjacency[nid].begin;
                for (size_t step = 0; step < m_adjacency[nid].size; ++step)
                {
                    increase_perf_boruvka();

                    // find next adjacent edge in the list
                    size_type parsed_edge_id = m_adjacency_list[adjacency_data_ptr].link_id;
                    adjacency_data_ptr = m_adjacency_list[adjacency_data_ptr].next;

                    // check if the edge is valid (connected to a existing node)
                    // and if the weight is better than the previously found one
                    size_type opp_node = m_link_basins[parsed_edge_id][0];
                    if (opp_node == nid)
                        opp_node = m_link_basins[parsed_edge_id][1];

                    if (opp_node != nid && m_adjacency[opp_node].size > 0
                        && m_edges[parsed_edge_id].pass_elevation < found_edge_weight)
                    {
                        found_edge = parsed_edge_id;
                        found_edge_weight = m_edges[parsed_edge_id].pass_elevation;
                        node_B_id = opp_node;
                    }
                }

                if (found_edge == init_idx)
                    continue;  // TODO does it happens?

                // add edge to the tree
                check_capacity(m_tree);
                m_tree.push_back(found_edge);

                //  and collapse it toward opposite node

                // rename all A to B in adjacency of A
                adjacency_data_ptr = m_adjacency[nid].begin;
                for (size_t step = 0; step < m_adjacency[nid].size; ++step)
                {
                    increase_perf_boruvka();

                    // find next adjacent edge in the list
                    size_type edge_AC_id = m_adjacency_list[adjacency_data_ptr].link_id;

                    // TODO optimize that out?
                    if (step != m_adjacency[nid].size - 1)
                        adjacency_data_ptr = m_adjacency_list[adjacency_data_ptr].next;

                    // avoid self loop. A doesn't exist anymore, so edge AB
                    // will be discarded

                    if (m_link_basins[edge_AC_id][0] == nid)
                        m_link_basins[edge_AC_id][0] = node_B_id;
                    else
                        m_link_basins[edge_AC_id][1] = node_B_id;
                }

                // Append adjacency of B at the end of A
                m_adjacency_list[adjacency_data_ptr].next = m_adjacency[node_B_id].begin;

                // And collapse A into B
                m_adjacency[node_B_id].begin = m_adjacency[nid].begin;
                m_adjacency[node_B_id].size += m_adjacency[nid].size;

                // Remove the node from the graph
                m_adjacency[nid].size = 0;
            }

            m_low_degrees.clear();

            // Clean up graph (many edges are duplicates or self loops).
            size_type cur_large_degree = 0;
            for (size_type node_A_id : m_large_degrees)
            {
                // we will store all edges from A in the bucket, so that each edge
                // can appear only once
                m_edge_in_bucket.clear();
                size_type adjacency_data_ptr = m_adjacency[node_A_id].begin;

                for (size_t step = 0; step < m_adjacency[node_A_id].size; ++step)
                {
                    increase_perf_boruvka();

                    size_type edge_AB_id = m_adjacency_list[adjacency_data_ptr].link_id;
                    adjacency_data_ptr = m_adjacency_list[adjacency_data_ptr].next;

                    // find node B
                    size_type node_B_id = m_link_basins[edge_AB_id][0];
                    if (node_B_id == node_A_id)
                        node_B_id = m_link_basins[edge_AB_id][1];

                    if (m_adjacency[node_B_id].size > 0 && node_B_id != node_A_id)
                    {
                        // edge_bucket contain the edge_id connecting to opp_node_id
                        // or NodeId(-1)) if this is the first time we see it
                        size_type edge_AB_id_in_bucket = m_edge_bucket[node_B_id];

                        // first time we see
                        if (edge_AB_id_in_bucket == init_idx)
                        {
                            m_edge_bucket[node_B_id] = edge_AB_id;
                            m_edge_in_bucket.push_back(node_B_id);
                        }
                        else
                        {
                            // get weight of AB and of previously stored weight
                            data_type weight_in_bucket
                                = m_edges[edge_AB_id_in_bucket].pass_elevation;
                            data_type weight_AB = m_edges[edge_AB_id].pass_elevation;

                            // if both weight are the same, we choose edge
                            // with min id
                            if (weight_in_bucket == weight_AB)
                                m_edge_bucket[node_B_id]
                                    = std::min(edge_AB_id_in_bucket, edge_AB_id);
                            else if (weight_AB < weight_in_bucket)
                                m_edge_bucket[node_B_id] = edge_AB_id;
                        }
                    }
                }

                // recompute connectivity of node A
                size_type cur_ptr = m_adjacency[node_A_id].begin;
                m_adjacency[node_A_id].size = m_edge_in_bucket.size();

                for (size_type node_B_id : m_edge_in_bucket)
                {
                    increase_perf_boruvka();

                    m_adjacency_list[cur_ptr].link_id = m_edge_bucket[node_B_id];
                    cur_ptr = m_adjacency_list[cur_ptr].next;

                    // reset occupency of edge_bucket for latter use
                    m_edge_bucket[node_B_id] = init_idx;
                }


                // update low degree information, if node A has low degree
                if (m_adjacency[node_A_id].size <= m_max_low_degree)
                {
                    check_capacity(m_low_degrees);
                    // add the node in low degree list
                    if (m_adjacency[node_A_id].size > 0)
                        m_low_degrees.push_back(node_A_id);
                }
                else
                    m_large_degrees[cur_large_degree++] = node_A_id;
            }
            m_large_degrees.resize(cur_large_degree);
            check_capacity(m_large_degrees);
        }
    }

    /*
     * Orient the edges of the basin graph (tree) in the upward (or counter)
     * flow direction.
     *
     * If needed, swap the link (i.e., basin node indices) and pass (i.e., grid
     * node indices) of the edges.
     */
    template <class FG>
    void basin_graph<FG>::orient_edges()
    {
        const auto nbasins = basins_count();

        // used to (re)initalize container index / position
        size_type init_idx = static_cast<size_type>(-1);

        // nodes connections
        m_nodes_connects_size.resize(nbasins);
        std::fill(m_nodes_connects_size.begin(), m_nodes_connects_size.end(), size_t(0));
        m_nodes_connects_ptr.resize(nbasins);

        // parse the edges to compute the number of edges per node
        for (size_type l_id : m_tree)
        {
            m_nodes_connects_size[m_edges[l_id].link[0]] += 1;
            m_nodes_connects_size[m_edges[l_id].link[1]] += 1;
        }

        // compute the id of first edge in adjacency table
        m_nodes_connects_ptr[0] = 0;
        for (size_t i = 1; i < nbasins; ++i)
        {
            m_nodes_connects_ptr[i] = (m_nodes_connects_ptr[i - 1] + m_nodes_connects_size[i - 1]);
            m_nodes_connects_size[i - 1] = 0;
        }

        // create the adjacency table
        m_nodes_adjacency.resize(m_nodes_connects_ptr.back() + m_nodes_connects_size.back());
        m_nodes_connects_size.back() = 0;

        // parse the edges to update the adjacency
        for (size_type l_id : m_tree)
        {
            size_type n0 = m_edges[l_id].link[0];
            size_type n1 = m_edges[l_id].link[1];
            m_nodes_adjacency[m_nodes_connects_ptr[n0] + m_nodes_connects_size[n0]] = l_id;
            m_nodes_adjacency[m_nodes_connects_ptr[n1] + m_nodes_connects_size[n1]] = l_id;
            m_nodes_connects_size[n0] += 1;
            m_nodes_connects_size[n1] += 1;
        }

        // depth-first parse of the tree, starting from basin0
        // stack of node, parent
        m_reorder_stack.reserve(nbasins);
        m_reorder_stack.clear();

        m_reorder_stack.push_back({ m_root,
                                    m_root,
                                    std::numeric_limits<data_type>::min(),
                                    std::numeric_limits<data_type>::min() });

        if (m_keep_order)
        {
            m_pass_stack.clear();
            m_parent_basins.resize(nbasins);
            m_parent_basins[m_root] = m_root;
        }

        while (m_reorder_stack.size())
        {
            size_type node, parent;
            data_type pass_elevation, parent_pass_elevation;
            std::tie(node, parent, pass_elevation, parent_pass_elevation) = m_reorder_stack.back();
            m_reorder_stack.pop_back();


            for (size_t i = m_nodes_connects_ptr[node];
                 i < m_nodes_connects_ptr[node] + m_nodes_connects_size[node];
                 ++i)
            {
                edge& edg = m_edges[m_nodes_adjacency[i]];

                // the edge comming from the parent node has already been updated.
                // in this case, the edge is (parent, node)
                if (edg.link[0] == parent && node != parent)
                {
                    if (m_keep_order)
                    {
                        // force children of base nodes to be parsed
                        if (pass_elevation <= parent_pass_elevation && edg.pass[0] != init_idx)
                            // the pass is bellow the water level of the parent basin
                            m_parent_basins[edg.link[1]] = m_parent_basins[edg.link[0]];
                        else
                        {
                            m_parent_basins[edg.link[1]] = edg.link[1];
                            m_pass_stack.push_back(m_nodes_adjacency[i]);
                        }
                    }
                }
                else
                {
                    // we want the edge to be (node, next), where next is upper in flow order
                    // we check if the first node of the edge is not "node"
                    if (node != edg.link[0])
                    {
                        std::swap(edg.link[0], edg.link[1]);
                        std::swap(edg.pass[0], edg.pass[1]);
                    }

                    m_reorder_stack.push_back({ edg.link[1],
                                                node,
                                                std::max(edg.pass_elevation, pass_elevation),
                                                pass_elevation });
                }
            }
        }
    }
}

#endif
