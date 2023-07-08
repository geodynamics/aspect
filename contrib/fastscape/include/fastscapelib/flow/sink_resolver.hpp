/**
 * @brief Class used to resolve local depressions of
 * the topographic surface.
 *
 * ``sink_resolver`` is meant to be used in combination
 * with ``flow_router`` through a ``flow_graph``.
 *
 */
#ifndef FASTSCAPELIB_FLOW_SINK_RESOLVER_H
#define FASTSCAPELIB_FLOW_SINK_RESOLVER_H

#include <iostream>
#include <limits>
#include <memory>

#include "fastscapelib/algo/pflood.hpp"
#include "fastscapelib/flow/flow_graph_impl.hpp"
#include "fastscapelib/flow/flow_operator.hpp"
#include "fastscapelib/flow/basin_graph.hpp"


namespace fastscapelib
{

    /**
     * Priority-flood sink resolver operator.
     *
     * This flow operator fills the closed depressions in the topographic
     * surface using the priority flood algorithm +epsilon variant (Barnes et
     * al., 2014). This variant prevents flat surfaces and hence ensure
     * that the flow can be routed towards the outlets without disruption.
     *
     * \rst
     * See :cite:t:`Barnes2014` for a more detailed description of the algorithm.
     *
     * \endrst
     *
     */
    struct pflood_sink_resolver : public flow_operator
    {
        inline std::string name() const noexcept override
        {
            return "pflood_sink_resolver";
        }

        static constexpr bool elevation_updated = true;
    };


    namespace detail
    {

        /**
         * Priority-flood sink resolver operator implementation.
         */
        template <class FG, class Tag>
        class flow_operator_impl<FG, pflood_sink_resolver, Tag>
            : public flow_operator_impl_base<FG, pflood_sink_resolver>
        {
        public:
            using graph_impl_type = FG;
            using base_type = flow_operator_impl_base<graph_impl_type, pflood_sink_resolver>;

            using data_array_type = typename graph_impl_type::data_array_type;

            flow_operator_impl(std::shared_ptr<pflood_sink_resolver> ptr)
                : base_type(std::move(ptr)){};

            void apply(graph_impl_type& graph_impl, data_array_type& elevation)
            {
                fill_sinks_sloped(graph_impl.grid(), elevation);
            }
        };
    }


    /**
     * Method used by ``mst_sink_resolver`` to route flow within each closed
     * depressions.
     *
     * The ``basic`` method is the most efficient one but does not result in a
     * "realistic", planar flow graph.
     *
     * The ``carve`` method mimics carving a narrow canyon within the
     * depression.
     */
    enum class mst_route_method
    {
        basic, /**< Connect the pit node directly to the depression spill node. */
        carve  /**< Revert the (unique) flow path between the spill and pit nodes */
    };

    inline std::ostream& operator<<(std::ostream& os, const mst_route_method meth)
    {
        switch (meth)
        {
            case mst_route_method::basic:
                os << "basic";
                break;
            case mst_route_method::carve:
                os << "carve";
                break;
        }
        return os;
    }

    /**
     * Minimum Spanning Tree (MST) sink resolver operator.
     *
     * This flow operator re-routes the flow trapped in closed depressions
     * towards their spill, using an efficient algorithm that explicitly
     * computes a graph of (inner and outer) basins and reduces it as a tree
     * (Cordonnier et al., 2019).
     *
     * It requires a single flow graph as input.
     *
     * This operator also use the updated routes in closed depressions to
     * fill these with nearly flat surfaces (a tiny slope ensure natural
     * flow routing for the operators applied after this one).
     *
     * \rst
     * See :cite:t:`Cordonnier2019` for a more detailed description of the algorithm.
     *
     * \endrst
     *
     * @see fastscapelib::basin_graph
     *
     */
    class mst_sink_resolver : public flow_operator
    {
    public:
        inline std::string name() const noexcept override
        {
            return "mst_sink_resolver";
        }

        static constexpr bool graph_updated = true;
        static constexpr bool elevation_updated = true;
        static constexpr flow_direction in_flowdir = flow_direction::single;
        static constexpr flow_direction out_flowdir = flow_direction::single;

        mst_sink_resolver(mst_method basin_method = mst_method::kruskal,
                          mst_route_method route_method = mst_route_method::carve)
            : m_basin_method(basin_method)
            , m_route_method(route_method)
        {
        }

        mst_method m_basin_method = mst_method::kruskal;
        mst_route_method m_route_method = mst_route_method::carve;
    };

    inline std::ostream& operator<<(std::ostream& os, const mst_sink_resolver resolver)
    {
        os << resolver.m_basin_method << "-" << resolver.m_route_method;
        return os;
    }


    namespace detail
    {

        /**
         * Minimum Spanning Tree (MST) sink resolver operator implementation.
         */
        template <class FG>
        class flow_operator_impl<FG, mst_sink_resolver, flow_graph_fixed_array_tag>
            : public flow_operator_impl_base<FG, mst_sink_resolver>
        {
        public:
            using graph_impl_type = FG;
            using base_type = flow_operator_impl_base<graph_impl_type, mst_sink_resolver>;

            using size_type = typename graph_impl_type::size_type;
            using data_type = typename graph_impl_type::data_type;
            using data_array_type = typename graph_impl_type::data_array_type;

            using basin_graph_type = basin_graph<graph_impl_type>;

            flow_operator_impl(std::shared_ptr<mst_sink_resolver> ptr)
                : base_type(std::move(ptr)){};

            void apply(graph_impl_type& graph_impl, data_array_type& elevation)
            {
                // make sure the basins are up-to-date
                graph_impl.compute_basins();

                if (graph_impl.pits().empty())
                {
                    // return early, no sink to resolve
                    return;
                }

                get_basin_graph(graph_impl).update_routes(elevation);

                if (this->m_op_ptr->m_route_method == mst_route_method::basic)
                {
                    update_routes_sinks_basic(graph_impl, elevation);
                }
                else
                {
                    // carve
                    update_routes_sinks_carve(graph_impl);
                }

                // finalize flow route update (donors and dfs graph traversal indices)
                graph_impl.compute_donors();
                graph_impl.compute_dfs_indices_bottomup();

                // fill sinks with tiny tilted surface
                fill_sinks_sloped(graph_impl, elevation);
            };

        private:
            std::unique_ptr<basin_graph_type> m_basin_graph_ptr;

            // basin graph edges are oriented in the counter flow direction
            static constexpr std::uint8_t outflow = 0;
            static constexpr std::uint8_t inflow = 1;

            /*
             * Get the basin graph instance, create it if it doesn't exists.
             */
            basin_graph_type& get_basin_graph(const graph_impl_type& graph_impl)
            {
                if (!m_basin_graph_ptr)
                {
                    m_basin_graph_ptr = std::make_unique<basin_graph_type>(
                        graph_impl, this->m_op_ptr->m_basin_method);
                }

                return *m_basin_graph_ptr;
            }

            void update_routes_sinks_basic(graph_impl_type& graph_impl,
                                           const data_array_type& elevation);
            void update_routes_sinks_carve(graph_impl_type& graph_impl);

            void fill_sinks_sloped(graph_impl_type& graph_impl, data_array_type& elevation);
        };

        template <class FG>
        void flow_operator_impl<FG, mst_sink_resolver, flow_graph_fixed_array_tag>::
            update_routes_sinks_basic(graph_impl_type& graph_impl, const data_array_type& elevation)
        {
            const auto& basin_graph = get_basin_graph(graph_impl);
            auto& receivers = graph_impl.m_receivers;
            auto& dist2receivers = graph_impl.m_receivers_distance;

            // pits corresponds to outlets in inner basins
            const auto& pits = basin_graph.outlets();

            for (size_type edge_idx : basin_graph.tree())
            {
                auto& edge = basin_graph.edges()[edge_idx];

                // skip outer basins
                if (edge.pass[outflow] == static_cast<size_type>(-1))
                {
                    continue;
                }

                size_type pit_inflow = pits[edge.link[inflow]];

                // set infinite distance from the pit node to its receiver
                // (i.e., one of the two pass nodes).
                dist2receivers(pit_inflow, 0) = std::numeric_limits<data_type>::max();

                if (elevation.flat(edge.pass[inflow]) < elevation.flat(edge.pass[outflow]))
                {
                    receivers(pit_inflow, 0) = edge.pass[outflow];
                }
                else
                {
                    // also need to resolve the flow between the two
                    // nodes forming the pass, i.e., route the flow like this:
                    // pit -> pass (inflow) -> pass (outflow)
                    receivers(pit_inflow, 0) = edge.pass[inflow];
                    receivers(edge.pass[inflow], 0) = edge.pass[outflow];
                    dist2receivers(edge.pass[inflow], 0) = edge.pass_length;
                }
            }
        }

        template <class FG>
        void flow_operator_impl<FG, mst_sink_resolver, flow_graph_fixed_array_tag>::
            update_routes_sinks_carve(graph_impl_type& graph_impl)
        {
            const auto& basin_graph = get_basin_graph(graph_impl);
            auto& receivers = graph_impl.m_receivers;
            auto& dist2receivers = graph_impl.m_receivers_distance;

            // pits corresponds to outlets in inner basins
            const auto& pits = basin_graph.outlets();

            for (size_type edge_idx : basin_graph.tree())
            {
                auto& edge = basin_graph.edges()[edge_idx];

                // skip outer basins
                if (edge.pass[outflow] == static_cast<size_type>(-1))
                {
                    continue;
                }

                size_type pit_inflow = pits[edge.link[inflow]];

                // start at the pass (inflow) and then below follow the
                // receivers until the pit is found
                size_type current_node = edge.pass[inflow];
                size_type next_node = receivers(current_node, 0);
                data_type previous_dist = dist2receivers(current_node, 0);

                // re-route the pass inflow to the pass outflow
                receivers(current_node, 0) = edge.pass[outflow];
                dist2receivers(current_node, 0) = edge.pass_length;

                while (current_node != pit_inflow)
                {
                    auto rec_next_node = receivers(next_node, 0);
                    receivers(next_node, 0) = current_node;
                    std::swap(dist2receivers(next_node, 0), previous_dist);
                    current_node = next_node;
                    next_node = rec_next_node;
                }
            }
        }

        template <class FG>
        void
        flow_operator_impl<FG, mst_sink_resolver, flow_graph_fixed_array_tag>::fill_sinks_sloped(
            graph_impl_type& graph_impl, data_array_type& elevation)
        {
            const auto& dfs_indices = graph_impl.dfs_indices();
            const auto& receivers = graph_impl.receivers();

            for (const auto& idfs : dfs_indices)
            {
                const auto& irec = receivers(idfs, 0);

                if (idfs == irec)
                {
                    continue;
                }

                const auto& irec_elev = elevation.flat(irec);

                if (elevation.flat(idfs) <= elevation.flat(irec))
                {
                    auto tiny_step
                        = std::nextafter(irec_elev, std::numeric_limits<data_type>::infinity());
                    elevation.flat(idfs) = tiny_step;
                }
            }
        }
    }
}

#endif
