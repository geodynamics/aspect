#ifndef FASTSCAPELIB_FLOW_SNAPSHOT_H_
#define FASTSCAPELIB_FLOW_SNAPSHOT_H_

#include "fastscapelib/flow/flow_graph_impl.hpp"
#include "fastscapelib/flow/flow_operator.hpp"


namespace fastscapelib
{

    /**
     * Flow snapshot operator.
     *
     * A special flow operator used to save intermediate states
     * of the flow graph and/or topographic elevation values while
     * applying the other operators in chain.
     *
     * Those saved states are accessible from the ``flow_graph`` object.
     *
     */
    class flow_snapshot : public flow_operator
    {
    public:
        inline std::string name() const noexcept override
        {
            return "flow_snapshot";
        }

        flow_snapshot(std::string snapshot_name,
                      bool save_graph = true,
                      bool save_elevation = false)
            : m_snapshot_name(snapshot_name)
            , m_save_graph(save_graph)
            , m_save_elevation(save_elevation)
        {
        }

        const std::string& snapshot_name() const noexcept
        {
            return m_snapshot_name;
        }

        bool save_graph() const noexcept
        {
            return m_save_graph;
        }

        bool save_elevation() const noexcept
        {
            return m_save_elevation;
        }

    private:
        std::string m_snapshot_name;
        bool m_save_graph;
        bool m_save_elevation;
    };


    namespace detail
    {

        /*
         * Flow snapshot operator implementation (fixed array graph).
         */
        template <class FG>
        class flow_operator_impl<FG, flow_snapshot, flow_graph_fixed_array_tag>
            : public flow_operator_impl_base<FG, flow_snapshot>
        {
        public:
            using base_type = flow_operator_impl_base<FG, flow_snapshot>;
            using data_array_type = typename base_type::data_array_type;

            using graph_impl_map = std::map<std::string, FG&>;
            using elevation_map = std::map<std::string, std::unique_ptr<data_array_type>>;

            flow_operator_impl(std::shared_ptr<flow_snapshot> ptr)
                : base_type(std::move(ptr)){};

            void save(const FG& graph_impl,
                      graph_impl_map& graph_impl_snapshots,
                      const data_array_type& elevation,
                      elevation_map& elevation_snapshots) const
            {
                if (this->m_op_ptr->save_graph())
                {
                    _save(graph_impl, get_snapshot(graph_impl_snapshots));
                }
                if (this->m_op_ptr->save_elevation())
                {
                    _save(elevation, get_snapshot(elevation_snapshots));
                }
            }

        private:
            FG& get_snapshot(graph_impl_map& graph_impl_snapshots) const
            {
                return graph_impl_snapshots.at(this->m_op_ptr->snapshot_name());
            }

            data_array_type& get_snapshot(elevation_map& elevation_snapshots) const
            {
                return *(elevation_snapshots.at(this->m_op_ptr->snapshot_name()));
            }

            void _save(const FG& graph_impl, FG& graph_impl_snapshot) const
            {
                graph_impl_snapshot.m_receivers_count = graph_impl.m_receivers_count;
                graph_impl_snapshot.m_donors_count = graph_impl.m_donors_count;
                graph_impl_snapshot.m_dfs_indices = graph_impl.m_dfs_indices;

                if (graph_impl_snapshot.single_flow())
                {
                    auto receivers_col = xt::col(graph_impl_snapshot.m_receivers, 0);
                    receivers_col = xt::col(graph_impl.m_receivers, 0);
                    auto receivers_distance_col
                        = xt::col(graph_impl_snapshot.m_receivers_distance, 0);
                    receivers_distance_col = xt::col(graph_impl.m_receivers_distance, 0);
                    auto receivers_weight_col = xt::col(graph_impl_snapshot.m_receivers_weight, 0);
                    receivers_weight_col = xt::col(graph_impl.m_receivers_weight, 0);
                    auto donors_col = xt::col(graph_impl_snapshot.m_donors, 0);
                    donors_col = xt::col(graph_impl.m_donors, 0);
                }
                else
                {
                    graph_impl_snapshot.m_receivers = graph_impl.m_receivers;
                    graph_impl_snapshot.m_receivers_distance = graph_impl.m_receivers_distance;
                    graph_impl_snapshot.m_receivers_weight = graph_impl.m_receivers_weight;
                    graph_impl_snapshot.m_donors = graph_impl.m_donors;
                }
            }

            void _save(const data_array_type& elevation, data_array_type& elevation_snapshot) const
            {
                elevation_snapshot = elevation;
            }
        };
    }


    // Now we can implement flow_operator_sequence::update_snapshots
    template <class FG>
    void flow_operator_sequence<FG>::update_snapshots(const flow_snapshot& snapshot)
    {
        const auto& snapshot_name = snapshot.snapshot_name();

        if (snapshot.save_graph())
        {
            m_graph_snapshot_keys.push_back(snapshot_name);

            if (m_out_flowdir == flow_direction::undefined)
            {
                throw std::invalid_argument(
                    "no flow routing operator defined before graph snapshot.");
            }

            bool single_flow = m_out_flowdir == flow_direction::single ? true : false;
            m_graph_snapshot_single_flow.insert({ snapshot_name, single_flow });
        }
        if (snapshot.save_elevation())
        {
            m_elevation_snapshot_keys.push_back(snapshot_name);
        }
    }
}

#endif  // FASTSCAPELIB_FLOW_SNAPSHOT_H_
