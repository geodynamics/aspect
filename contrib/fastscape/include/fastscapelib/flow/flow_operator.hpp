#ifndef FASTSCAPELIB_FLOW_OPERATOR_H_
#define FASTSCAPELIB_FLOW_OPERATOR_H_

#include <map>
#include <memory>
#include <string>
#include <vector>


namespace fastscapelib
{

    /**
     * Flow direction enum class.
     *
     * Used as a tag to specify the kind of flow graph (i.e., single, multiple
     * or undefined) expected as input or computed as output of a flow operator.
     *
     * A ``single`` direction flow (directed acyclic) graph is a tree (or a
     * forest of trees) where the total amount of matter (e.g., water, sediment)
     * at each grid node is propagated to a unique (downslope) node neighbor. In
     * a ``multiple`` direction flow graph, that amount is partionned among one
     * or more node neighbors.
     *
     * Use ``undefined`` for a flow operator that works on both single and
     * multiple flow direction graphs or that do not use the graph (e.g., sink
     * resolver only correcting the topographic elevation).
     *
     */
    enum class flow_direction
    {
        undefined,
        single,
        multi
    };

    /**
     * Flow operator.
     *
     * It represents a logical unit that can read and/or modify in-place the flow
     * graph and topographic elevation. It is a common type for flow routers, sink
     * resolvers and snapshots (save intermediate graph and/or elevation states).
     *
     * A flow operator may have one or more implementations, each relative to a
     * specific flow graph implementation. Note: this class and its derived classes
     * only contain the operators' (static) properties and parameters
     * (implementations are in separate classes).
     *
     * Do not use this class directly (it has no implementation). Use instead its
     * derived classes.
     *
     */
    class flow_operator
    {
    public:
        virtual ~flow_operator()
        {
        }

        /**
         * Returns the name of the operator.
         */
        inline virtual std::string name() const noexcept = 0;

        /**
         * Whether the operator updates topographic elevation.
         */
        static constexpr bool elevation_updated = false;

        /**
         * Whether the operator updates the flow graph.
         */
        static constexpr bool graph_updated = false;

        /**
         * The type of flow direction required as input of this operator.
         */
        static constexpr flow_direction in_flowdir = flow_direction::undefined;

        /**
         * The output flow direction type of this operator.
         */
        static constexpr flow_direction out_flowdir = flow_direction::undefined;
    };


    namespace detail
    {

        /**
         * Flow operator implementation base class.
         *
         * It exposes two methods:
         *
         * - ``apply``: may update in-place a flow graph implementation and/or
         *   topographic elevation
         * - ``save`` : only used by the flow_snapshot operator.
         *
         * By default, these two methods do nothing. A flow operator implementation
         * should re-implement only the ``apply`` method.
         *
         * The methods are not declared virtual (no polymorphism). Template
         * specialization + type erasure is used instead.
         *
         * A flow operator implementation is decoupled from its flow operator
         * corresponding instance (shared pointer). While flow operators may be
         * instantiated outside of the flow_graph class, flow operator implementations
         * are instantiated inside the flow_graph class as they need the
         * (grid-dependent) flow graph implementation and topographic elevation types.
         *
         * @tparam FG The flow graph implementation type (grid-dependent)
         * @tparam OP The flow operator type
         *
         */
        template <class FG, class OP>
        class flow_operator_impl_base
        {
        public:
            using graph_impl_type = FG;
            using data_array_type = typename graph_impl_type::data_array_type;
            using graph_impl_map = std::map<std::string, FG&>;
            using elevation_map = std::map<std::string, std::unique_ptr<data_array_type>>;

            void apply(FG& /*graph_impl*/, data_array_type& /*elevation*/)
            {
            }

            void save(const FG& /*graph_impl*/,
                      graph_impl_map& /*graph_impl_snapshots*/,
                      const data_array_type& /*elevation*/,
                      elevation_map& /*elevation_snapshots*/) const
            {
            }

        protected:
            flow_operator_impl_base(std::shared_ptr<OP> ptr)
                : m_op_ptr(std::move(ptr))
            {
            }

            ~flow_operator_impl_base() = default;

            std::shared_ptr<const OP> m_op_ptr;
        };

        /**
         * Flow operator implementation.
         *
         * This template class is not directly constructible. Template specialized
         * implementation classes must be provided for each operator (OP) for at least
         * one of the available flow graph implementations (Tag).
         *
         * @tparam FG The flow graph implementation type (grid-dependent)
         * @tparam OP The flow operator type
         * @tparam Tag The flow graph implementation tag
         *
         */
        template <class FG, class OP, class Tag>
        class flow_operator_impl : public flow_operator_impl_base<FG, OP>
        {
        public:
            flow_operator_impl(std::shared_ptr<OP> ptr) = delete;
        };

        /**
         * Flow operator implementation facade.
         *
         * This facade class implements type erasure so that multiple flow operators
         * of different types can be applied in chain in a flow graph.
         *
         * It takes a flow operator instance as input and creates a wrapper to the
         * corresponding flow operator implementation instance (if such implementation
         * exists for the flow graph implementation given as template argument).
         *
         * @tparam FG The flow graph implementation type (grid-dependent)
         *
         */
        template <class FG>
        class flow_operator_impl_facade
        {
        public:
            using data_array_type = typename FG::data_array_type;
            using graph_impl_map = std::map<std::string, FG&>;
            using elevation_map = std::map<std::string, std::unique_ptr<data_array_type>>;

            template <class OP>
            flow_operator_impl_facade(std::shared_ptr<OP> ptr)
                : m_wrapper_ptr(std::make_unique<flow_operator_impl_wrapper<OP>>(std::move(ptr)))
            {
            }

            // implement move semantics only (entity of flow_sequence)
            flow_operator_impl_facade(flow_operator_impl_facade<FG>& op_impl_facade) = delete;
            explicit flow_operator_impl_facade(flow_operator_impl_facade<FG>&& op_impl_facade)
                : m_wrapper_ptr(std::move(op_impl_facade.m_wrapper_ptr))
            {
            }

            void apply(FG& graph_impl, data_array_type& elevation)
            {
                return m_wrapper_ptr->apply(graph_impl, elevation);
            }

            void save(const FG& graph_impl,
                      graph_impl_map& graph_impl_snapshots,
                      const data_array_type& elevation,
                      elevation_map& elevation_snapshots) const
            {
                m_wrapper_ptr->save(
                    graph_impl, graph_impl_snapshots, elevation, elevation_snapshots);
            }

            struct flow_operator_impl_wrapper_base
            {
                virtual ~flow_operator_impl_wrapper_base()
                {
                }
                virtual void apply(FG& graph_impl, data_array_type& elevation) = 0;
                virtual void save(const FG& graph_impl,
                                  graph_impl_map& graph_impl_snapshots,
                                  const data_array_type& elevation,
                                  elevation_map& elevation_snapshots) const
                    = 0;
            };

            template <class OP>
            class flow_operator_impl_wrapper : public flow_operator_impl_wrapper_base
            {
            public:
                flow_operator_impl_wrapper(std::shared_ptr<OP> ptr)
                    : m_op_impl(std::move(ptr))
                {
                }

                void apply(FG& graph_impl, data_array_type& elevation) override
                {
                    return m_op_impl.apply(graph_impl, elevation);
                }

                void save(const FG& graph_impl,
                          graph_impl_map& graph_impl_snapshots,
                          const data_array_type& elevation,
                          elevation_map& elevation_snapshots) const override
                {
                    m_op_impl.save(
                        graph_impl, graph_impl_snapshots, elevation, elevation_snapshots);
                }

            private:
                flow_operator_impl<FG, OP, typename FG::tag> m_op_impl;
            };

            std::unique_ptr<flow_operator_impl_wrapper_base> m_wrapper_ptr;
        };
    }


    // forward delcarations
    class flow_snapshot;

    template <class G, class S, class Tag>
    class flow_graph;

    /**
     * Immutable container of flow operators (e.g., flow routers, sink resolvers,
     * flow snapshots) that are applied in chain when updating a flow graph.
     *
     * This class is not intended to be used as a stand-alone container. It is
     * used internally as an entity of flow_graph and can (should) be created
     * implicitly in the flow_graph constructor.
     *
     * @tparam FG The flow graph implementation type.
     */
    template <class FG>
    class flow_operator_sequence
    {
    public:
        using impl_type = FG;
        using operator_impl_type = detail::flow_operator_impl_facade<impl_type>;
        using iterator_type = typename std::vector<operator_impl_type>::iterator;

        template <class... OPs>
        flow_operator_sequence(OPs&&... operators)
        {
            int i = 0;
            (
                [&]
                {
                    ++i;
                    add_operator(std::forward<OPs>(operators));
                }(),
                ...);
        }

        // implement move semantics only (entity of flow_graph)
        flow_operator_sequence(flow_operator_sequence<FG>& operators) = delete;
        flow_operator_sequence(flow_operator_sequence<FG>&& operators)
            : m_op_vec(std::move(operators.m_op_vec))
            , m_op_impl_vec(std::move(operators.m_op_impl_vec))
            , m_graph_snapshot_keys(operators.graph_snapshot_keys())
            , m_graph_snapshot_single_flow(operators.m_graph_snapshot_single_flow)
            , m_elevation_snapshot_keys(operators.elevation_snapshot_keys())
            , m_elevation_updated(operators.elevation_updated())
            , m_graph_updated(operators.graph_updated())
            , m_out_flowdir(operators.out_flowdir())
            , m_all_single_flow(operators.all_single_flow())
        {
        }

        flow_operator_sequence<FG>& operator=(flow_operator_sequence<FG>& operators) = delete;
        flow_operator_sequence<FG>& operator=(flow_operator_sequence<FG>&& operators)
        {
            m_op_vec = std::move(operators.m_op_vec);
            m_op_impl_vec = std::move(operators.m_op_impl_vec);
            m_graph_snapshot_keys = std::move(operators.graph_snapshot_keys());
            m_graph_snapshot_single_flow = std::move(operators.m_graph_snapshot_single_flow);
            m_elevation_snapshot_keys = std::move(operators.elevation_snapshot_keys());
            m_elevation_updated = operators.elevation_updated();
            m_graph_updated = operators.graph_updated();
            m_out_flowdir = operators.out_flowdir();
            m_all_single_flow = operators.all_single_flow();
            return *this;
        }

        /**
         * STL-compatible iterator pointing to the 1st operator.
         */
        iterator_type begin()
        {
            return m_op_vec.begin();
        }

        /**
         * STL-compatible iterator pointing to the last operator.
         */
        iterator_type end()
        {
            return m_op_vec.end();
        }

        /**
         * Returns true if at least one operator in the sequence
         * updates topographic elevation.
         */
        bool elevation_updated() const
        {
            return m_elevation_updated;
        }

        /**
         * Returns true if at least one operator in the sequence
         * updates the flow graph.
         */
        bool graph_updated() const
        {
            return m_graph_updated;
        }

        /**
         * Returns the flow direction type (single, multiple or undefined)
         * of the final state of the flow graph, after having applied all operators.
         */
        flow_direction out_flowdir() const
        {
            return m_out_flowdir;
        }

        /**
         * Returns true if all intermediate states of the flow graph
         * have single flow directions.
         */
        bool all_single_flow() const
        {
            return m_all_single_flow;
        }

        /**
         * Returns true if the graph snapshot given by ``name`` is single flow.
         */
        bool snapshot_single_flow(std::string name)
        {
            return m_graph_snapshot_single_flow.at(name);
        }

        /**
         * Returns the names of all flow graph snapshots to create.
         */
        const std::vector<std::string>& graph_snapshot_keys() const
        {
            return m_graph_snapshot_keys;
        }

        /**
         * Returns the names of all topographic elevation snapshots to create.
         */
        const std::vector<std::string>& elevation_snapshot_keys() const
        {
            return m_elevation_snapshot_keys;
        }

    private:
        std::vector<const flow_operator*> m_op_vec;
        std::vector<operator_impl_type> m_op_impl_vec;

        std::vector<std::string> m_graph_snapshot_keys;
        std::map<std::string, bool> m_graph_snapshot_single_flow;
        std::vector<std::string> m_elevation_snapshot_keys;

        bool m_elevation_updated = false;
        bool m_graph_updated = false;
        flow_direction m_out_flowdir = flow_direction::undefined;
        bool m_all_single_flow = true;

        template <class OP>
        void add_operator(OP&& op)
        {
            auto op_ptr = std::make_shared<OP>(std::forward<OP>(op));
            add_operator(std::move(op_ptr));
        }

        template <class OP>
        void add_operator(std::shared_ptr<OP> ptr);

        void update_snapshots(const flow_snapshot& snapshot);

        // Only for bindings (in case the variadic templates constructor cannot be used).
        template <class _FG, class OPs>
        friend flow_operator_sequence<_FG> make_flow_operator_sequence(OPs&& operators);

        // iterate through operator implementations (in flow_graph)
        iterator_type impl_begin()
        {
            return m_op_impl_vec.begin();
        }

        iterator_type impl_end()
        {
            return m_op_impl_vec.end();
        }

        template <class G, class S, class Tag>
        friend class flow_graph;
    };


    // Add an operator of an abitrary type to the container, create an instance
    // of the corresponding implementation and update the sequence properties
    // accordingly.
    //
    // The current (final state) flow direction is updated only if the operator
    // updates the flow graph and explicitly defines an output flow direction
    // type (i.e., single or multiple).
    //
    // Also checks consistency between the current flow direction of the
    // sequence and the expected input flow direction of the operator to add.
    //
    template <class FG>
    template <class OP>
    void flow_operator_sequence<FG>::add_operator(std::shared_ptr<OP> ptr)
    {
        static_assert(std::is_base_of_v<flow_operator, OP>, "not a flow_operator type");

        if constexpr (std::is_same_v<OP, flow_snapshot>)
        {
            update_snapshots(*ptr);
        }
        if (ptr->in_flowdir != flow_direction::undefined && ptr->in_flowdir != m_out_flowdir)
        {
            throw std::invalid_argument("flow operator " + ptr->name()
                                        + " has incompatible input flow directions");
        }
        if (ptr->elevation_updated)
        {
            m_elevation_updated = true;
        }
        if (ptr->graph_updated)
        {
            m_graph_updated = true;

            if (ptr->out_flowdir != flow_direction::undefined)
            {
                m_out_flowdir = ptr->out_flowdir;

                if (ptr->out_flowdir != flow_direction::single)
                {
                    m_all_single_flow = false;
                }
            }
        }

        m_op_vec.push_back(ptr.get());
        m_op_impl_vec.push_back(operator_impl_type(std::move(ptr)));
    }
}


#endif  // FASTSCAPELIB_FLOW_OPERATOR_H_
