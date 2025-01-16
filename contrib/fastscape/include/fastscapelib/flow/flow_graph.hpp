#ifndef FASTSCAPELIB_FLOW_FLOW_GRAPH_HPP
#define FASTSCAPELIB_FLOW_FLOW_GRAPH_HPP

#include "fastscapelib/flow/flow_graph_impl.hpp"
#include "fastscapelib/flow/flow_kernel.hpp"
#include "fastscapelib/flow/flow_operator.hpp"
#include "fastscapelib/utils/thread_pool.hpp"

#include "xtensor/xstrided_view.hpp"

#include <map>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <string>
#include <vector>


namespace fastscapelib
{

    /**
     * Main class used to compute or follow flow routes on
     * the topographic surface.
     *
     * @tparam G The grid type.
     * @tparam S The xtensor container selector for data array members.
     * @tparam Tag The flow graph implementation tag.
     */
    template <class G,
              class S = typename G::container_selector,
              class Tag = flow_graph_fixed_array_tag>
    class flow_graph
    {
    public:
        using self_type = flow_graph<G, S, Tag>;
        using grid_type = G;
        using container_selector = S;
        using impl_type = detail::flow_graph_impl<G, S, Tag>;
        using operators_type = flow_operator_sequence<impl_type>;

        using size_type = typename grid_type::size_type;

        using data_type = typename grid_type::grid_data_type;
        using data_array_type = dynamic_shape_container_t<container_selector, data_type>;
        using shape_type = typename data_array_type::shape_type;
        using data_array_size_type = dynamic_shape_container_t<container_selector, size_type>;

        using graph_map = std::map<std::string, std::unique_ptr<self_type>>;
        using graph_impl_map = std::map<std::string, std::shared_ptr<impl_type>>;
        using elevation_map = std::map<std::string, std::unique_ptr<data_array_type>>;

        /**
         * Flow graph constructor.
         *
         * @param grid The grid object.
         * @param operators The sequence of flow operators to use for updating the graph.
         */
        flow_graph(G& grid, operators_type operators);

        /**
         * Returns a vector of pointers to the flow operators objects used in
         * this graph.
         *
         * Note: the ``flow_operator`` base class has very limited polymorphism
         * (i.e., only its ``name`` method is virtual). It is still possible to
         * downcast the returned pointers (e.g., using the visitor pattern).
         */
        const std::vector<const flow_operator*>& operators() const;

        /**
         * Returns ``true`` if the graph has single flow directions.
         */
        bool single_flow() const;

        /*
         * Returns the names of the graph snapshots (if any).
         */
        const std::vector<std::string>& graph_snapshot_keys() const;

        /**
         * Graph snapshot getter.
         *
         * @param name Name of the snapshot.
         */
        self_type& graph_snapshot(std::string name) const;

        /**
         * Returns the names of the elevation snapshots.
         */
        const std::vector<std::string>& elevation_snapshot_keys() const;

        /**
         * Elevation snapshot getter.
         *
         * @param name Name of the snapshot.
         */
        const data_array_type& elevation_snapshot(std::string name) const;

        /**
         * Update flow routes from the input topographic surface.
         *
         * This applies in chain the flow operators and takes snapshots (if any).
         *
         * @param elevation The input topographic surface elevation.
         *
         * @return Either the input elevation unchanged or the elevation that
         * has been updated in order to route flow accross closed depressions.
         *
         */
        const data_array_type& update_routes(const data_array_type& elevation);

        /**
         * Returns a reference to the corresponding grid object.
         */
        grid_type& grid() const;

        /**
         * Returns the total number of nodes in the graph (equals to the size of
         * the grid).
         */
        size_type size() const;

        /**
         * Returns the shape of the grid node arrays.
         */
        shape_type grid_shape() const;

        /**
         * Returns a reference to the flow graph implementation.
         */
        const impl_type& impl() const;

        /**
         * Returns a shared pointer to the flow graph implementation.
         *
         * While this might be useful in some cases (e.g., Python bindings), in
         * general accessing the implementation through ``impl()`` (const
         * reference) should be preferred.
         */
        std::shared_ptr<impl_type> impl_ptr() const;

        /**
         * Return the indices of the base level nodes.
         */
        std::vector<size_type> base_levels() const;

        /**
         * Clear all existing base level nodes and set new ones.
         *
         * @tparam C Any stl-compatible container type.
         * @param levels The indices of the new base level nodes.
         */
        template <class C>
        void set_base_levels(C&& levels);

        /**
         * Return a mask of where elements with a value ``true`` correspond
         * to grid nodes that are not included in the flow graph.
         */
        dynamic_shape_container_t<container_selector, bool> mask() const;

        /**
         * Set a new grid mask.
         *
         * @tparam C Any xtensor-compatible container or expression type.
         * @param mask The new mask where elements with a value ``true``
         * correspond to grid nodes that are not included in the flow graph.
         */
        template <class C>
        void set_mask(C&& mask);

        /**
         * Traverse the flow graph in the top->down direction and accumulate
         * locally produced quantities or fluxes.
         *
         * The local quantitites or fluxes (i.e., source) may for example
         * correspond to precipitation, surface water runoff, sediment flux,
         * etc.
         *
         * The accumulated values represent at each node of the graph the
         * integral of the source over the node upslope contributing area.
         *
         * For example, if the source units are meters (height), the units of
         * the output accumulated values are cubic meters (volume).
         *
         * @param src The source, must be given per area unit and have the same
         * shape than ``grid_shape``.
         *
         * @return The output accumulated values (array of same shape than ``grid_shape``).
         */
        data_array_type accumulate(const data_array_type& src) const;

        /**
         * Traverse the flow graph in the top->down direction and accumulate
         * locally produced quantities or fluxes.
         *
         * This version allows reusing a pre-allocated array to store the output.
         *
         * @param acc The output accumulated values.
         * @param src The source, must be given per area unit.
         *
         * Both `acc` and `src` must of the same shape then ``grid_shape``.
         */
        void accumulate(data_array_type& acc, const data_array_type& src) const;

        /**
         * Traverse the flow graph in the top->down direction and accumulate
         * locally produced quantities or fluxes.
         *
         * @param src The spatially uniform source (scalar).
         *
         * @return The output accumulated values (array of same shape than ``grid_shape``).
         */
        data_array_type accumulate(data_type src) const;

        /**
         * Traverse the flow graph in the top->down direction and accumulate
         * locally produced quantities or fluxes.
         *
         * @param acc The output accumulated values.
         * @param src The spatially uniform source (scalar).
         *
         * Both `acc` and `src` must of the same shape then ``grid_shape``.
         */
        void accumulate(data_array_type& acc, data_type src) const;

        /**
         * Delineate catchments (or basins).
         *
         * A catchment is defined by all adjacent nodes that flow towards the
         * same outlet (or pit) graph node.
         *
         * @return Catchment ids (array of same shape than ``grid_shape``).
         *
         * @note
         * Results may be cached.
         * All masked grid nodes have the same assigned catchment id set by the maximum
         * limit of the integer value range. It is therefore preferable to mask the
         * results prior to, e.g., plotting it.
         */
        data_array_size_type basins();

        /**
         * Apply a given kernel along the flow graph.
         *
         * Visit the graph nodes in the direction and order given in the kernel
         * object, call the kernel function and fill the output variables
         * referenced in the kernel data.
         *
         * @tparam FK The flow kernel type.
         * @tparam FKD The flow kernel data type.
         * @param kernel The flow kernel object to apply along the graph.
         * @param data The object holding or referencing input and output data
         * used by the flow kernel.
         *
         */
        template <class FK, class FKD>
        int apply_kernel(FK& kernel, FKD& data);

    private:
        using thread_pool_type = thread_pool<size_type>;

        bool m_writeable = true;
        grid_type& m_grid;
        std::shared_ptr<impl_type> m_impl_ptr;
        data_array_type m_elevation_copy;

        graph_map m_graph_snapshots;
        graph_impl_map m_graph_impl_snapshots;
        elevation_map m_elevation_snapshots;

        operators_type m_operators;

        thread_pool_type m_thread_pool;

        // used internally for creating graph snapshots
        explicit flow_graph(grid_type& grid, bool single_flow);

        template <typename T,
                  typename F,
                  typename R = std::invoke_result_t<std::decay_t<F>, T, T, T>>
        void run_blocks(const T first_index, const T index_after_last, F&& func);

        template <class FK, class FKD>
        int apply_kernel_seq(FK& kernel, FKD& data);

        template <class FK, class FKD>
        int apply_kernel_par(FK& kernel, FKD& data);
    };
}

#include "./impl/flow_graph_inl.hpp"

#endif  // FASTSCAPELIB_FLOW_FLOW_GRAPH_HPP
