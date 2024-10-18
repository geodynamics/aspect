#ifndef FASTSCAPELIB_FLOW_IMPL_FLOW_GRAPH_INL_HPP
#define FASTSCAPELIB_FLOW_IMPL_FLOW_GRAPH_INL_HPP

#include <map>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <string>
#include <vector>

#include "fastscapelib/utils/containers.hpp"
#include "xtensor/xstrided_view.hpp"

#include "fastscapelib/flow/flow_graph_impl.hpp"
#include "fastscapelib/flow/flow_operator.hpp"


namespace fastscapelib
{
    template <class G, class S, class Tag>
    flow_graph<G, S, Tag>::flow_graph(G& grid, operators_type operators)
        : m_grid(grid)
        , m_operators(std::move(operators))
        , m_thread_pool(10)
    {
        m_impl_ptr = std::make_shared<impl_type>(grid, m_operators.all_single_flow());

        // sanity checks
        if (!m_operators.graph_updated())
        {
            throw std::invalid_argument(
                "must have at least one operator that updates the flow graph");
        }
        if (m_operators.out_flowdir() == flow_direction::undefined)
        {
            throw std::invalid_argument(
                "must have at least one operator that defines the output flow direction type");
        }

        // initialize default base levels at fixed value nodes
        m_impl_ptr->set_base_levels(m_grid.nodes_indices(node_status::fixed_value));

        // pre-allocate graph and elevation snapshots
        for (const auto& key : m_operators.graph_snapshot_keys())
        {
            bool single_flow = m_operators.snapshot_single_flow(key);
            auto graph = new self_type(grid, single_flow);

            m_graph_snapshots.insert({ key, std::unique_ptr<self_type>(std::move(graph)) });
            m_graph_impl_snapshots.insert({ key, (*m_graph_snapshots.at(key)).m_impl_ptr });
        }
        for (const auto& key : m_operators.elevation_snapshot_keys())
        {
            auto snapshot = data_array_type::from_shape(grid.shape());
            m_elevation_snapshots.insert(
                { key, std::make_unique<data_array_type>(std::move(snapshot)) });
        }

        // pre-allocate hydrologically corrected elevation
        if (m_operators.elevation_updated())
        {
            m_elevation_copy = xt::empty<data_type>(grid.shape());
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    const std::vector<const flow_operator*>& flow_graph<G, S, Tag>::operators() const
    {
        return m_operators.m_op_vec;
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    bool flow_graph<G, S, Tag>::single_flow() const
    {
        return m_operators.out_flowdir() == flow_direction::single;
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    const std::vector<std::string>& flow_graph<G, S, Tag>::graph_snapshot_keys() const
    {
        return m_operators.graph_snapshot_keys();
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    auto flow_graph<G, S, Tag>::graph_snapshot(std::string name) const -> self_type&
    {
        return *(m_graph_snapshots.at(name));
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    const std::vector<std::string>& flow_graph<G, S, Tag>::elevation_snapshot_keys() const
    {
        return m_operators.elevation_snapshot_keys();
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    auto flow_graph<G, S, Tag>::elevation_snapshot(std::string name) const -> const data_array_type&
    {
        return *(m_elevation_snapshots.at(name));
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    auto flow_graph<G, S, Tag>::update_routes(const data_array_type& elevation)
        -> const data_array_type&
    {
        if (!m_writeable)
        {
            throw std::runtime_error("cannot update routes (graph is read-only)");
        }

        data_array_type* elevation_ptr;

        if (m_operators.elevation_updated())
        {
            // reset and use hydrologically corrected elevation
            m_elevation_copy = elevation;
            elevation_ptr = &m_elevation_copy;
        }
        else
        {
            // pretty safe to remove the const qualifier (shouldn't be updated)
            elevation_ptr = const_cast<data_array_type*>(&elevation);
        }

        // loop over flow operator implementations
        for (auto op = m_operators.impl_begin(); op != m_operators.impl_end(); ++op)
        {
            op->apply(*m_impl_ptr, *elevation_ptr, m_thread_pool);
            op->save(*m_impl_ptr, m_graph_impl_snapshots, *elevation_ptr, m_elevation_snapshots);
        }

        return *elevation_ptr;
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    auto flow_graph<G, S, Tag>::grid() const -> grid_type&
    {
        return m_grid;
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    auto flow_graph<G, S, Tag>::size() const -> size_type
    {
        return m_grid.size();
    }

    template <class G, class S, class Tag>
    auto flow_graph<G, S, Tag>::grid_shape() const -> shape_type
    {
        // grid shape may have a different type (e.g., from xtensor containers)
        auto shape = m_grid.shape();
        shape_type data_array_shape(shape.begin(), shape.end());
        return data_array_shape;
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    auto flow_graph<G, S, Tag>::impl() const -> const impl_type&
    {
        return *m_impl_ptr;
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    auto flow_graph<G, S, Tag>::impl_ptr() const -> std::shared_ptr<impl_type>
    {
        return m_impl_ptr;
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    auto flow_graph<G, S, Tag>::base_levels() const -> std::vector<size_type>
    {
        const auto& impl_levels = m_impl_ptr->base_levels();
        std::vector<size_type> indices(impl_levels.begin(), impl_levels.end());
        return indices;
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    template <class C>
    void flow_graph<G, S, Tag>::set_base_levels(C&& levels)
    {
        if (!m_writeable)
        {
            throw std::runtime_error("cannot set base levels (graph is read-only)");
        }

        m_impl_ptr->set_base_levels(levels);
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    auto flow_graph<G, S, Tag>::mask() const -> dynamic_shape_container_t<container_selector, bool>
    {
        return m_impl_ptr->mask();
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    template <class C>
    void flow_graph<G, S, Tag>::set_mask(C&& mask)
    {
        if (!m_writeable)
        {
            throw std::runtime_error("cannot set mask (graph is read-only)");
        }
        if (!xt::same_shape(mask.shape(), m_grid.shape()))
        {
            throw std::runtime_error("cannot set mask (shape mismatch with grid shape)");
        }

        m_impl_ptr->set_mask(std::forward<C>(mask));
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    auto flow_graph<G, S, Tag>::accumulate(const data_array_type& src) const -> data_array_type
    {
        return m_impl_ptr->accumulate(src);
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    void flow_graph<G, S, Tag>::accumulate(data_array_type& acc, const data_array_type& src) const
    {
        return m_impl_ptr->accumulate(acc, src);
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    auto flow_graph<G, S, Tag>::accumulate(data_type src) const -> data_array_type
    {
        return m_impl_ptr->accumulate(src);
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    void flow_graph<G, S, Tag>::accumulate(data_array_type& acc, data_type src) const
    {
        return m_impl_ptr->accumulate(acc, src);
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    auto flow_graph<G, S, Tag>::basins() -> data_array_size_type
    {
        data_array_size_type basins = data_array_size_type::from_shape(m_grid.shape());
        auto basins_flat = xt::flatten(basins);

        m_impl_ptr->compute_basins();
        basins_flat = m_impl_ptr->basins();

        return basins;
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    template <class FK, class FKD>
    int flow_graph<G, S, Tag>::apply_kernel_seq(FK& kernel, FKD& data)
    {
        using nodes_indices_type = typename impl_type::nodes_indices_type;
        const nodes_indices_type* indices;

        switch (kernel.apply_dir)
        {
            case flow_graph_traversal_dir::any:
                indices = &impl().storage_indices();
                break;
            case flow_graph_traversal_dir::breadth_upstream:
                indices = &impl().bfs_indices();
                break;
            case flow_graph_traversal_dir::depth_upstream:
                indices = &impl().dfs_indices();
                break;
            default:
                throw std::runtime_error(
                    "Unsupported kernel application order for sequential execution");
                break;
        }

        auto new_node_data = kernel.node_data_create();
        if (kernel.node_data_init)
            kernel.node_data_init(new_node_data, data.data);

        for (std::size_t i : *indices)
        {
            if (kernel.node_data_getter(i, data.data, new_node_data))
            {
                throw std::runtime_error("Invalid index encountered in node_data getter "
                                         "function\n"
                                         "Please check if you are using dynamic receivers count "
                                         "('max_receivers=-1') or adjust this setting in the "
                                         "'Kernel' "
                                         "specification");
            };
            kernel.func(new_node_data);
            kernel.node_data_setter(i, new_node_data, data.data);
        }

        kernel.node_data_free(new_node_data);
        return 0;
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    template <class FK, class FKD>
    int flow_graph<G, S, Tag>::apply_kernel_par(FK& kernel, FKD& data)
    {
        using nodes_indices_type = typename impl_type::nodes_indices_type;
        const nodes_indices_type *indices, *levels;

        switch (kernel.apply_dir)
        {
            case flow_graph_traversal_dir::any:
                indices = &impl().storage_indices();
                levels = &impl().any_order_levels();
                break;
            case flow_graph_traversal_dir::breadth_upstream:
                indices = &impl().bfs_indices();
                levels = &impl().bfs_levels();
                break;
            default:
                throw std::runtime_error(
                    "Unsupported kernel application order for parallel execution");
                break;
        }

        auto n_threads = kernel.n_threads;

        m_thread_pool.resume();
        m_thread_pool.resize(n_threads);

        std::vector<decltype(kernel.node_data_create())> node_data(n_threads);
        for (auto i = 0; i < n_threads; ++i)
        {
            node_data[i] = kernel.node_data_create();

            if (kernel.node_data_init)
                kernel.node_data_init(node_data[i], data.data);
        }

        auto run = [&kernel, &data, indices, node_data](
                       std::size_t runner, std::size_t start, std::size_t end)
        {
            for (auto i = start; i < end; ++i)
            {
                auto node_idx = (*indices)[i];
                auto n_data = node_data[runner];
                if (kernel.node_data_getter(node_idx, data.data, n_data))
                {
                    throw std::runtime_error(
                        "Invalid index encountered in node_data getter "
                        "function\n"
                        "Please check if you are using dynamic receivers count "
                        "('max_receivers=-1') or adjust this setting in the "
                        "'Kernel' "
                        "specification");
                };
                kernel.func(n_data);
                kernel.node_data_setter(node_idx, n_data, data.data);
            }
        };

        for (std::size_t i = 1; i < levels->size(); ++i)
        {
            const size_type first_idx = (*levels)[i - 1];
            const size_type after_last_idx = (*levels)[i];
            const size_type level_size = after_last_idx - first_idx;

            if (level_size < kernel.min_level_size)
                run(0, first_idx, after_last_idx);
            else
                m_thread_pool.run_blocks(first_idx, after_last_idx, run, kernel.min_block_size);
        }

        for (std::size_t i = 0; i < n_threads; ++i)
            kernel.node_data_free(node_data[i]);

        m_thread_pool.pause();

        return 0;
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    template <class FK, class FKD>
    int flow_graph<G, S, Tag>::apply_kernel(FK& kernel, FKD& data)
    {
        int ret;
        if (kernel.n_threads > 1)
            ret = apply_kernel_par(kernel, data);
        else
            ret = apply_kernel_seq(kernel, data);

        return ret;
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    flow_graph<G, S, Tag>::flow_graph(grid_type& grid, bool single_flow)
        : m_writeable(false)
        , m_grid(grid)
        , m_thread_pool(10)
    {
        m_impl_ptr = std::make_shared<impl_type>(grid, single_flow);
    }
}  // namespace fastscapelib

#endif  // FASTSCAPELIB_FLOW_IMPL_FLOW_GRAPH_INL_HPP
