#ifndef FASTSCAPELIB_FLOW_KERNEL_HPP
#define FASTSCAPELIB_FLOW_KERNEL_HPP

#include <cstddef>
#include <cstdint>
#include <thread>
#include <iostream>
#include <vector>
#include <functional>

namespace fastscapelib
{
    /**
     * Flow graph traversal direction and order
     *
     * This enum class is used to specify the direction and order in which to
     * visit each node of a flow graph.
     *
     */
    enum class flow_graph_traversal_dir
    {
        any,                /**< unspecified direction */
        depth_downstream,   /**< from up to downstream in the depth-first order */
        depth_upstream,     /**< from down to upstream in the depth-first order */
        breadth_downstream, /**< from up to downstream in the breath-first order */
        breadth_upstream    /**< from down to upstream in the breath-first order */
    };

    /////////////////////////////////////////////////////////////////////////////////////////

    namespace detail
    {
        struct flow_kernel
        {
            std::function<int(void*)> func; /**< the kernel function to be applied on each node */
            std::function<int(std::size_t, void*, void*)>
                node_data_getter; /**< gets the node data from the kernel data at specified index
                                     before a kernel function call */
            std::function<int(std::size_t, void*, void*)>
                node_data_setter; /**< sets the kernel data back from the node data at specified
                                     index after a kernel function call */
            std::function<void*()> node_data_create; /**< creates a new node data */
            std::function<void(void*, void*)>
                node_data_init;                        /**< init a node data with kernel data*/
            std::function<void(void*)> node_data_free; /**< frees a node data */
            int n_threads;      /**< specifies the number of threads for parallel application of the
                                   kernel function */
            int min_block_size; /**< specifies the minimum block size (number of nodes) to dispatch
                                   to each thread */
            int min_level_size; /**< specifies the minimum level size to trigger parallel dispatch
                                   instead of sequential one */
            flow_graph_traversal_dir apply_dir
                = flow_graph_traversal_dir::any; /**<  order for kernel application*/
        };

        /////////////////////////////////////////////////////////////////////////////////////////

        struct flow_kernel_data
        {
            void* data; /**<  kernel data*/
        };
    }
}  // namespace fastscapelib

#endif  // FASTSCAPELIB_FLOW_KERNEL_HPP
