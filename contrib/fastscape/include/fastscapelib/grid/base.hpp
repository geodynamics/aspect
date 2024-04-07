/**
 * Common grid element types (node status, neighbors).
 * grid base (abstract) class.
 */
#ifndef FASTSCAPELIB_GRID_BASE_H
#define FASTSCAPELIB_GRID_BASE_H

#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <map>
#include <vector>

#include "xtensor/xadapt.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xview.hpp"

#include "fastscapelib/utils/iterators.hpp"
#include "fastscapelib/utils/xtensor_utils.hpp"


namespace fastscapelib
{

    //*****************
    //* Grid boundaries
    //*****************

    /**
     * Node status values.
     *
     * Node status is a label assigned to each grid/mesh node. It is used to
     * define the structure and boundaries of the modeled domain.
     *
     * The description of each value is given for indicative pruposes only.
     * Exact semantics may differ depending on the grid type or the objects that
     * are consuming grids.
     *
     * For example, ``looped`` is used only for structured grids.
     *
     */
    enum class node_status : std::uint8_t
    {
        core = 0,           /**< Inner grid node */
        fixed_value = 1,    /**< Dirichlet boundary condition */
        fixed_gradient = 2, /**< Neumann boundary condition */
        looped = 3          /**< Reflective boundaries */
    };

    namespace detail
    {

        inline bool node_status_cmp(node_status a, node_status b)
        {
            static std::map<node_status, int> priority{ { node_status::core, 0 },
                                                        { node_status::looped, 1 },
                                                        { node_status::fixed_gradient, 2 },
                                                        { node_status::fixed_value, 3 } };

            return priority[a] < priority[b];
        }
    }  // namespace detail


    //***************
    //* Grid elements
    //***************

    /**
     * Represents a grid/mesh node.
     */
    struct node
    {
        std::size_t idx;    /**< Node index */
        node_status status; /**< Node status */
    };

    /**
     * Represents a grid/mesh node neighbor.
     */
    struct neighbor
    {
        std::size_t idx;    /**< Index of the neighbor node */
        double distance;    /**< Distance to the neighbor node */
        node_status status; /**< Status at the neighbor node */

        bool operator==(const neighbor& rhs) const
        {
            return (idx == rhs.idx) && (distance == rhs.distance) && (status == rhs.status);
        }
    };

    namespace detail
    {
        /**
         * Add a given offset (positive, negative or zero) to a grid node index.
         *
         * This utility function is mainly to avoid signed/unsigned conversion
         * warnings. Use it carefully.
         */
        inline std::size_t add_offset(std::size_t idx, std::ptrdiff_t offset)
        {
            return static_cast<std::size_t>(static_cast<std::ptrdiff_t>(idx) + offset);
        }
    }

    /**********************************
     * Grid neighbors indices caching *
     **********************************/

    /**
     * Provides a cache for grid neighbor indices look-up.
     *
     * @tparam N The size of each array entry in the cache (should correspond to
     * the fixed maximum number of neighbors).
     */
    template <std::uint8_t N>
    class neighbors_cache
    {
    public:
        static constexpr unsigned int cache_width = N;

        template <class T>
        using storage_type = std::array<T, N>;

        using neighbors_indices_type = storage_type<std::size_t>;

        neighbors_cache(std::size_t size)
            : m_cache(cache_shape_type({ size }))
        {
            for (std::size_t i = 0; i < size; ++i)
            {
                m_cache[i].fill(std::numeric_limits<std::size_t>::max());
            }
        }

        bool has(const std::size_t& idx) const
        {
            return m_cache[idx][0] == std::numeric_limits<std::size_t>::max() ? false : true;
        }

        neighbors_indices_type& get(const std::size_t& idx)
        {
            return m_cache[idx];
        }

        neighbors_indices_type& get_storage(const std::size_t& idx)
        {
            return m_cache[idx];
        }

        void store(const std::size_t& idx, const neighbors_indices_type& neighbors_indices)
        {
            m_cache[idx] = neighbors_indices;
        }

        std::size_t cache_size() const
        {
            return m_cache.size();
        };

        std::size_t cache_used() const
        {
            std::size_t count = 0;

            for (std::size_t i = 0; i < m_cache.size(); ++i)
            {
                if (m_cache[i][0] != std::numeric_limits<std::size_t>::max())
                {
                    count += 1;
                }
            }

            return count;
        };

        void reset()
        {
            for (std::size_t i = 0; i < m_cache.size(); ++i)
            {
                m_cache[i].fill(std::numeric_limits<std::size_t>::max());
            }
        }

        void remove(const std::size_t& idx)
        {
            m_cache[idx].fill(std::numeric_limits<std::size_t>::max());
        }

    protected:
        using cache_type = xt::xtensor<neighbors_indices_type, 1>;
        using cache_shape_type = typename cache_type::shape_type;

        cache_type m_cache;
    };


    /**
     * A simple pass-through for grid neighbor indices look-up.
     *
     * @tparam N The size of temporary storage of indices (should correspond to
     * the fixed maximum number of neighbors). If set to zero, the cache will
     * use a variable size container (``std::vector``) storage type.
     */
    template <unsigned int N>
    class neighbors_no_cache
    {
    public:
        static constexpr unsigned int cache_width = N;

        template <class T>
        using storage_type = std::conditional_t<N == 0, std::vector<T>, std::array<T, N>>;

        using neighbors_indices_type = storage_type<std::size_t>;

        neighbors_no_cache(std::size_t /*size*/)
        {
        }

        bool has(const std::size_t& /*idx*/) const
        {
            return false;
        }

        neighbors_indices_type& get(const std::size_t& /*idx*/)
        {
            return m_node_neighbors;
        }

        neighbors_indices_type& get_storage(const std::size_t& /*idx*/)
        {
            return m_node_neighbors;
        }

        void store(const std::size_t& /*idx*/, const neighbors_indices_type neighbors_indices)
        {
            m_node_neighbors = neighbors_indices;
        }

        std::size_t cache_size() const
        {
            return 0;
        }

        std::size_t cache_used() const
        {
            return 0;
        }

        void reset()
        {
        }

        void remove(const std::size_t& /*idx*/)
        {
        }

    protected:
        neighbors_indices_type m_node_neighbors;
    };

    //****************
    //* Grid interface
    //****************

    // clang-format off
    /**
     * Small template class that holds a few types (static members) specialized
     * for each grid type.
     *
     * It is used for accessing those types from within the grid base classes.
     *
     * \rst
     * .. seealso::
     *   :cpp:class:`~template\<class S, class C> fastscapelib::grid_inner_types\<profile_grid_xt\<S, C>>`,
     *   :cpp:class:`~template\<class S, raster_connect RC, class C> fastscapelib::grid_inner_types\<raster_grid_xt\<S, RC, C>>`,
     *   :cpp:class:`~template\<class S, unsigned int N> fastscapelib::grid_inner_types\<trimesh_xt\<S, N>>`
     * \endrst
     */
    // clang-format on
    template <class G>
    struct grid_inner_types
    {
    };

    /**
     * Base class for all grid or mesh types.
     *
     * This class defines grid properties as well as a common API for iterating
     * through grid nodes and their neighbors.
     *
     * @tparam G The derived grid type.
     */
    template <class G>
    class grid
    {
    public:
        using derived_grid_type = G;
        using inner_types = grid_inner_types<derived_grid_type>;

        static constexpr bool is_structured();
        static constexpr bool is_uniform();
        static constexpr std::size_t xt_ndims();
        static constexpr std::uint8_t n_neighbors_max();

        using grid_data_type = typename inner_types::grid_data_type;
        using xt_selector = typename inner_types::xt_selector;
        using xt_type = xt_tensor_t<xt_selector, grid_data_type, inner_types::xt_ndims>;

        using size_type = typename xt_type::size_type;
        using shape_type = typename xt_type::shape_type;

        using neighbors_cache_type = typename inner_types::neighbors_cache_type;

        static_assert(neighbors_cache_type::cache_width == 0
                          || neighbors_cache_type::cache_width >= n_neighbors_max(),
                      "Cache width is too small!");

        using neighbors_type = std::vector<neighbor>;

        // using xt:xtensor for indices as not all containers support resizing
        // (e.g., using pyarray may cause segmentation faults with Python)
        using neighbors_indices_type = xt::xtensor<size_type, 1>;
        using neighbors_distances_type = xt::xtensor<grid_data_type, 1>;

        using nodes_status_type = xt_tensor_t<xt_selector, node_status, inner_types::xt_ndims>;

        size_type size() const noexcept;
        shape_type shape() const noexcept;

        const nodes_status_type& nodes_status() const;
        node_status nodes_status(const size_type& idx) const;

        inline grid_nodes_indices<G> nodes_indices() const;
        inline grid_nodes_indices<G> nodes_indices(node_status status) const;

        xt_type nodes_areas() const;
        grid_data_type nodes_areas(const size_type& idx) const noexcept;

        size_type neighbors_count(const size_type& idx) const;
        neighbors_distances_type neighbors_distances(const size_type& idx) const;

        // no const since it may update the cache internally
        neighbors_indices_type neighbors_indices(const size_type& idx);
        neighbors_indices_type& neighbors_indices(const size_type& idx,
                                                  neighbors_indices_type& neighbors_indices);

        neighbors_type neighbors(const size_type& idx);
        neighbors_type& neighbors(const size_type& idx, neighbors_type& neighbors);

        const neighbors_cache_type& neighbors_indices_cache();

    protected:
        using neighbors_indices_impl_type =
            typename neighbors_cache_type::template storage_type<size_type>;
        using neighbors_distances_impl_type =
            typename neighbors_cache_type::template storage_type<grid_data_type>;

        grid(std::size_t size)
            : m_neighbors_indices_cache(neighbors_cache_type(size)){};
        ~grid() = default;

        const derived_grid_type& derived_grid() const noexcept;
        derived_grid_type& derived_grid() noexcept;

        inline xt_type nodes_areas_impl() const;
        inline grid_data_type nodes_areas_impl(const size_type& idx) const noexcept;

        inline size_type neighbors_count_impl(const size_type& idx) const;

        void neighbors_indices_impl(neighbors_indices_impl_type& neighbors,
                                    const size_type& idx) const;

        inline const neighbors_indices_impl_type& get_nb_indices_from_cache(const size_type& idx);

        inline const neighbors_distances_impl_type& neighbors_distances_impl(
            const size_type& idx) const;

        neighbors_cache_type m_neighbors_indices_cache;
    };

    /**
     * @name Grid static properties
     */
    //@{
    /**
     * True if the grid is strutured or False if the grid is an unstructured mesh.
     */
    template <class G>
    constexpr bool grid<G>::is_structured()
    {
        return inner_types::is_structured;
    }

    /**
     * True if the grid is uniform (i.e., constant spacing along each
     * dimension).
     */
    template <class G>
    constexpr bool grid<G>::is_uniform()
    {
        return inner_types::is_uniform;
    }

    /**
     * Number of dimensions of the grid field arrays.
     */
    template <class G>
    constexpr std::size_t grid<G>::xt_ndims()
    {
        return inner_types::xt_ndims;
    }

    /**
     * Maximum number of grid node neighbors.
     *
     * For strutured grids this corresponds to the actual (fixed) number of
     * neighbors.
     */
    template <class G>
    constexpr std::uint8_t grid<G>::n_neighbors_max()
    {
        return inner_types::n_neighbors_max;
    }
    //@}

    template <class G>
    inline auto grid<G>::derived_grid() const noexcept -> const derived_grid_type&
    {
        return *static_cast<const derived_grid_type*>(this);
    }

    template <class G>
    inline auto grid<G>::derived_grid() noexcept -> derived_grid_type&
    {
        return *static_cast<derived_grid_type*>(this);
    }

    /**
     * @name Grid properties
     */
    //@{
    /**
     * Returns the total number of grid nodes.
     */
    template <class G>
    inline auto grid<G>::size() const noexcept -> size_type
    {
        return derived_grid().m_size;
    }

    /**
     * Returns the shape of the grid node arrays.
     */
    template <class G>
    inline auto grid<G>::shape() const noexcept -> shape_type
    {
        return derived_grid().m_shape;
    }

    /**
     * Returns the cache used to store node neighbor indices.
     */
    template <class G>
    auto grid<G>::neighbors_indices_cache() -> const neighbors_cache_type&
    {
        return m_neighbors_indices_cache;
    };
    //@}

    /**
     * @name Node methods
     */
    /**
     * Returns a virtual container that may be used to iterate over all grid
     * nodes.
     */
    template <class G>
    inline grid_nodes_indices<G> grid<G>::nodes_indices() const
    {
        const auto& derived = derived_grid();
        return grid_nodes_indices<G>(derived);
    };

    /**
     * Returns a virtual container that may be used to iterate over all grid
     * nodes of a given status.
     *
     * @param status The node status.
     */
    template <class G>
    inline grid_nodes_indices<G> grid<G>::nodes_indices(node_status status) const
    {
        const auto& derived = derived_grid();
        return grid_nodes_indices<G>(derived,
                                     [=](const grid& grid, size_type idx)
                                     { return grid.nodes_status().flat(idx) == status; });
    };

    /**
     * Returns a constant reference to the array of status at grid nodes.
     */
    template <class G>
    inline auto grid<G>::nodes_status() const -> const nodes_status_type&
    {
        return derived_grid().m_nodes_status;
    }

    /**
     * Returns the status at a given grid node.
     *
     * @param idx The grid node flat index.
     */
    template <class G>
    inline auto grid<G>::nodes_status(const size_type& idx) const -> node_status
    {
        return derived_grid().m_nodes_status.flat(idx);
    }

    /**
     * Returns the areas of the direct vicinity of each grid node as an array.
     *
     * Note: this creates a new container or returns a copy.
     */
    template <class G>
    inline auto grid<G>::nodes_areas() const -> xt_type
    {
        return std::move(nodes_areas_impl());
    }

    /**
     * Returns the area of the direct vicinity of a grid node.
     *
     * @param idx The grid node flat index.
     */
    template <class G>
    inline auto grid<G>::nodes_areas(const size_type& idx) const noexcept -> grid_data_type
    {
        return nodes_areas_impl(idx);
    }
    //@}

    /**
     * @name Neighbor methods
     */
    /**
     * Returns the number of neighbors of a given grid node.
     *
     * @param idx The grid node flat index.
     */
    template <class G>
    inline auto grid<G>::neighbors_count(const size_type& idx) const -> size_type
    {
        return neighbors_count_impl(idx);
    }

    /**
     * Returns an array of the indices of the neighbors of a given grid node.
     *
     * Follows looped boundary conditions, if any.
     *
     * @param idx The grid node flat index.
     */
    template <class G>
    inline auto grid<G>::neighbors_indices(const size_type& idx) -> neighbors_indices_type
    {
        neighbors_indices_type indices = xt::adapt(get_nb_indices_from_cache(idx));
        auto view = xt::view(indices, xt::range(0, neighbors_count(idx)));

        return view;
    }

    /**
     * Resize and fills an array with the neighbors indices of a given grid node.
     *
     * Follows looped boundary conditions, if any.
     *
     * This method prevents allocating a new container for better performance.
     *
     * @param idx The grid node flat index.
     * @param neighbors_indices Reference to the container to be updated with the neighbors indices.
     */
    template <class G>
    inline auto grid<G>::neighbors_indices(const size_type& idx,
                                           neighbors_indices_type& neighbors_indices)
        -> neighbors_indices_type&
    {
        const auto& n_count = neighbors_count(idx);
        const auto& n_indices = get_nb_indices_from_cache(idx);

        if (neighbors_indices.size() != n_count)
        {
            neighbors_indices.resize({ n_count });
        }

        for (size_type i = 0; i < n_count; ++i)
        {
            neighbors_indices[i] = n_indices[i];
        }

        return neighbors_indices;
    }

    /**
     * Returns an array of the distances to the neighbors of a given grid node.
     *
     * Follows looped boundary conditions, if any.
     *
     * @param idx The grid node flat index.
     */
    template <class G>
    inline auto grid<G>::neighbors_distances(const size_type& idx) const -> neighbors_distances_type
    {
        neighbors_distances_type distances = xt::adapt(neighbors_distances_impl(idx));
        auto view = xt::view(distances, xt::range(0, neighbors_count(idx)));

        return view;
    }

    /**
     * Returns a vector of the neighbors of a given grid node.
     *
     * Follows looped boundary conditions, if any.
     *
     * @param idx The grid node flat index.
     * @return A vector of neighbor node objects.
     */
    template <class G>
    inline auto grid<G>::neighbors(const size_type& idx) -> neighbors_type
    {
        neighbors_type nb;
        neighbors(idx, nb);

        return nb;
    }

    /**
     * Resize and fills a vactor with the neighbors of a given grid node.
     *
     * Follows looped boundary conditions, if any.
     *
     * This method prevents allocating a new container for better performance.
     *
     * @param idx The grid node flat index.
     * @param neighbors Reference to the vector to be updated with the neighbor objects.
     */
    template <class G>
    inline auto grid<G>::neighbors(const size_type& idx, neighbors_type& neighbors)
        -> neighbors_type&
    {
        size_type n_idx;
        const auto& n_count = neighbors_count(idx);
        const auto& n_indices = get_nb_indices_from_cache(idx);
        const auto& n_distances = neighbors_distances_impl(idx);

        if (neighbors.size() != n_count)
        {
            neighbors.resize({ n_count });
        }

        for (size_type i = 0; i < n_count; ++i)
        {
            n_idx = n_indices[i];
            neighbors[i] = neighbor({ n_idx, n_distances[i], nodes_status()[n_idx] });
        }

        return neighbors;
    }
    //@}

    template <class G>
    inline auto grid<G>::nodes_areas_impl() const -> xt_type
    {
        return derived_grid().nodes_areas_impl();
    }

    template <class G>
    inline auto grid<G>::nodes_areas_impl(const size_type& idx) const noexcept -> grid_data_type
    {
        return derived_grid().nodes_areas_impl(idx);
    }

    template <class G>
    inline auto grid<G>::neighbors_count_impl(const size_type& idx) const -> size_type
    {
        return derived_grid().neighbors_count_impl(idx);
    }

    template <class G>
    inline auto grid<G>::get_nb_indices_from_cache(const size_type& idx)
        -> const neighbors_indices_impl_type&
    {
        if (m_neighbors_indices_cache.has(idx))
        {
            neighbors_indices_impl_type& n_indices = m_neighbors_indices_cache.get(idx);
            return n_indices;
        }
        else
        {
            neighbors_indices_impl_type& n_indices = m_neighbors_indices_cache.get_storage(idx);
            this->derived_grid().neighbors_indices_impl(n_indices, idx);
            return n_indices;
        }
    }

    template <class G>
    inline auto grid<G>::neighbors_indices_impl(neighbors_indices_impl_type& neighbors,
                                                const size_type& idx) const -> void
    {
        return derived_grid().neighbors_indices_impl(neighbors, idx);
    }

    template <class G>
    inline auto grid<G>::neighbors_distances_impl(const size_type& idx) const
        -> const neighbors_distances_impl_type&
    {
        return derived_grid().neighbors_distances_impl(idx);
    }

}

#endif  // FASTSCAPELIB_GRID_BASE_H
