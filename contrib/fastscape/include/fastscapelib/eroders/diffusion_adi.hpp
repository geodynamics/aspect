/**
 * Functions to compute hillslope erosion.
 */
#ifndef FASTSCAPELIB_ERODERS_DIFFUSION_ADI_H
#define FASTSCAPELIB_ERODERS_DIFFUSION_ADI_H

#include <array>
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include "xtensor/xbroadcast.hpp"
#include "xtensor/xbuilder.hpp"
#include "xtensor/xnoalias.hpp"
#include "xtensor/xview.hpp"
#include "xtensor/xmanipulation.hpp"

#include "fastscapelib/grid/structured_grid.hpp"
#include "fastscapelib/utils/utils.hpp"
#include "fastscapelib/utils/xtensor_utils.hpp"


namespace fastscapelib
{
    template <class G>
    struct is_raster_grid
    {
        static constexpr bool value = G::is_structured() && G::is_uniform() && G::xt_ndims() == 2;
    };


    /**
     * Hillslope erosion using linear diffusion.
     *
     * It numerically solves the diffusion equation using an Alternating
     * Direction Implicit (ADI) scheme.
     *
     * The equation is given by:
     *
     * @f[
     * \frac{\partial h}{\partial t} = K \nabla^2 h
     * @f]
     *
     * where \f$K\f$ is an erosion coefficient (diffusivity) and \f$\nabla^2
     * h\f$ is the local curvature of the topographic surface.
     *
     * This equation implies that the amount of sediment eroded is linearly
     * proportional to the local gradient of the topographic surface.
     *
     * \rst
     * .. warning::
     *    Only raster grids are supported.
     *
     * .. note::
     *    This eroder assumes Dirichlet boundary conditions at the border
     *    nodes of the raster grid.
     * \endrst
     *
     * @tparam FG The raster grid type.
     * @tparam S The xtensor container selector for data array members.
     *
     */
    template <class G, class S = typename G::xt_selector>
    class diffusion_adi_eroder
    {
    public:
        using grid_type = G;
        using xt_selector = S;

        using data_type = typename grid_type::grid_data_type;
        using data_array_type = xt_array_t<xt_selector, data_type>;
        using data_tensor_type = xt_tensor_t<xt_selector, data_type, grid_type::xt_ndims()>;

        /**
         * Create a new diffusion ADI eroder.
         *
         * @param grid The input raster grid.
         * @param k_coef The \f$K\f$ value (spatially uniform or variable).
         *
         * @tparam K Either a scalar float type or an array (xtensor expression) type.
         */
        template <class K>
        diffusion_adi_eroder(G& grid,
                             K&& k_coef,
                             typename std::enable_if_t<is_raster_grid<G>::value>* = 0)
            : m_grid(grid)
        {
            auto grid_shape = m_grid.shape();
            m_nrows = grid_shape[0];
            m_ncols = grid_shape[1];

            m_erosion.resize({ m_nrows, m_ncols });

            set_k_coef(k_coef);
        };

        /**
         * Return the \f$K\f$ value.
         */
        data_array_type k_coef()
        {
            if (m_k_coef_is_scalar)
            {
                return xt::ones<data_type>(m_grid.shape()) * m_k_coef_scalar;
            }
            else
            {
                return m_k_coef_array;
            }
        }

        /**
         * Set a spatially uniform, scalar \f$K\f$ value.
         */
        template <class T>
        void set_k_coef(T value, typename std::enable_if_t<std::is_floating_point<T>::value>* = 0)
        {
            m_k_coef_is_scalar = true;
            m_k_coef_scalar = value;
            m_k_coef_array.resize({ 0 });
            set_factors();
        };

        /**
         * Set a spatially variable, array \f$K\f$ value.
         *
         * The array expression or container must have the same shape than the
         * raster grid arrays.
         */
        template <class K>
        void set_k_coef(K&& value, typename std::enable_if_t<xt::is_xexpression<K>::value>* = 0)
        {
            if (!xt::same_shape(value.shape(), m_grid.shape()))
            {
                throw std::runtime_error("cannot set k_coef value: shape mismatch");
            }

            m_k_coef_is_scalar = false;
            m_k_coef_array = value;
            set_factors();
        };

        /**
         * Solve diffusion for one time step.
         *
         * @param elevation The elevation of surface topography at each grid node.
         * @param dt The duration of the time step.
         */
        const data_array_type& erode(const data_array_type& elevation, double dt);

    private:
        grid_type& m_grid;

        using size_type = typename grid_type::size_type;
        using shape_type = typename grid_type::shape_type;
        size_type m_nrows, m_ncols;

        data_array_type m_erosion;

        bool m_k_coef_is_scalar;
        data_type m_k_coef_scalar;
        data_tensor_type m_k_coef_array;

        using factors_type = xt::xtensor<data_type, 3>;
        factors_type m_factors_row;
        factors_type m_factors_col;

        void set_factors();

        // tridiagonal matrix system arrays
        xt::xtensor<data_type, 1> m_vec;
        xt::xtensor<data_type, 1> m_lower;
        xt::xtensor<data_type, 1> m_diag;
        xt::xtensor<data_type, 1> m_upper;

        void resize_tridiagonal(size_type size);
        xt::xtensor<data_type, 1> solve_tridiagonal();

        template <class E, class F1, class F2>
        auto solve_adi_row(E&& elevation,
                           F1&& factors_row,
                           F2&& factors_col,
                           size_type nrows,
                           size_type ncols,
                           double dt);
    };

    /*
     * Set constant ADI factors
     */
    template <class G, class S>
    void diffusion_adi_eroder<G, S>::set_factors()
    {
        auto spacing = m_grid.spacing();
        data_type dx = spacing[1];
        data_type dy = spacing[0];

        std::array<size_type, 3> factors_shape({ { 3, m_nrows, m_ncols } });
        data_type fr, fc;

        if (m_k_coef_is_scalar)
        {
            fr = m_k_coef_scalar * 0.5 / (dy * dy);
            fc = m_k_coef_scalar * 0.5 / (dx * dx);

            m_factors_row = xt::ones<data_type>(factors_shape) * fr;
            m_factors_col = xt::ones<data_type>(factors_shape) * fc;
        }
        else
        {
            fr = 0.25 / (dy * dy);
            fc = 0.25 / (dx * dx);

            m_factors_row = xt::empty<data_type>(factors_shape);
            m_factors_col = xt::empty<data_type>(factors_shape);

            const auto& k = m_k_coef_array;

            for (size_type r = 1; r < m_nrows - 1; ++r)
            {
                for (size_type c = 1; c < m_ncols - 1; ++c)
                {
                    m_factors_row(0, r, c) = fr * (k(r - 1, c) + k(r, c));
                    m_factors_row(1, r, c) = fr / 2 * (k(r - 1, c) + 2 * k(r, c) + k(r + 1, c));
                    m_factors_row(2, r, c) = fr * (k(r, c) + k(r + 1, c));

                    m_factors_col(0, r, c) = fc * (k(r, c - 1) + k(r, c));
                    m_factors_col(1, r, c) = fc / 2 * (k(r, c - 1) + 2 * k(r, c) + k(r, c + 1));
                    m_factors_col(2, r, c) = fc * (k(r, c) + k(r, c + 1));
                }
            }
        }
    }


    template <class G, class S>
    void diffusion_adi_eroder<G, S>::resize_tridiagonal(size_type size)
    {
        m_vec.resize({ size });
        m_upper.resize({ size });
        m_diag.resize({ size });
        m_lower.resize({ size });
    }


    /*
     * Solve tri-diagonal system of equations using Thomas' algorithm (TDMA).
     */
    template <class G, class S>
    auto diffusion_adi_eroder<G, S>::solve_tridiagonal() -> xt::xtensor<data_type, 1>
    {
        size_type n = m_vec.size();

        auto result = xt::empty_like(m_vec);
        auto gam = xt::empty_like(m_vec);

        if (m_diag(0) == 0)
        {
            throw std::runtime_error("division by zero while solving tri-diagonal system");
        }

        auto bet = m_diag(0);
        result(0) = m_vec(0) / bet;

        for (size_type i = 1; i < n; ++i)
        {
            gam(i) = m_upper(i - 1) / bet;
            bet = m_diag(i) - m_lower(i) * gam(i);

            if (bet == 0)
            {
                throw std::runtime_error("division by zero while solving tri-diagonal system");
            }

            result(i) = (m_vec(i) - m_lower(i) * result(i - 1)) / bet;
        }

        for (int i = static_cast<int>(n) - 2; i > -1; --i)
        {
            result(i) -= gam(i + 1) * result(i + 1);
        }

        return result;
    }


    /*
     * Solve ADI linear diffusion for the row direction.
     */
    template <class G, class S>
    template <class E, class F1, class F2>
    auto diffusion_adi_eroder<G, S>::solve_adi_row(E&& elevation,
                                                   F1&& factors_row,
                                                   F2&& factors_col,
                                                   size_type nrows,
                                                   size_type ncols,
                                                   double dt)
    {
        xt::xtensor<double, 2> elevation_out = elevation;

        for (size_type r = 1; r < nrows - 1; ++r)
        {
            m_lower = -1 * xt::view(factors_col, 0, r, xt::all()) * dt;
            m_diag = 1 + 2 * xt::view(factors_col, 1, r, xt::all()) * dt;
            m_upper = -1 * xt::view(factors_col, 2, r, xt::all()) * dt;

            for (size_type c = 1; c < ncols - 1; ++c)
            {
                m_vec(c) = ((1 - 2 * factors_row(1, r, c) * dt) * elevation(r, c)
                            + factors_row(0, r, c) * elevation(r - 1, c) * dt
                            + factors_row(2, r, c) * elevation(r + 1, c) * dt);
            }

            // ignore grid boundary status and assume fixed value boundary conditions
            auto ilast = ncols - 1;

            m_lower(0) = 0;
            m_lower(ilast) = 0;
            m_diag(0) = 1;
            m_diag(ilast) = 1;
            m_upper(0) = 0;
            m_upper(ilast) = 0;
            m_vec(0) = elevation(r, 0);
            m_vec(ilast) = elevation(r, ilast);

            auto elevation_out_r = xt::view(elevation_out, r, xt::all());

            elevation_out_r = solve_tridiagonal();
        }

        return elevation_out;
    }


    template <class G, class S>
    auto diffusion_adi_eroder<G, S>::erode(const data_array_type& elevation, double dt)
        -> const data_array_type&
    {
        // solve for rows
        resize_tridiagonal(m_ncols);
        auto elevation_tmp
            = solve_adi_row(elevation, m_factors_row, m_factors_col, m_nrows, m_ncols, dt);

        // solve for cols (i.e., transpose)
        resize_tridiagonal(m_nrows);
        auto tranposed_dims = std::array<std::size_t, 3>{ 0, 2, 1 };

        auto elevation_next = solve_adi_row(xt::transpose(elevation_tmp),
                                            xt::transpose(m_factors_col, tranposed_dims),
                                            xt::transpose(m_factors_row, tranposed_dims),
                                            m_ncols,
                                            m_nrows,
                                            dt);

        auto erosion_v = xt::view(m_erosion, xt::all(), xt::all());
        erosion_v = elevation - xt::transpose(elevation_next);

        return m_erosion;
    }


    /**
     * Helper to create a new diffusion ADI eroder.
     *
     * @param grid The input raster grid.
     * @param k_coef The \f$K\f$ value (spatially uniform or variable).
     *
     * @tparam G The raster grid type.
     * @tparam K Either a scalar float type or an array (xtensor expression) type.
     */
    template <class G, class K>
    diffusion_adi_eroder<G> make_diffusion_adi_eroder(G& grid, K&& k_coef)
    {
        return diffusion_adi_eroder<G>(grid, k_coef);
    }

}  // namespace fastscapelib

#endif
