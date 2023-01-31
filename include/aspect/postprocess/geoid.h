/*
 Copyright (C) 2015 - 2022 by the authors of the ASPECT code.

 This file is part of ASPECT.

 ASPECT is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2, or (at your option)
 any later version.

 ASPECT is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with ASPECT; see the file LICENSE.  If not see
 <http://www.gnu.org/licenses/>.
 */


#ifndef _aspect_postprocess_geoid_h
#define _aspect_postprocess_geoid_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace Postprocess
  {
    /**
     * A postprocessor that computes the geoid anomaly at the surface.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class Geoid : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Evaluate the solution for the geoid in spherical harmonics and then transfer it to grid output.
         */
        std::pair<std::string,std::string>
        execute (TableHandler &statistics) override;

        /**
         * Register with the simulator the other postprocessors that we need
         * (namely: dynamic topography).
         */
        std::list<std::string>
        required_other_postprocessors() const override;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

        /**
         * Find if the top or bottom boundaries are free surfaces.
         */
        void initialize() override;

        /**
         * Evaluate the geoid solution at a point. The evaluation point
         * must be outside of the model domain, and it must be called
         * after execute().
         */
        double
        evaluate (const Point<dim> &) const;

      private:
        /**
         * Parameters to set the maximum and minimum degree when computing geoid from spherical harmonics
         */
        unsigned int max_degree;
        unsigned int min_degree;

        /**
         * A parameter to control whether to output the data in geographical coordinates.
         * If true, output the data in longitudes and latitudes; if false, output data in x y z.
         */
        bool output_in_lat_lon;

        /**
         * A parameter to control whether to output the spherical harmonic coefficients of the geoid anomaly
         */
        bool output_geoid_anomaly_SH_coes;

        /**
         * A parameter to control whether to output the spherical harmonic coefficients of the surface topography contribution
         */
        bool output_surface_topo_contribution_SH_coes;

        /**
         * A parameter to control whether to output the spherical harmonic coefficients of the CMB topography contribution
         */
        bool output_CMB_topo_contribution_SH_coes;

        /**
         * A parameter to control whether to output the spherical harmonic coefficients of the density anomaly
         */
        bool output_density_anomaly_contribution_SH_coes;

        /**
         * A parameter to control whether to output the free-air gravity anomaly
         */
        bool output_gravity_anomaly;

        /**
         * Parameters to set the density value out of the surface and CMB boundary
         */
        double density_above;
        double density_below;

        /**
         * A parameter to control whether to include the surface topography contribution on geoid
         */
        bool include_surface_topo_contribution;

        /**
         * A parameter to control whether to include the CMB topography contribution on geoid
         */
        bool include_CMB_topo_contribution;

        /**
         * A parameter to specify if the top boundary is an active free surface
         */
        bool use_free_surface_topography;

        /**
         * A parameter to specify if the bottom boundary is an active free surface
         */
        bool use_free_CMB_topography;

        /**
         * Function to compute the real spherical harmonic coefficients (cos and sin part) from min degree to max degree
         * The input spherical_function is a vector of vectors.
         * The inner vector stores theta, phi, spherical infinitesimal, and function value on a spherical surface.
         * The outer vector stores the inner vector associated with each quadrature point on a spherical surface.
         */
        std::pair<std::vector<double>,std::vector<double>>
        to_spherical_harmonic_coefficients(const std::vector<std::vector<double>> &spherical_function) const;

        /**
         * Function to compute the density contribution in spherical harmonic expansion throughout the mantle
         * The input outer radius is needed to evaluate the density integral contribution of whole model domain at the surface
         * This function returns a pair containing real spherical harmonics of density integral (cos and sin part) from min degree to max degree.
         */
        std::pair<std::vector<double>,std::vector<double>>
        density_contribution (const double &outer_radius) const;

        /**
         * Function to compute the surface and CMB topography contribution in spherical harmonic expansion
         * The input inner and outer radius are used to calculate spherical infinitesimal area, i.e., sin(theta)*d_theta*d_phi
         * associated with each quadrature point on surface and bottom respectively.
         * This function returns a pair containing surface and CMB topography's real spherical harmonic coefficients (cos and sin part)
         * from min degree to max degree. The surface and CMB average density are also included as the first single element of each subpair respectively.
         * The topography is based on the dynamic topography postprocessor in case of no free surface,
         * and based on the real surface from the geometry model in case of a free surface.
         */
        std::pair<std::pair<double, std::pair<std::vector<double>,std::vector<double>>>, std::pair<double, std::pair<std::vector<double>,std::vector<double>>>>
        topography_contribution(const double &outer_radius,
                                const double &inner_radius) const;

        /**
         * A vector to store the cosine terms of the geoid anomaly spherical harmonic coefficients.
         */
        std::vector<double> geoid_coecos;
        /**
         * A vector to store the sine terms of the geoid anomaly spherical harmonic coefficients.
         */
        std::vector<double> geoid_coesin;
    };
  }
}


#endif
