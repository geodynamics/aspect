/*
 Copyright (C) 2015 - 2017 by the authors of the ASPECT code.

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
 along with ASPECT; see the file doc/COPYING.  If not see
 <http://www.gnu.org/licenses/>.
 */


#ifndef __aspect__postprocess_geoid_h
#define __aspect__postprocess_geoid_h

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
        virtual
        std::pair<std::string,std::string>
        execute (TableHandler &statistics);

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

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
        bool also_output_geoid_anomaly_SH_coes;

        /**
         * A parameter to control whether to output the spherical harmonic coefficients of the surface dynamic topography contribution
         */
        bool also_output_surface_dynamic_topo_contribution_SH_coes;

        /**
         * A parameter to control whether to output the spherical harmonic coefficients of the CMB dynamic topography contribution
         */
        bool also_output_CMB_dynamic_topo_contribution_SH_coes;

        /**
         * A parameter to control whether to output the spherical harmonic coefficients of the density anomaly
         */
        bool also_output_density_anomaly_contribution_SH_coes;

        /**
         * Parameters to set the density value out of the surface and CMB boundary
         */
        double density_above;
        double density_below;

        /**
         * Function to compute the real spherical harmonic coefficients (cos and sin part) from min degree to max degree
         * The input spherical_function is a vector of vectors.
         * The inner vector stores theta, phi, spherical infinitesimal, and function value on a spherical surface.
         * The outer vector stores the inner vector associated with each quadrature point on a spherical surface.
         */
        std::pair<std::vector<double>,std::vector<double> >
        to_spherical_harmonic_coefficients(const std::vector<std::vector<double> > &spherical_function) const;

        /**
         * Function to compute the density contribution in spherical harmonic expansion throughout the mantle
         * The input outer radius is needed to evaluate the density integral contribution of whole model domain at the surface
         * This function returns a pair containing real spherical harmonics of density integral (cos and sin part) from min degree to max degree.
         */
        std::pair<std::vector<double>,std::vector<double> >
        density_contribution (const double &outer_radius) const;

        /**
         * Function to compute the surface and CMB dynamic topography contribution in spherical harmonic expansion
         * The input inner and outer radius are used to calculate spherical infinitesimal area, i.e., sin(theta)*d_theta*d_phi
         * associated with each quadrature point on surface and bottom respectively.
         * This function returns a pair containing surface and CMB dynamic topography's real spherical harmonic coefficients (cos and sin part)
         * from min degree to max degree. The surface and CMB average density are also included as the first single element of each subpair respectively.
         */
        std::pair<std::pair<double, std::pair<std::vector<double>,std::vector<double> > >, std::pair<double, std::pair<std::vector<double>,std::vector<double> > > >
        dynamic_topography_contribution(const double &outer_radius,
                                        const double &inner_radius) const;
    };
  }
}


#endif
