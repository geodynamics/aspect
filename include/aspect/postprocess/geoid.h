/*
 Copyright (C) 2011, 2012 by the authors of the ASPECT code.
 
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


#ifndef __aspect__postprocess_geoid_anomaly_h
#define __aspect__postprocess_geoid_anomaly_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace Postprocess
  {
      
    /**
     * A postprocessor that computes geoid anomaly at the surface.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class Geoid : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Evaluate the solution for the geoid anomaly.
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
        int max_degree;
        int min_degree;
        
        /**
         * A parameter to control the resolution of the density anomaly depth integral
         */
        int depth_integral_refinement;
        
        /**
         * A parameter to control whether to output the data in geographical coordinates.
         * If true, output the data in longitudes and latitudes; if false, output data in x y z.
         */
        bool output_in_lat_lon;
        
        /**
         * Parameters to set the density value out of the surface and CMB boundary
         */
        double density_above;
        double density_below;
        
        /**
         * Function to compute the spherical harmonic coefficients
         */
        std::pair<std::vector<double>,std::vector<double> >
        sph_fun2coes(const std::vector<std::vector<double> > &spherical_function) const;
        
        /**
         * Function to compute the surface *normal stress over gravity* in spherical harmonic expansion
         */
        std::pair<std::vector<double>,std::vector<double> >
        get_surface_dynamic_topography() const;
        
        /**
         * Function to compute the CMB *normal stress over gravity* in spherical harmonic expansion
         */
        std::pair<std::vector<double>,std::vector<double> >
        get_CMB_dynamic_topography() const;
        
        /**
         * Function to compute the density anomaly on a global sphere within the 
         * mantle in spherical harmonic expansion
         */
        std::pair<double,std::pair<std::vector<double>,std::vector<double> > >
        get_sphere_density_anomaly (const double &r_lower_sphere,
                                    const double &r_upper_sphere) const;
    };
  }
}


#endif
