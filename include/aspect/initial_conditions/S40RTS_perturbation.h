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


#ifndef __aspect__initial_conditions_S40RTS_perturbation_h
#define __aspect__initial_conditions_S40RTS_perturbation_h

#include <aspect/initial_conditions/interface.h>
#include <aspect/simulator.h>


namespace aspect
{
  namespace InitialConditions
  {
    using namespace dealii;

    namespace internal
    {

       class SphericalHarmonicsLookup;
       class SplineDepthsLookup;
    }

    /**
     * A class that describes a perturbed initial temperature field for a spherical
     * shell geometry model. The perturbation is based on the S20RTS / S40RTS 
     * global shear wave velocity model by Ritsema et al.
     * http://www.earth.lsa.umich.edu/~jritsema/research.html
     *
     * @ingroup InitialConditionsModels
     */

    template <int dim>
    class S40RTSPerturbation : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Initialization function. Loads the material data and sets up
         * pointers.
         */
        void
        initialize ();

         /**
         * Return the initial temperature as a function of position.
         */
        virtual
        double initial_temperature (const Point<dim> &position) const;

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
         * Returns spherical coordinates of a cartesian position.
         */
        const Tensor<1,dim>
        spherical_surface_coordinates(const Tensor<1,dim> &position) const;
        
        /**
         * File directory and names
         */
        std::string datadirectory;
        std::string spline_depth_file_name;
    
        /**
         * This parameter allows setting the input file for the shear-wave perturbation. Options so far
         * are S20RTS.sph and S40RTS.sph. For S40RTS there are different versions available that differ
         * by the degree of damping in the seismic inversion. These models could be downloaded and used
         * as well. 
         */
        std::string harmonics_coeffs_file_name;

        /**
         * The parameters below describe the perturbation of shear wave velocity into a temperatures perturbation
         * The first parameter is constant so far but could be made depth dependent as constraint
         * by e.g. Forte, A.M. & Woodward, R.L., 1997. Seismic-geodynamic constraints on three-
         * dimensional structure, vertical flow, and heat transfer in the mantle, J. Geophys. Res.
         * 102 (B8), 17,981-17,994. 
         */
        double vs_to_density;
        double thermal_alpha;

        /**
         * This parameter allows to set the degree 0 component of the shear wave velocity perturbation to 
         * zero, which guarantees that average temperature at a certain depth is the background temperature.
         */
        bool zero_out_degree_0;

        /**
         * This parameter gives the reference temperature, which will be perturbed. In the compressional case
         * the background temperature will be the adiabat.
         */
        double reference_temperature;

        std_cxx1x::shared_ptr<internal::SphericalHarmonicsLookup> spherical_harmonics_lookup;
        std_cxx1x::shared_ptr<internal::SplineDepthsLookup> spline_depths_lookup;


    };



    // tk does the cubic spline interpolation.
    // This interpolation is based on the script spline.h, which was downloaded from 
    // http://kluge.in-chemnitz.de/opensource/spline/spline.h   //
    // Copyright (C) 2011, 2014 Tino Kluge (ttk448 at gmail.com)
        
    namespace tk {
        
        // band matrix solver
        class band_matrix {
            private:
            std::vector< std::vector<double> > m_upper;  // upper band
            std::vector< std::vector<double> > m_lower;  // lower band
            public:
            band_matrix() {};                             // constructor
            band_matrix(int dim, int n_u, int n_l);       // constructor
            ~band_matrix() {};                            // destructor
            void resize(int dim, int n_u, int n_l);      // init with dim,n_u,n_l
            int dim() const;                             // matrix dimension
            int num_upper() const {
                return m_upper.size()-1;
            }
            int num_lower() const {
                return m_lower.size()-1;
            }
            // access operator
            double & operator () (int i, int j);            // write
            double   operator () (int i, int j) const;      // read
            // we can store an additional diogonal (in m_lower)
            double& saved_diag(int i);
            double  saved_diag(int i) const;
            void lu_decompose();
            std::vector<double> r_solve(const std::vector<double>& b) const;
            std::vector<double> l_solve(const std::vector<double>& b) const;
            std::vector<double> lu_solve(const std::vector<double>& b,
                                         bool is_lu_decomposed=false);
            
        };
        
        
        // spline interpolation
         class spline {
            private:
            std::vector<double> m_x,m_y;           // x,y coordinates of points
            // interpolation parameters
            // f(x) = a*(x-x_i)^3 + b*(x-x_i)^2 + c*(x-x_i) + y_i
            std::vector<double> m_a,m_b,m_c,m_d;
            public:
            void set_points(const std::vector<double>& x,
                            const std::vector<double>& y, bool cubic_spline=true);
            double operator() (double x) const;
        };
     }
  }
}

#endif
