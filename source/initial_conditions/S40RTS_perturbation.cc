/*
  Copyright (C) 2011 - 2014 by the authors of the ASPECT code.

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


#include <aspect/initial_conditions/S40RTS_perturbation.h>
#include <aspect/initial_conditions/spline.h>
// The cubic spline interpolation script was downloaded from http://kluge.in-chemnitz.de/opensource/spline/spline.h

#include <fstream>
#include <iostream>

#include <boost/math/special_functions/spherical_harmonic.hpp>

namespace aspect
{
  namespace InitialConditions
  { 
    
    namespace internal
    {
       // Read in the spherical harmonics that are located in data/initial-conditions/S40RTS
       // and were downloaded from http://www.earth.lsa.umich.edu/~jritsema/research.html
       // Ritsema et al. choose real sine and cosine coefficients that follow the normalization
       // by Dahlen & Tromp, Theoretical Global Seismology (equations B.58 and B.99). 

       class SphericalHarmonicsLookup
       {
         public:
         SphericalHarmonicsLookup(const std::string &filename)
         {
           std::string temp;
           std::ifstream in(filename.c_str(), std::ios::in);
           AssertThrow (in,
                        ExcMessage (std::string("Couldn't open file <") + filename));

           in >> order;
           getline(in,temp);  // throw away the rest of the line    
          
           const int num_splines = 21;
           const int maxnumber = num_splines * (order+1)*(order+1);

           // read in all coefficients as a single data vector
           for (int i=0; i<maxnumber; i++)
           {
              double new_val;
              in >> new_val;
              coeffs.push_back(new_val);
           }

           // reorder the coefficients into sin and cos coefficients. a_lm will be the cos coefficients
           // and b_lm the sin coefficients.
           int ind = 0;
           int ind_degree;

           for (int j=0; j<num_splines; j++)

             for (int i=0; i<order+1; i++)
             {
               a_lm.push_back(coeffs[ind]);
               b_lm.push_back(0.0);
               ind += 1;
               
               ind_degree = 0;
               while (ind_degree < i)
               {
                 a_lm.push_back(coeffs[ind]);
                 ind += 1;
                 b_lm.push_back(coeffs[ind]);
                 ind += 1;
                 ind_degree +=1;
               }
             } 
           }
         std::vector<double> cos_coeffs()
         {
           return a_lm;
         }

         std::vector<double> sin_coeffs()
         {
           return b_lm;
         }

         int maxdegree()
         {
           return order;
         }

         private:
           int order;
           std::vector<double> coeffs;
           std::vector<double> a_lm;
           std::vector<double> b_lm;

       };

      // Read in the knot points for the spline interpolation. They are located in data/
      // initial-conditions/S40RTS and were taken from the plotting script 
      // lib/libS20/splhsetup.f which is part of the plotting package downloadable at
      // http://www.earth.lsa.umich.edu/~jritsema/research.html
      class SplineDepthsLookup
      {
         public:
         SplineDepthsLookup(const std::string &filename)
         {
           std::string temp;
           std::ifstream in(filename.c_str(), std::ios::in);
           AssertThrow (in,
       		ExcMessage (std::string("Couldn't open file <") + filename));

           getline(in,temp);  // throw away the rest of the line 
           getline(in,temp);  // throw away the rest of the line

           int num_splines = 21;

           for (int i=0; i<num_splines; i++)
           {
              double new_val;
              in >> new_val;

              depths.push_back(new_val);
           } 
         }

         std::vector<double> spline_depths()
         { 
           return depths;
         }

         private:
         std::vector<double> depths;
       };

     }

    template <int dim>
    void
    S40RTSPerturbation<dim>::initialize()
    {
      spherical_harmonics_lookup.reset(new internal::SphericalHarmonicsLookup(datadirectory+harmonics_coeffs_file_name));
      spline_depths_lookup.reset(new internal::SplineDepthsLookup(datadirectory+spline_depth_file_name));
    }

    // NOTE: this module uses the Boost spherical harmonics package which is not designed
    // for very high order (> 100) spherical harmonics computation. If you use harmonic
    // perturbations of a high order be sure to confirm the accuracy first.
    // For more information, see:
    // http://www.boost.org/doc/libs/1_49_0/libs/math/doc/sf_and_dist/html/math_toolkit/special/sf_poly/sph_harm.html
    
    template <int dim>
    double
    S40RTSPerturbation<dim>::
    initial_temperature (const Point<dim> &position) const
    {

      // use either the user-input reference temperature as background temperature
      // (incompressible model) or the adiabatic temperature profile (compressible model)
      const double background_temperature = this->get_material_model().is_compressible() ?
                                            this->get_adiabatic_conditions().temperature(position) :
                                            reference_temperature;
        
        //get the degree from the input file (20 or 40)
        const int maxdegree = spherical_harmonics_lookup->maxdegree();

        const int num_spline_knots = 21; // The tomography models are parameterized by 21 layers

        const int num_coeffs = (maxdegree+1) * (maxdegree+2) / 2 * num_spline_knots;

        // get the spherical harmonics coefficients
        const std::vector<double> a_lm = spherical_harmonics_lookup->cos_coeffs();
        const std::vector<double> b_lm = spherical_harmonics_lookup->sin_coeffs();
       
        // get spline knots and rescale them from [-1 1] to [CMB moho]
        const std::vector<double> r = spline_depths_lookup->spline_depths();
        const double rmoho = 6346e3;
        const double rcmb = 3480e3;
        std::vector<double> depth_values(num_spline_knots,0); 

        for (int i = 0; i<num_spline_knots; i++)
           depth_values[i] = rcmb+(rmoho-rcmb)*0.5*(r[i]+1);       

        // convert coordinates from [x,y,z] to [r, phi, theta] 
        const Tensor<1,dim> scoord = spherical_surface_coordinates(position);

        // iterate over all degrees and orders at each depth and sum them all up.
        std::vector<double> spline_values(num_spline_knots,0);
        double prefact;
        int ind = 0;

        for (int depth_interp = 0; depth_interp < num_spline_knots; depth_interp++)
        {
          for (int degree_l = 0; degree_l < maxdegree+1; degree_l++)
          {
            for (int order_m = 0; order_m < degree_l+1; order_m++)
            {
              const double cos_component = boost::math::spherical_harmonic_r(degree_l,order_m,scoord[2],scoord[1]); //real / cos part
              const double sin_component = boost::math::spherical_harmonic_i(degree_l,order_m,scoord[2],scoord[1]); //imaginary / sine part
                if (order_m == 0) {
                  // option to zero out degree 0, i.e. make sure that the average of the perturbation
                  // is 0 and the average of the temperature is the background temperature 
                  prefact = (zero_out_degree_0
                             ?
                             0.
                             :
                             1.);}
                else {
                  prefact = sqrt(2.);}
		spline_values[depth_interp] += prefact * (a_lm[ind]*cos_component + b_lm[ind]*sin_component);

             ind += 1;
           }
         }
       }

     // We need to reorder the spline_values because the coefficients are given from 
     // the surface down to the CMB and the interpolation knots range from the CMB up to
     // the surface.
     std::vector<double> spline_values_inv(num_spline_knots,0);
     for (int i=0; i<num_spline_knots; i++)
         spline_values_inv[i] = spline_values[num_spline_knots-1 - i];

     // The boundary condition for the cubic spline interpolation is that the function is linear 
     // at the boundary (i.e. moho and CMB). Values outside the range are linearly 
     // extrapolated.
     tk::spline s;
     s.set_points(depth_values,spline_values_inv);

     // Get value at specific depth
     const double perturbation = s(scoord[0]);

     // scale the perturbation in seismic velocity into a density perturbation
     // vs_to_density is an input parameter
     const double density_perturbation = vs_to_density * perturbation;

     // scale the density perturbation into a temperature perturbation
     const double temperature_perturbation =  -1./thermal_alpha * density_perturbation;
     const double temperature = background_temperature + temperature_perturbation;
     return temperature;

    }

    template <int dim>
    const Tensor<1,dim>
    S40RTSPerturbation<dim>::
    spherical_surface_coordinates(const Tensor<1,dim> &position) const
    {
      Tensor<1,dim> scoord;

      scoord[0] = std::sqrt(position.norm_square()); // R
      scoord[1] = std::atan2(position[1],position[0]); // Phi
      if (scoord[1] < 0.0) scoord[1] = 2*numbers::PI + scoord[1]; // correct phi to [0,2*pi]
      if (dim==3)
        scoord[2] = std::acos(position[2]/std::sqrt(position.norm_square())); // Theta

      return scoord;
    }  

    template <int dim>
    void
    S40RTSPerturbation<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial conditions");
      { 
          prm.enter_subsection("S40RTS perturbation");
          {
          prm.declare_entry("Data directory", "$ASPECT_SOURCE_DIR/data/initial-conditions/S40RTS/",
                            Patterns::DirectoryName (),
                             "The path to the model data. ");
          prm.declare_entry ("Initial condition file name", "S40RTS.sph",
                            Patterns::Anything(),
                             "The file name of the spherical harmonics coefficients"
                             "from Ritsema et al.");
          prm.declare_entry ("Spline knots depth file name", "Spline_knots.txt",
                            Patterns::Anything(),
                             "The file name of the spline knot locations from"
                             "Ritsema et al.");
          prm.declare_entry ("vs to density scaling", "0.25",
                             Patterns::Double (0),
                             " ");
          prm.declare_entry ("Thermal expansion coefficient in initial temperature scaling", "2e-5",
                             Patterns::Double (0),
                             "The value of the thermal expansion coefficient $\\beta$. "
                             "Units: $1/K$.");
          prm.declare_entry ("Remove degree 0 from perturbation","true",
                             Patterns::Bool (),
                             "Option to remove the degree zero component from the perturbation,"
                             "which will ensure that the depth-average temperature is equal to"
                             "the background temperature.");
          prm.declare_entry ("Reference temperature", "1600.0",
                             Patterns::Double (0),
                             "The reference temperature that is perturbed by the"
                             "harmonic function. Only used in incompressible models.");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    S40RTSPerturbation<dim>::parse_parameters (ParameterHandler &prm)
    {

      prm.enter_subsection("Initial conditions");
      {
        prm.enter_subsection("S40RTS perturbation");
        { 
          datadirectory           = prm.get ("Data directory");
          {
            const std::string      subst_text = "$ASPECT_SOURCE_DIR";
            std::string::size_type position;
            while (position = datadirectory.find (subst_text),  position!=std::string::npos)
              datadirectory.replace (datadirectory.begin()+position,
                                      datadirectory.begin()+position+subst_text.size(),
                                      ASPECT_SOURCE_DIR);
          }
          harmonics_coeffs_file_name = prm.get ("Initial condition file name");
          spline_depth_file_name  = prm.get ("Spline knots depth file name");
          vs_to_density           = prm.get_double ("vs to density scaling");
          thermal_alpha           = prm.get_double ("Thermal expansion coefficient in initial temperature scaling");
          zero_out_degree_0       = prm.get_bool ("Remove degree 0 from perturbation");
          reference_temperature   = prm.get_double ("Reference temperature");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    
      initialize ();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace InitialConditions
  {
    ASPECT_REGISTER_INITIAL_CONDITIONS(S40RTSPerturbation,
                                       "S40RTS perturbation",
                                       "An initial temperature field in which the temperature "
                                       "is perturbed following the S20RTS or S40RTS shear wave "
                                       "velocity model by Ritsema and others, which can be downloaded" 
                                       "here http://www.earth.lsa.umich.edu/~jritsema/research.html"
                                       "Information on the vs model can be found in Ritsema, J., Deuss," 
                                       "A., van Heijst, H.J. & Woodhouse, J.H., 2011. S40RTS: a" 
                                       "degree-40 shear-velocity model for the mantle from new Rayleigh" 
                                       "wave dispersion, teleseismic traveltime and normal-mode" 
                                       "splitting function measurements, Geophys. J. Int. 184, 1223-1236." 
                                       "The scaling between the shear wave perturbation and the"
                                       "temperature perturbation can be set by the user with the" 
                                       "'vs to density scaling' parameter and the 'Thermal" 
                                       "expansion coefficient in initial temperature scaling'" 
                                       "parameter. The scaling is as follows: $\\delta ln \\rho"
                                       "(r,\\theta,\\phi) = \\xi \\cdot \\delta ln v_s(r,\\theta,"
                                       "\\phi)$ and $\\delta T(r,\\theta,\\phi) = - \\frac{1}{\\alpha}"
                                       "\\delta ln \\rho(r,\\theta,\\phi)$. $\\xi$ is the 'vs to"
                                       "density scaling' parameter and $\\alpha$ is the 'Thermal" 
                                       "expansion coefficient in initial temperature scaling'" 
                                       "parameter. The temperature perturbation is added to an" 
                                       "otherwise constant temperature (incompressible model) or" 
                                       "adiabatic reference profile (compressible model).")
  }
}
