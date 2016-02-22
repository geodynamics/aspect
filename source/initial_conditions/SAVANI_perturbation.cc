/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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


#include <aspect/initial_conditions/SAVANI_perturbation.h>
#include <aspect/utilities.h>
#include <fstream>
#include <iostream>
#include <deal.II/base/std_cxx11/array.h>

#include <boost/math/special_functions/spherical_harmonic.hpp>

namespace aspect
{
  namespace InitialConditions
  {
    namespace internal
    {
      namespace SAVANI
      {
        // Read in the spherical harmonics that are located in data/initial-conditions/SAVANI
        // and were downloaded from http://n.ethz.ch/~auerl/research.html
        // choose real sine and cosine coefficients that follow the normalization
        // by Dahlen & Tromp, Theoretical Global Seismology (equations B.58 and B.99).

        class SphericalHarmonicsLookup
        {
          public:
            SphericalHarmonicsLookup(const std::string &filename,
                                     const MPI_Comm &comm)
            {
              std::string temp;
              // Read data from disk and distribute among processes
              std::istringstream in(Utilities::read_and_distribute_file_content(filename, comm));

              in >> order;
              getline(in,temp);  // throw away the rest of the line

              const int num_layers = 28;
              // const int maxnumber = num_layers * (order+1)*(order+2);

              // read in all coefficients as a single data vector
              for (int i=0; i<num_layers; i++)
                {
                  for (int j=0; j<(order+1)*(order+2); j++)
                    {
                      double new_val;
                      in >> new_val;
                      coeffs.push_back(0.01*new_val);
                    }
                  getline(in,temp);
                }
              // reorder the coefficients into sin and cos coefficients. a_lm will be the cos coefficients
              // and b_lm the sin coefficients.
              int ind = 0;
              int ind_degree;

              for (int j=0; j<num_layers; j++)

                for (int i=0; i<order+1; i++)
                  {
                    //a_lm.push_back(coeffs[ind]);
                    //b_lm.push_back(0.0);
                    //ind += 1;

                    ind_degree = 0;
                    while (ind_degree <= i)
                      {
                        a_lm.push_back(coeffs[ind]);
                        ind += 1;
                        b_lm.push_back(coeffs[ind]);
                        ind += 1;
                        ind_degree +=1;
                      }
                  }
            }

            // Declare a function that returns the cosine coefficients
            const std::vector<double> &cos_coeffs() const
            {
              return a_lm;
            }

            // Declare a function that returns the sine coefficients
            const std::vector<double> &sin_coeffs() const
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
        // initial-conditions/SAVANI and were taken from the 28 spherical layers of SAVANI
        // tomography model by a matlab script convert_to_knots.m located in the same directory.
        class SplineDepthsLookup
        {
          public:
            SplineDepthsLookup(const std::string &filename,
                               const MPI_Comm &comm)
            {
              std::string temp;
              // Read data from disk and distribute among processes
              std::istringstream in(Utilities::read_and_distribute_file_content(filename, comm));

              getline(in,temp);  // throw away the rest of the line
              getline(in,temp);  // throw away the rest of the line

              int num_splines = 28;

              for (int i=0; i<num_splines; i++)
                {
                  double new_val;
                  in >> new_val;
                  depths.push_back(new_val);
                }
            }

            const std::vector<double> &spline_depths() const
            {
              return depths;
            }

          private:
            std::vector<double> depths;
        };
      }
    }


    template <int dim>
    void
    SAVANIPerturbation<dim>::initialize()
    {
      spherical_harmonics_lookup.reset(new internal::SAVANI::SphericalHarmonicsLookup(datadirectory+harmonics_coeffs_file_name,this->get_mpi_communicator()));
      spline_depths_lookup.reset(new internal::SAVANI::SplineDepthsLookup(datadirectory+spline_depth_file_name,this->get_mpi_communicator()));
    }

    // NOTE: this module uses the Boost spherical harmonics package which is not designed
    // for very high order (> 100) spherical harmonics computation. If you use harmonic
    // perturbations of a high order be sure to confirm the accuracy first.
    // For more information, see:
    // http://www.boost.org/doc/libs/1_49_0/libs/math/doc/sf_and_dist/html/math_toolkit/special/sf_poly/sph_harm.html

    template <>
    double
    SAVANIPerturbation<2>::
    initial_temperature (const Point<2> &) const
    {
      // we shouldn't get here but instead should already have been
      // kicked out by the assertion in the parse_parameters()
      // function
      Assert (false, ExcNotImplemented());
      return 0;
    }


    template <>
    double
    SAVANIPerturbation<3>::
    initial_temperature (const Point<3> &position) const
    {
      const unsigned int dim = 3;

      // use either the user-input reference temperature as background temperature
      // (incompressible model) or the adiabatic temperature profile (compressible model)
      const double background_temperature = this->get_material_model().is_compressible() ?
                                            this->get_adiabatic_conditions().temperature(position) :
                                            reference_temperature;

      //get the degree from the input file (60)
      const int maxdegree = spherical_harmonics_lookup->maxdegree();
      // const int maxdegree = 60;

      const int num_spline_knots = 28; // The tomography models are parameterized by 28 layers

      // get the spherical harmonics coefficients
      const std::vector<double> a_lm = spherical_harmonics_lookup->cos_coeffs();
      const std::vector<double> b_lm = spherical_harmonics_lookup->sin_coeffs();

      // get spline knots and rescale them from [-1 1], i.e., CMB to moho.
      const std::vector<double> r = spline_depths_lookup->spline_depths();
      const double rmoho = 6346e3;
      const double rcmb = 3480e3;
      std::vector<double> depth_values(num_spline_knots,0);

      for (int i = 0; i<num_spline_knots; i++)
        depth_values[i] = rcmb+(rmoho-rcmb)*0.5*(r[i]+1);

      // convert coordinates from [x,y,z] to [r, phi, theta]
      std_cxx11::array<double,dim> scoord = aspect::Utilities::spherical_coordinates(position);

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

                  // normalization after Dahlen and Tromp, 1986, Appendix B.6
                  if (degree_l == 0)
                    prefact = (zero_out_degree_0
                               ?
                               0.
                               :
                               1.);
                  else if (order_m == 0)
                    prefact = 1.;
                  else
                    prefact = sqrt(2.);

                  spline_values[depth_interp] += prefact * (a_lm[ind]*cos_component + b_lm[ind]*sin_component);

                  ind += 1;
                }
            }
        }


      // The boundary condition for the cubic spline interpolation is that the function is linear
      // at the boundary (i.e. moho and CMB). Values outside the range are linearly
      // extrapolated.
      aspect::Utilities::tk::spline s;
      s.set_points(depth_values,spline_values);

      // Get value at specific depth
      const double perturbation = s(scoord[0]);

      // scale the perturbation in seismic velocity into a density perturbation
      // vs_to_density is an input parameter
      const double density_perturbation = vs_to_density * perturbation;

      const double depth = this->get_geometry_model().depth(position);
      double temperature_perturbation;
      if (depth > no_perturbation_depth)
        // scale the density perturbation into a temperature perturbation
        temperature_perturbation =  -1./thermal_alpha * density_perturbation;
      else
        // set heterogeneity to zero down to a specified depth
        temperature_perturbation = 0.0;

      // add the temperature perturbation to the background temperature
      return background_temperature + temperature_perturbation;
    }


    template <int dim>
    void
    SAVANIPerturbation<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial conditions");
      {
        prm.enter_subsection("SAVANI perturbation");
        {
          prm.declare_entry("Data directory", "$ASPECT_SOURCE_DIR/data/initial-conditions/SAVANI/",
                            Patterns::DirectoryName (),
                            "The path to the model data. ");
          prm.declare_entry ("Initial condition file name", "savani.dlnvs.60.m.ab",
                             Patterns::Anything(),
                             "The file name of the spherical harmonics coefficients "
                             "from Auer et al.");
          prm.declare_entry ("Spline knots depth file name", "Spline_knots.txt",
                             Patterns::Anything(),
                             "The file name of the spline knots taken from the 28 spherical layers"
                             "of SAVANI tomography model.");
          prm.declare_entry ("vs to density scaling", "0.25",
                             Patterns::Double (0),
                             "This parameter specifies how the perturbation in shear wave velocity "
                             "as prescribed by SAVANI is scaled into a density perturbation. "
                             "See the general description of this model for more detailed information.");
          prm.declare_entry ("Thermal expansion coefficient in initial temperature scaling", "2e-5",
                             Patterns::Double (0),
                             "The value of the thermal expansion coefficient $\\beta$. "
                             "Units: $1/K$.");
          prm.declare_entry ("Remove degree 0 from perturbation","true",
                             Patterns::Bool (),
                             "Option to remove the degree zero component from the perturbation, "
                             "which will ensure that the laterally averaged temperature for a fixed "
                             "depth is equal to the background temperature.");
          prm.declare_entry ("Reference temperature", "1600.0",
                             Patterns::Double (0),
                             "The reference temperature that is perturbed by the spherical "
                             "harmonic functions. Only used in incompressible models.");
          prm.declare_entry ("Remove temperature heterogeneity down to specified depth", boost::lexical_cast<std::string>(-std::numeric_limits<double>::max()),
                             Patterns::Double (),
                             "This will set the heterogeneity prescribed by SAVANI to zero "
                             "down to the specified depth (in meters). Note that your resolution has "
                             "to be adquate to capture this cutoff. For example if you specify a depth "
                             "of 660km, but your closest spherical depth layers are only at 500km and "
                             "750km (due to a coarse resolution) it will only zero out heterogeneities "
                             "down to 500km. Similar caution has to be taken when using adaptive meshing.");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    SAVANIPerturbation<dim>::parse_parameters (ParameterHandler &prm)
    {
      AssertThrow (dim == 3,
                   ExcMessage ("The 'S40RTS perturbation' model for the initial "
                               "temperature is only available for 3d computations."));

      prm.enter_subsection("Initial conditions");
      {
        prm.enter_subsection("SAVANI perturbation");
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
          no_perturbation_depth   = prm.get_double ("Remove temperature heterogeneity down to specified depth");
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
    ASPECT_REGISTER_INITIAL_CONDITIONS(SAVANIPerturbation,
                                       "SAVANI perturbation",
                                       "An initial temperature field in which the temperature "
                                       "is perturbed following the SAVANI shear wave "
                                       "velocity model by Auer and others, which can be downloaded "
                                       "here \\url{http://n.ethz.ch/~auerl/savani.tar.bz2}. "
                                       "Information on the vs model can be found in Auer, L., Boschi, "
                                       "L., Becker, T.W., Nissen-Meyer, T. \\& Giardini, D., 2014. Savani:"
                                       "A variable resolution whole‐mantle model of anisotropic shear velocity"
                                       "variations based on multiple data sets. Journal of Geophysical"
                                       "Research: Solid Earth 119.4 (2014): 3006-3034. "
                                       "The scaling between the shear wave perturbation and the "
                                       "temperature perturbation can be set by the user with the "
                                       "'vs to density scaling' parameter and the 'Thermal "
                                       "expansion coefficient in initial temperature scaling' "
                                       "parameter. The scaling is as follows: $\\delta ln \\rho "
                                       "(r,\\theta,\\phi) = \\xi \\cdot \\delta ln v_s(r,\\theta, "
                                       "\\phi)$ and $\\delta T(r,\\theta,\\phi) = - \\frac{1}{\\alpha} "
                                       "\\delta ln \\rho(r,\\theta,\\phi)$. $\\xi$ is the 'vs to "
                                       "density scaling' parameter and $\\alpha$ is the 'Thermal "
                                       "expansion coefficient in initial temperature scaling' "
                                       "parameter. The temperature perturbation is added to an "
                                       "otherwise constant temperature (incompressible model) or "
                                       "adiabatic reference profile (compressible model).")
  }
}
