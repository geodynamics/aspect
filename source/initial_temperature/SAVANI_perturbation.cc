/*
  Copyright (C) 2011 - 2022 by the authors of the ASPECT code.

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


#include <aspect/initial_temperature/SAVANI_perturbation.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/utilities.h>
#include <aspect/simulator_access.h>
#include <aspect/initial_composition/interface.h>
#include <aspect/material_model/interface.h>
#include <fstream>
#include <iostream>
#include <array>

#include <boost/lexical_cast.hpp>

namespace aspect
{
  namespace InitialTemperature
  {
    namespace internal
    {
      namespace SAVANI
      {
        // Read in the spherical harmonics that are located in data/initial-temperature/SAVANI
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
              std::getline(in,temp);  // throw away the rest of the line

              const unsigned int num_layers = 28;

              // read in all coefficients as a single data vector
              std::vector<double> coeffs;
              for (unsigned int i=0; i<num_layers; ++i)
                {
                  for (unsigned int j=0; j<(order+1)*(order+2); ++j)
                    {
                      double new_val;
                      in >> new_val;
                      coeffs.push_back(0.01*new_val);
                    }
                  std::getline(in,temp);
                }

              // reorder the coefficients into sin and cos coefficients. a_lm will be the cos coefficients
              // and b_lm the sin coefficients.
              unsigned int ind = 0;

              for (unsigned int j=0; j<num_layers; ++j)
                for (unsigned int i=0; i<order+1; ++i)
                  {
                    unsigned int ind_degree = 0;
                    while (ind_degree <= i)
                      {
                        a_lm.push_back(coeffs[ind]);
                        ++ind;
                        b_lm.push_back(coeffs[ind]);
                        ++ind;
                        ++ind_degree;
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

            unsigned int maxdegree()
            {
              return order;
            }

          private:
            unsigned int order;
            std::vector<double> a_lm;
            std::vector<double> b_lm;

        };

        // Read in the knot points for the spline interpolation. They are located in data/
        // initial-temperature/SAVANI and were taken from the 28 spherical layers of SAVANI
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

              std::getline(in,temp);  // throw away the rest of the line
              std::getline(in,temp);  // throw away the rest of the line

              const unsigned int num_splines = 28;
              depths.resize(num_splines);
              for (unsigned int i=0; i<num_splines; ++i)
                {
                  in >> depths[i];
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
    SAVANIPerturbation<dim>::SAVANIPerturbation()
      :
      vs_to_density_index(numbers::invalid_unsigned_int)
    {}


    template <int dim>
    void
    SAVANIPerturbation<dim>::initialize()
    {
      spherical_harmonics_lookup
        = std::make_unique<internal::SAVANI::SphericalHarmonicsLookup>(data_directory+harmonics_coeffs_file_name,
                                                                       this->get_mpi_communicator());
      spline_depths_lookup
        = std::make_unique<internal::SAVANI::SplineDepthsLookup>(data_directory+spline_depth_file_name,
                                                                 this->get_mpi_communicator());

      if (vs_to_density_method == file)
        {
          profile.initialize(this->get_mpi_communicator());
          vs_to_density_index = profile.get_column_index_from_name("vs_to_density");
        }

    }


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
    SAVANIPerturbation<2>::
    get_Vs (const Point<2> &/*position*/) const
    {
      Assert (false, ExcNotImplemented());
      return 0;
    }

    template <>
    double
    SAVANIPerturbation<3>::
    get_Vs (const Point<3> &position) const
    {
      const unsigned int dim = 3;

      // get the degree from the input file (60)
      unsigned int max_degree = spherical_harmonics_lookup->maxdegree();

      // lower the maximum order if needed
      if (lower_max_order)
        {
          AssertThrow(max_order <= max_degree, ExcMessage("Specifying a maximum order higher than the order of spherical harmonic data is not allowed"));
          max_degree = max_order;
        }

      const int num_spline_knots = 28; // The tomography models are parameterized by 28 layers

      // get the spherical harmonics coefficients
      const std::vector<double> &a_lm = spherical_harmonics_lookup->cos_coeffs();
      const std::vector<double> &b_lm = spherical_harmonics_lookup->sin_coeffs();

      // get spline knots and rescale them from [-1 1], i.e., CMB to Moho.
      const std::vector<double> &r = spline_depths_lookup->spline_depths();
      const double rmoho = 6346e3;
      const double rcmb = 3480e3;
      std::vector<double> depth_values(num_spline_knots, 0.);

      for (unsigned int i = 0; i<num_spline_knots; ++i)
        depth_values[i] = rcmb+(rmoho-rcmb)*0.5*(r[i]+1.);

      // convert coordinates from [x,y,z] to [r, phi, theta]
      std::array<double,dim> scoord = aspect::Utilities::Coordinates::cartesian_to_spherical_coordinates(position);

      // Evaluate the spherical harmonics at this position. Since they are the
      // same for all depth splines, do it once to avoid multiple evaluations.
      std::vector<std::vector<double>> cosine_components(max_degree+1, std::vector<double>(max_degree+1, 0.0));
      std::vector<std::vector<double>> sine_components(max_degree+1, std::vector<double>(max_degree+1, 0.0));

      for (unsigned int degree_l = 0; degree_l < max_degree+1; ++degree_l)
        {
          for (unsigned int order_m = 0; order_m < degree_l+1; ++order_m)
            {
              const std::pair<double,double> sph_harm_vals = Utilities::real_spherical_harmonic( degree_l, order_m, scoord[2], scoord[1] );
              cosine_components[degree_l][order_m] = sph_harm_vals.first;
              sine_components[degree_l][order_m] = sph_harm_vals.second;
            }
        }

      // iterate over all degrees and orders at each depth and sum them all up.
      std::vector<double> spline_values(num_spline_knots, 0.);
      double prefact;
      unsigned int ind = 0;

      for (unsigned int depth_interp = 0; depth_interp < num_spline_knots; ++depth_interp)
        {
          for (unsigned int degree_l = 0; degree_l < max_degree+1; ++degree_l)
            {
              for (unsigned int order_m = 0; order_m < degree_l+1; ++order_m)
                {
                  // normalization after Dahlen and Tromp, 1986, Appendix B.6
                  if (degree_l == 0)
                    prefact = (zero_out_degree_0
                               ?
                               0.
                               :
                               1.);
                  else
                    prefact = 1.0;

                  spline_values[depth_interp] += prefact * (a_lm[ind] * cosine_components[degree_l][order_m]
                                                            + b_lm[ind] * sine_components[degree_l][order_m]);

                  ++ind;
                }
            }
        }


      // The boundary condition for the cubic spline interpolation is that the function is linear
      // at the boundary (i.e. Moho and CMB). Values outside the range are linearly
      // extrapolated.
      aspect::Utilities::tk::spline s;
      s.set_points(depth_values, spline_values);

      // Get value at specific depth
      return s(scoord[0]);

    }


    template <>
    double
    SAVANIPerturbation<3>::
    initial_temperature (const Point<3> &position) const
    {
      // use either the user-input reference temperature as background temperature
      // (incompressible model) or the adiabatic temperature profile (compressible model)
      const double background_temperature = this->get_material_model().is_compressible() ?
                                            this->get_adiabatic_conditions().temperature(position) :
                                            reference_temperature;

      //Read in Vs perturbation data using function above
      const double perturbation = get_Vs (position);

      // Get the vs to density conversion
      const double depth = this->get_geometry_model().depth(position);

      double vs_to_density = 0.0;
      if (vs_to_density_method == file)
        vs_to_density = profile.get_data_component(Point<1>(depth), vs_to_density_index);
      else if (vs_to_density_method == constant)
        vs_to_density = vs_to_density_constant;
      else
        // we shouldn't get here but instead should already have been
        // kicked out by the assertion in the parse_parameters()
        // function
        Assert (false, ExcNotImplemented());

      // scale the perturbation in seismic velocity into a density perturbation
      // vs_to_density is an input parameter
      const double density_perturbation = vs_to_density * perturbation;

      double temperature_perturbation;
      if (depth > no_perturbation_depth)
        {
          // scale the density perturbation into a temperature perturbation
          // see if we need to ask material model for the thermal expansion coefficient
          if (use_material_model_thermal_alpha)
            {
              MaterialModel::MaterialModelInputs<3> in(1, this->n_compositional_fields());
              MaterialModel::MaterialModelOutputs<3> out(1, this->n_compositional_fields());
              in.position[0] = position;
              in.temperature[0] = background_temperature;
              in.pressure[0] = this->get_adiabatic_conditions().pressure(position);
              in.velocity[0] = Tensor<1,3> ();
              for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                in.composition[0][c] = this->get_initial_composition_manager().initial_composition(position, c);
              in.requested_properties = MaterialModel::MaterialProperties::thermal_expansion_coefficient;

              this->get_material_model().evaluate(in, out);

              temperature_perturbation = -1./(out.thermal_expansion_coefficients[0]) * density_perturbation;
            }
          else
            temperature_perturbation = -1./thermal_alpha * density_perturbation;
        }
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
      prm.enter_subsection ("Initial temperature model");
      {
        prm.enter_subsection("SAVANI perturbation");
        {
          prm.declare_entry ("Data directory", "$ASPECT_SOURCE_DIR/data/initial-temperature/SAVANI/",
                             Patterns::DirectoryName (),
                             "The path to the model data.");
          prm.declare_entry ("Initial condition file name", "savani.dlnvs.60.m.ab",
                             Patterns::Anything(),
                             "The file name of the spherical harmonics coefficients "
                             "from Auer et al.");
          prm.declare_entry ("Spline knots depth file name", "Spline_knots.txt",
                             Patterns::Anything(),
                             "The file name of the spline knots taken from the 28 spherical layers "
                             "of SAVANI tomography model.");
          prm.declare_entry ("Vs to density scaling method", "constant",
                             Patterns::Selection("file|constant"),
                             "Method that is used to specify how the vs-to-density scaling varies "
                             "with depth.");
          prm.declare_entry ("Vs to density scaling", "0.25",
                             Patterns::Double (0.),
                             "This parameter specifies how the perturbation in shear wave velocity "
                             "as prescribed by SAVANI is scaled into a density perturbation. "
                             "See the general description of this model for more detailed information.");
          prm.declare_entry ("Thermal expansion coefficient in initial temperature scaling", "2e-5",
                             Patterns::Double (0.),
                             "The value of the thermal expansion coefficient $\\beta$. "
                             "Units: \\si{\\per\\kelvin}.");
          prm.declare_entry ("Use thermal expansion coefficient from material model", "false",
                             Patterns::Bool (),
                             "Option to take the thermal expansion coefficient from the "
                             "material model instead of from what is specified in this "
                             "section.");
          prm.declare_entry ("Remove degree 0 from perturbation","true",
                             Patterns::Bool (),
                             "Option to remove the degree zero component from the perturbation, "
                             "which will ensure that the laterally averaged temperature for a fixed "
                             "depth is equal to the background temperature.");
          prm.declare_entry ("Reference temperature", "1600.0",
                             Patterns::Double (0.),
                             "The reference temperature that is perturbed by the spherical "
                             "harmonic functions. Only used in incompressible models.");
          prm.declare_entry ("Remove temperature heterogeneity down to specified depth",
                             boost::lexical_cast<std::string>(std::numeric_limits<double>::lowest()),
                             Patterns::Double (),
                             "This will set the heterogeneity prescribed by SAVANI to zero "
                             "down to the specified depth (in meters). Note that your resolution has "
                             "to be adequate to capture this cutoff. For example if you specify a depth "
                             "of 660km, but your closest spherical depth layers are only at 500km and "
                             "750km (due to a coarse resolution) it will only zero out heterogeneities "
                             "down to 500km. Similar caution has to be taken when using adaptive meshing.");
          prm.declare_entry ("Specify a lower maximum order","false",
                             Patterns::Bool (),
                             "Option to use a lower maximum order when reading the data file of spherical "
                             "harmonic coefficients. This is probably used for the faster tests or when the "
                             "users only want to see the spherical harmonic pattern up to a certain order.");
          prm.declare_entry ("Maximum order","20",
                             Patterns::Integer (0),
                             "The maximum order the users specify when reading the data file of spherical harmonic "
                             "coefficients, which must be smaller than the maximum order the data file stored. "
                             "This parameter will be used only if 'Specify a lower maximum order' is set to true.");
          aspect::Utilities::AsciiDataProfile<dim>::declare_parameters(prm,
                                                                       "$ASPECT_SOURCE_DIR/data/initial-temperature/S40RTS/",
                                                                       "vs_to_density_Steinberger.txt",
                                                                       "Ascii data vs to density model");
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
                   ExcMessage ("The 'SAVANI perturbation' model for the initial "
                               "temperature is only available for 3d computations."));

      prm.enter_subsection ("Initial temperature model");
      {
        prm.enter_subsection("SAVANI perturbation");
        {
          data_directory = Utilities::expand_ASPECT_SOURCE_DIR(prm.get ("Data directory"));
          if ((data_directory.size() > 0) && (data_directory[data_directory.size()-1] != '/'))
            data_directory += "/";
          harmonics_coeffs_file_name = prm.get ("Initial condition file name");
          spline_depth_file_name  = prm.get ("Spline knots depth file name");
          vs_to_density_constant           = prm.get_double ("Vs to density scaling");
          thermal_alpha           = prm.get_double ("Thermal expansion coefficient in initial temperature scaling");
          use_material_model_thermal_alpha = prm.get_bool ("Use thermal expansion coefficient from material model");
          zero_out_degree_0       = prm.get_bool ("Remove degree 0 from perturbation");
          reference_temperature   = prm.get_double ("Reference temperature");
          no_perturbation_depth   = prm.get_double ("Remove temperature heterogeneity down to specified depth");
          lower_max_order         = prm.get_bool ("Specify a lower maximum order");
          max_order               = prm.get_integer ("Maximum order");

          if (prm.get("Vs to density scaling method") == "file")
            vs_to_density_method = file;
          else if (prm.get("Vs to density scaling method") == "constant")
            vs_to_density_method = constant;
          else
            {
              AssertThrow(false, ExcMessage("Unknown method for vs to density scaling."));
            }

          profile.parse_parameters(prm,"Ascii data vs to density model");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace InitialTemperature
  {
    ASPECT_REGISTER_INITIAL_TEMPERATURE_MODEL(SAVANIPerturbation,
                                              "SAVANI perturbation",
                                              "An initial temperature field in which the temperature "
                                              "is perturbed following the SAVANI shear wave "
                                              "velocity model by Auer and others, which can be downloaded "
                                              "here \\url{http://n.ethz.ch/~auerl/savani.tar.bz2}. "
                                              "Information on the vs model can be found in Auer, L., Boschi, "
                                              "L., Becker, T.W., Nissen-Meyer, T. \\& Giardini, D., 2014. "
                                              "Savani: A variable resolution whole-mantle model of anisotropic "
                                              "shear velocity variations based on multiple data sets. Journal "
                                              "of Geophysical Research: Solid Earth 119.4 (2014): 3006-3034. "
                                              "The scaling between the shear wave perturbation and the "
                                              "density perturbation can be constant and set by the user with the "
                                              "'Vs to density scaling' parameter or depth-dependent and "
                                              "read in from a file. To convert density the user can specify "
                                              "the 'Thermal expansion coefficient in initial temperature scaling' "
                                              "parameter. The scaling is as follows: $\\delta \\ln \\rho "
                                              "(r,\\theta,\\phi) = \\xi \\cdot \\delta \\ln v_s(r,\\theta, "
                                              "\\phi)$ and $\\delta T(r,\\theta,\\phi) = - \\frac{1}{\\alpha} "
                                              "\\delta \\ln \\rho(r,\\theta,\\phi)$. $\\xi$ is the `vs to "
                                              "density scaling' parameter and $\\alpha$ is the 'Thermal "
                                              "expansion coefficient in initial temperature scaling' "
                                              "parameter. The temperature perturbation is added to an "
                                              "otherwise constant temperature (incompressible model) or "
                                              "adiabatic reference profile (compressible model).If a depth "
                                              "is specified in 'Remove temperature heterogeneity down to "
                                              "specified depth', there is no temperature perturbation "
                                              "prescribed down to that depth."
                                              "\n"
                                              "Note the required file format if the vs to density scaling is read in "
                                              "from a file: The first lines may contain any number of comments "
                                              "if they begin with '#', but one of these lines needs to "
                                              "contain the number of points in the reference state as "
                                              "for example '# POINTS: 3'. "
                                              "Following the comment lines there has to be a single line "
                                              "containing the names of all data columns, separated by arbitrarily "
                                              "many spaces. Column names are not allowed to contain spaces. "
                                              "The file can contain unnecessary columns, but for this plugin it "
                                              "needs to at least provide the columns named `depth' and "
                                              "`vs\\_to\\_density'. "
                                              "Note that the data lines in the file need to be sorted in order "
                                              "of increasing depth from 0 to the maximal depth in the model "
                                              "domain. Points in the model that are outside of the provided "
                                              "depth range will be assigned the maximum or minimum depth values, "
                                              "respectively. Points do not need to be equidistant, "
                                              "but the computation of properties is optimized in speed "
                                              "if they are.")
  }
}
