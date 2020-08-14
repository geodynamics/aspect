/*
  Copyright (C) 2011 - 2020 by the authors of the ASPECT code.

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


#include <aspect/initial_temperature/lithosphere_mask.h>
#include <aspect/utilities.h>
#include <limits>
#include <aspect/geometry_model/interface.h>


namespace aspect
{
  namespace InitialTemperature
  {
    namespace LABDepth
    {
      template <int dim>
      LABDepthLookup<dim>::LABDepthLookup ()
        :
        lab_depths(1, 1.0)
      {}



      template <int dim>
      void
      LABDepthLookup<dim>::initialize ()
      {
        if (LAB_depth_source == File)
          {
            const std::string filename = data_directory+LAB_file_name;
            this->get_pcout() << "   Loading Ascii data lookup file " << filename << "." << std::endl;

            lab_depths.load_file(filename,this->get_mpi_communicator());
          }
      }



      template <int dim>
      double
      LABDepthLookup<dim>::get_lab_depth (const Point<dim> &position) const
      {
        if (LAB_depth_source == File)
          {
            //Get spherical coordinates for model
            Assert (dim == 3, ExcNotImplemented());
            std::array<double,dim> scoord      = Utilities::Coordinates::cartesian_to_spherical_coordinates<dim>(position);
            const double phi = scoord[1];
            const double theta = scoord[2 % dim]; // work-around to compile without warnings for dim==2
            const Point<2> phi_theta (phi, theta);

            //Get lab depth for specific phi and theta
            const double lab_depth = lab_depths.get_data(phi_theta,0);
            return lab_depth;
          }

        else if (LAB_depth_source == Value)
          {
            return max_depth;
          }
        else
          {
            Assert( false, ExcMessage("Invalid method for depth specification method") );
            return 0.0;
          }

        return 0.0;
      }



      template <int dim>
      void
      LABDepthLookup<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Depth specification method", "Value",
                           Patterns::Selection("File|Value"),
                           "Method that is used to specify the depth of the lithosphere-asthenosphere boundary.");
        prm.declare_entry ("Maximum lithosphere depth", "200000.0",
                           Patterns::Double (0.),"Units: \\si{\\meter}."
                           "The maximum depth of the lithosphere. The model will be "
                           "NaNs below this depth.");
        prm.declare_entry ("Data directory", "$ASPECT_SOURCE_DIR/data/initial-temperature/lithosphere-mask/",
                           Patterns::DirectoryName (),
                           "The path to the LAB depth data file");
        prm.declare_entry ("LAB depth filename",
                           "LAB_CAM2016.txt",
                           Patterns::FileName (),
                           "File from which the lithosphere-asthenosphere boundary depth data is read.");
      }

      template <int dim>
      void
      LABDepthLookup<dim>::parse_parameters (ParameterHandler &prm)
      {
        if ( prm.get("Depth specification method") == "File" )
          {
            LAB_depth_source = File;
            data_directory = Utilities::expand_ASPECT_SOURCE_DIR (prm.get("Data directory"));
            LAB_file_name = prm.get("LAB depth filename");
          }
        else if ( prm.get("Depth specification method") == "Value" )
          {
            LAB_depth_source = Value;
            max_depth = prm.get_double ("Maximum lithosphere depth");
          }
      }
    }



    template <int dim>
    void
    LithosphereMask<dim>::initialize ()
    {
      lab_depth_lookup.initialize();
    }



    template <int dim>
    double
    LithosphereMask<dim>::initial_temperature (const Point<dim> &position) const
    {
      double temperature;
      const double depth = this->SimulatorAccess<dim>::get_geometry_model().depth(position);
      const double lab_depth = lab_depth_lookup.get_lab_depth(position);

      if (depth <= lab_depth)
        temperature = lithosphere_temperature;
      else
        temperature = std::numeric_limits<double>::quiet_NaN();

      return temperature;
    }

    template <int dim>
    void
    LithosphereMask<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial temperature model");
      {
        prm.enter_subsection ("Lithosphere Mask");
        {
          LABDepth::LABDepthLookup<dim>::declare_parameters(prm);

          prm.declare_entry ("Lithosphere temperature", "1600.",
                             Patterns::Double (0.),
                             "The initial temperature within lithosphere, applied above"
                             "the maximum lithosphere depth.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    LithosphereMask<dim>::parse_parameters (ParameterHandler &prm)
    {
      AssertThrow (dim == 3,
                   ExcMessage ("The 'Lithosphere mask' model for the initial "
                               "temperature is only available for 3d computations."));

      prm.enter_subsection("Initial temperature model");
      {
        prm.enter_subsection ("Lithosphere Mask");
        {
          lab_depth_lookup.initialize_simulator(this->get_simulator());
          lab_depth_lookup.parse_parameters(prm);
          lithosphere_temperature = prm.get_double ("Lithosphere temperature");
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
    template class LABDepth::LABDepthLookup<2>;
    template class LABDepth::LABDepthLookup<3>;


    ASPECT_REGISTER_INITIAL_TEMPERATURE_MODEL(LithosphereMask,
                                              "lithosphere mask",
                                              "Implementation of a model in which the initial "
                                              "temperature is set to a specified lithosphere temperature above the "
                                              "lithosphere-asthenosphere boundary (specified by an ascii file "
                                              "or maximum lithosphere depth value). Below this the initial temperature is set as "
                                              "NaN.  Note the required format of the input data file: The first lines may "
                                              "contain any number of comments if they begin with '#', but one of these lines "
                                              "needs to contain the number of grid points in each dimension as for example "
                                              "'# POINTS: 3 3'. For a spherical model, the order of the data columns has to be "
                                              "'phi', 'theta', 'depth (m)', where phi is the azimuth angle and theta is the "
                                              "polar angle measured positive from the north pole. This plug-in can be combined "
                                              "with another using the 'replace if valid' operator. ")
  }
}
