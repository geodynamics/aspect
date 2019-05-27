/*
  Copyright (C) 2011 - 2018 by the authors of the ASPECT code.

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
    template <int dim>
    LithosphereMask<dim>::LithosphereMask ()
      :
      lab_depths(1, 1.0)
    {}


    template <int dim>
    void
    LithosphereMask<dim>::initialize ()
    {
      if (LAB_depth_source == File)
        {
          const std::string filename = data_directory+LAB_file_name;

          std::cout << "   Loading Ascii data lookup file " << filename << "." << std::endl;

          lab_depths.load_file(filename,this->get_mpi_communicator());
        }
    }

    //Read in ascii data - two dimensions so the third column is treated as data
    template <int dim>
    double
    LithosphereMask<dim>::ascii_lab (const Point<2> &position) const
    {
      const double lab = lab_depths.get_data(position,0)*1000; //In km
      return lab;
    }

    template <int dim>
    double
    LithosphereMask<dim>::initial_temperature (const Point<dim> &position) const
    {
      double temperature;
      double depth = this->SimulatorAccess<dim>::get_geometry_model().depth(position);

      if (LAB_depth_source == File)
        {
          //Get spherical coordinates for model
          std::array<double,dim> scoord      = Utilities::Coordinates::cartesian_to_spherical_coordinates(position);
          const double phi = scoord[1];
          const double theta = scoord[2];
          const Point<2> phi_theta (phi, theta);

          //Get lab depth for specific phi and theta
          const double lab_depth = ascii_lab(phi_theta);

          if (depth <= lab_depth)
            temperature = lithosphere_temperature;
          else
            temperature = std::numeric_limits<double>::quiet_NaN();
        }

      else if (LAB_depth_source == Value)
        {
          if (depth <= max_depth*1000)//in km
            temperature = lithosphere_temperature;
          else
            temperature = std::numeric_limits<double>::quiet_NaN();
        }

      else
        {
          Assert( false, ExcMessage("Invalid method for depth specification method") );
          return 0.0;
        }

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
          prm.declare_entry ("Lithosphere temperature", "1600",
                             Patterns::Double (0),
                             "The initial temperature within lithosphere, applied above"
                             "the maximum lithosphere depth.");
          prm.declare_entry ("Depth specification method", "Value",
                             Patterns::Selection("File|Value"),
                             "Method that is used to specify the depth of the lithosphere-asthenosphere boundary.");
          prm.declare_entry ("Maximum lithosphere depth", "200.0",
                             Patterns::Double (0),
                             "The maximum depth of the lithosphere. The model will be "
                             "NaNs below this depth.");
          prm.declare_entry ("Data directory", "$ASPECT_SOURCE_DIR/data/initial-temperature/lithosphere_mask/",
                             Patterns::DirectoryName (),
                             "The path to the LAB depth data file");
          prm.declare_entry ("LAB depth filename",
                             "LAB_CAM2016.txt",
                             Patterns::FileName (),
                             "File from which the lithosphere-asthenosphere boundary depth data is read.");
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
          lithosphere_temperature   = prm.get_double ("Lithosphere temperature");

          if ( prm.get("Depth specification method") == "File" )
            {
              LAB_depth_source = File;
              data_directory = Utilities::expand_ASPECT_SOURCE_DIR (prm.get("Data directory"));
              LAB_file_name   = prm.get("LAB depth filename");
            }
          else if ( prm.get("Depth specification method") == "Value" )
            {
              LAB_depth_source = Value;
              max_depth           = prm.get_double ("Maximum lithosphere depth");
            }



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
    ASPECT_REGISTER_INITIAL_TEMPERATURE_MODEL(LithosphereMask,
                                              "lithosphere mask",
                                              "Implementation of a model in which the initial "
                                              "temperature is set to a specified lithosphere temperature above the "
                                              "above the lithsphere-asthenosphere boundary (specified by an ascii file "
                                              "or maximum lithosphere depth value). Below this the initial temperature is set as "
                                              "NaN.  Note the required format of the input data file: The first lines may "
                                              "contain any number of comments if they begin with ‘#’, but one of these lines "
                                              "needs to contain the number of grid points in each dimension as for example"
                                              "‘# POINTS: 3 3’. For a spherical model, the order of the data columns has to be"
                                              "'phi', 'theta','depth (km)', where phi is the  azimuth angle and theta is the "
                                              "polar angle measured positive from the north pole. This plug-in can be combined "
                                              "with another using the 'replace if valid' operator. ")
  }
}
