/*
  Copyright (C) 2017 - 2019 by the authors of the ASPECT code.

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


#include <aspect/initial_composition/crust1.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/material_model/crust1.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/utilities.h>
#include <aspect/simulator_access.h>
#include <fstream>
#include <iostream>
#include <array>

namespace aspect
{
  namespace InitialComposition
  {

    template <int dim>
    double
    Crust1<dim>::
    initial_composition (const Point<dim> &position,
                         const unsigned int compositional_index) const
    {
      AssertThrow(this->introspection().compositional_name_exists("crust1_density"),
                  ExcMessage("The initial composition plugin `crust1' did not find a "
                             "compositional field called `crust1_density' to initialize. Please add a "
                             "compositional field with this name."));

      AssertThrow(dynamic_cast<const GeometryModel::SphericalShell<dim>*> (&this->get_geometry_model()) != nullptr,
                  ExcMessage("The initial composition plugin `crust1' only works with "
                             "the 'SphericalShell' geometry model."));      

      AssertThrow(dynamic_cast<const MaterialModel::Crust1<dim>*> (&this->get_material_model()) != nullptr,
                  ExcMessage("The initial composition plugin `crust1' only works with "
                             "the 'Crust1' material model."));

      const unsigned int crust1_index = this->introspection().compositional_index_for_name("crust1_density");
      if (compositional_index == crust1_index)
        {
          
          // Get coordinates (surface elevation, longitude (0-360), latitude (-90 - 90))
          std::array<double,dim> wcoord = Utilities::Coordinates::WGS84_coordinates(position);
                
          // Define location of current point
          const double lon    = wcoord[1];       
          const double lat    = wcoord[2];
          const double colat  = 90. - lat;

          // The reference radius is hardwired in here and needs to be
          // an input parameter.
          const double depth = (position.norm() - 6371000.)/1000.;

          int ilat = std::min(int(colat),179);
          int ilon = std::min(int(lon),359);

          int cr1_index = (ilat+1)*360 + ilon;
          double cr1_rho = 0.;
          if (depth < bnds[cr1_index][8])
            {
              cr1_rho = 3.3;
            }
          else if (depth < bnds[cr1_index][7])
            {
              cr1_rho = rho[cr1_index][7]; 
            }
          else if (depth < bnds[cr1_index][6])
            {
              cr1_rho = rho[cr1_index][6];
            } 
          else if (depth < bnds[cr1_index][5])
            {
              cr1_rho = rho[cr1_index][5];
            }
          else if (depth < bnds[cr1_index][4])
            {
              cr1_rho = rho[cr1_index][4];
            }
          else if (depth < bnds[cr1_index][3])
            {
              cr1_rho = rho[cr1_index][3];
            }
          else if (depth < bnds[cr1_index][2])
            {
              cr1_rho = rho[cr1_index][2];
            }
          else if (depth < bnds[cr1_index][1])
            {
              cr1_rho = rho[cr1_index][1];
            }
          else 
            {
              cr1_rho = 0;
            }
           return cr1_rho*1000.;
        }
      return 0.0;
    }

    template <int dim>
    void
    Crust1<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial composition model");
      {
        prm.enter_subsection("Crust1");
        {
          prm.declare_entry ("Data directory", "$ASPECT_SOURCE_DIR/data/initial-composition/crust1/",
                             Patterns::DirectoryName (),
                             "The path to the model data. ");
          prm.declare_entry ("Crust1 boundaries file name", "crust1.bnds",
                             Patterns::Anything(),
                             "The file name of the CRUST 1.0 file containing layer boundaries.");
          prm.declare_entry ("Crust1 densities file name", "crust1.rho",
                             Patterns::Anything(),
                             "The file name of the CRUST 1.0 file containing layer densities.");
	}
        prm.leave_subsection ();
      } 
      prm.leave_subsection ();
    }

    template <int dim>
    void
    Crust1<dim>::parse_parameters (ParameterHandler &prm)
    {
      AssertThrow (dim == 3,
                   ExcMessage ("The 'Crust 1' model for the initial composition "
                               "is only available for 3d computations."));

      prm.enter_subsection ("Initial composition model");
      {
        prm.enter_subsection("Crust1");
        {
          data_directory = Utilities::expand_ASPECT_SOURCE_DIR(prm.get ("Data directory"));

          if ((data_directory.size() > 0) && (data_directory[data_directory.size()-1] != '/'))
            data_directory += "/";

          crust1_bnds_file_name = prm.get ("Crust1 boundaries file name");
          crust1_rho_file_name = prm.get ("Crust1 densities file name");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();

      const int num_layers = 9;

      const std::string bnds_filename = data_directory+crust1_bnds_file_name;

      // Read CRUST1.0 boundary layer data from file
      std::istringstream inbnds(Utilities::read_and_distribute_file_content(bnds_filename, this->get_mpi_communicator()));

      std::string bnds_line;
      for (unsigned int i=0; i<bnds.size(); ++i)
        {
          std::getline(inbnds,bnds_line);
          std::istringstream issbnds(bnds_line);
          for (unsigned int j=0; j<num_layers; ++j)
            issbnds >> bnds[i][j];
        }

      const std::string rho_filename = data_directory+crust1_rho_file_name;
      
      // Read CRUST1.0 density data from file
      std::istringstream inrho(Utilities::read_and_distribute_file_content(rho_filename, this->get_mpi_communicator()));
      
      std::string rho_line;
      for (unsigned int i=0; i<rho.size(); ++i) 
        { 
          std::getline(inrho,rho_line);
          std::istringstream issrho(rho_line);
          for (unsigned int j=0; j<num_layers; ++j)
            issrho >> rho[i][j]; 
        }

      // Fill variable with lat-lon coordinates for CRUST1.0. 
      unsigned int count = 0;
      const unsigned int lon_pts_cr1 = 360;
      const unsigned int lat_pts_cr1 = 180;
      const double lon_min_cr1 = -179.5;
      const double lat_min_cr1 = -89.5;
      for (unsigned int i=0; i<lat_pts_cr1; ++i)
        {
          for (unsigned int j=0; j<lon_pts_cr1; ++j)
            {
              crust1_lon_lat[count][0] = lon_min_cr1 + double(j);
              crust1_lon_lat[count][1] = lat_min_cr1 + double(i);
              count = count + 1;
             }
        }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace InitialComposition
  {
    ASPECT_REGISTER_INITIAL_COMPOSITION_MODEL(Crust1,
                                              "crust1",
                                              "A class that implements an initial compositional field containing "
                                              "densities from the CRUST 1.0 data set. Note that this plugin only "
                                              "works if there is a compositional field called 'crust1_densities'.")
  }
}
