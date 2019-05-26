/*
  Copyright (C) 2018 by the authors of the ASPECT code.

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

#include "common.h"
#include <aspect/utilities.h>
#include <aspect/initial_composition/interface.h>
#include <aspect/material_model/utilities.h>  // import the perplex and hefesto readers

#include <aspect/simulator_access.h>

#include <deal.II/numerics/kdtree.h>
#include <deal.II/base/point.h>
#include <deal.II/base/mpi.h>


/**
 * Some header functions and constructors:
 */
//------------------------------------------------------------------------------------
/**
 * Polynominal interpolation routines from:
 * Numerical Recipes. The Art of Scientific Computing. 3rd ed. Press et. al.
 *
 * Each algorithm is marked with the section number it is found in.
 * Comments are completely preserved.
 */
/*



//----------------------------------------------------------
/**
 *
 */
namespace aspect
{
  template <int dim> class SimulatorAccess;
  namespace Utilities
  {
    using namespace dealii;


    /**
     * UnstructuredLookup:
     *
     * This utility returns data values at requested coordinate
     * points starting from an n x m dimensional data table. The
     * KDTree functionality from dealii is used to query the closest
     * coordinates from the table to the requested coordinate,
     * and can do one of the following:
     *
     * - Return the specified number of closest coordinates and
     *   their associated data values.
     * - Return interpolated data values at the requested
     *   coordinate from the specified number of closest
     *   coordinates returned by KDTree.
     *
     * Interpolation routine must be specified by the user,
     * and be in the list of options available.
     */
    template <int querydim, int datadim>
    class UnstructuredLookup
    {
        /** Things and steps this function should do:
         * ------------------------------------------
         * -- Load the necessary n-dim coordinate and m-dim data from a table
         * -- Read in a target point that needs to be assessed
         * -- Read in the interpolation type to find the corresponding
         * -- Setup the KDTree from dealii
         *    -- Load only coordinate points into a KDTree
         *    -- Use the KDTree to return the n-number of closest points to a target point using:
         *    -- KDTree< dim >::get_closest_points (const Point< dim > &target, const unsigned int n_points) const
         *    data that goes with the target point
         * -- Return data values for use
         */


        //void set_data (std::vector<Point<querydim>> coordinates, vector of std::array<double,datadim> data)
        // from dealii::numerics::KDTree ...
        // void KDTree< dim >::set_points (const std::vector< Point< dim > > & pts)

        //void load_ascii_data (filename) ....call set_data

        //const std::array<double,datadim> & get_closest_point_data (Point<querydim> point_in)

        /*
        // A structure to contain the coordinate points in
        struct Coordinate
        {
          double x, y;
          Coordinate(double paramx, double paramy) : x(paramx), y(paramy) {}
        };
        */
        /*
        //std::vector<double> coords;
        const std::array<double,2> coords;
        std::vector<Point<2> > coordinate_points;
        std::vector<std::array<double,6> > data;

        std::double data1;
        std::double data2;
        std::double data3;
        std::double data4;
        std::double data5;
        std::double data6;

        std::vector<std::string> keys = {"T(K)","P(bar)","rho,kg/m3","alpha,1/K","cp,J/K/kg","vp,km/s","vs,km/s","h,J/kg"};
        */

        /*
          void set_data (std::vector<Point<2> > coordinates, std::vector<std::array<double,6> > data)
          {

          }
        */

      public:
        void load_ascii_data (const MPI_Comm &comm); // Load the data table into memory
        virtual void initialize ();
        void setup_kdtree (std::vector<std_cxx14::make_unique<double> > material_looup);
        //void closest_points; // Setup the KDTree and return the specified number of closest points
        //void data_interpolation; // Interpolate the data values based on data at the closest points


        /**
         * The following variables are properties of the material files
         * we read in.
         */
        std::string datadirectory;
        std::vector<std::string> material_file_names;
        std::vector<std::string> derivatives_file_names;
        unsigned int n_material_data;
        //bool use_table_properties;
        //bool use_enthalpy;
        bool use_bilinear_interpolation;
        std::vector<std_cxx14::make_unique<double> > material_lookup;


        /**
         * The format of the provided material files. Currently we support
         * the PERPLEX and HeFESTo data formats.
         */
        enum formats
        {
          perplex,
          hefesto
        } material_file_format;
    };


    template <int querydim, int datadim>
    void
    UnstructuredLookup<querydim,datadim>::initialize()
    {
      // !!! Need to add the file number to the front of the vectors !!!
      n_material_data = material_file_names.size();
      for (unsigned i = 0; i < n_material_data; i++)
        {
          if (material_file_format == perplex)
            material_lookup
            .push_back(std_cxx14::make_unique<MaterialModel::MaterialUtilities::Lookup::PerplexReader>(datadirectory+material_file_names[i],
                       use_bilinear_interpolation,
                       this->get_mpi_communicator()));
          else if (material_file_format == hefesto)
            material_lookup
            .push_back(std_cxx14::make_unique<MaterialModel::MaterialUtilities::Lookup::HeFESToReader>(datadirectory+material_file_names[i],
                       datadirectory+derivatives_file_names[i],
                       use_bilinear_interpolation,
                       this->get_mpi_communicator()));
          else
            AssertThrow (false, ExcNotImplemented());
        }
    }


    template <int querydim, int datadim>
    void
    UnstructuredLookup<querydim,datadim>::setup_kdtree (std::vector<std_cxx14::make_unique<double> > material_lookup)
    {
      std::vector<Point<querydim> > coords;
      //std::vector<std::vector <double> > data;
      const Table<datadim, double> data;

      for (unsigned int i=0; i<material_lookup.size(); ++i)
        {
          // read the coordinate components
          Point<querydim> coordinates;

          // !!! need to fix this -- it is hard-coded for querydim=3 !!!
          coordinates = (material_lookup[i][0],material_lookup[i][1],material_lookup[i][2]);
          coords.emplace_back(coordinates);

          // read the data components
          std::vector<double> values;
          for (unsigned int j=querydim; j<datadim; ++j)
            {
              values.emplace_back(material_lookup[i][j]);
            }
          data.emplace_back (values);
        }

      KDTree<querydim> kdtree;
      kdtree.set_points (coords);

      // !!! Need to return something or wrap it all in a "shared" object" !!!
    }


    template <int querydim, int datadim>
    void
    UnstructuredLookup<querydim,datadim>::load_ascii_data (const MPI_Comm &comm)
    {


      std::vector<Point<2> > coords;
      std::vector<std::vector <double> > data;

      std::vector<std::string> keys = {"T(K)","P(bar)","rho,kg/m3","alpha,1/K","cp,J/K/kg","vp,km/s","vs,km/s","h,J/kg"};


      std::string temp;
      // Read data from disk and distribute among processes
      std::string filename = "/home/paulbremner/data/ASPECT/development/devel_hackathon_2018/aspect/data/material-model/latent-heat-enthalpy-test/testdata.txt";
      std::istringstream in(Utilities::read_and_distribute_file_content(filename,comm));

      // Output to file to check the results
      std::ofstream outputFile;
      outputFile.open("CheckTheKDTree.txt");

      outputFile << "Starting output...\n\n";


      // Try using the Perplex Reader instead
      //**********************************************************
      std::string material_file_format = "perplex";


      //**********************************************************

      // Skip header
      for (unsigned int i=0; i<12; ++i)
        {
          std::getline(in,temp); // throw these lines away
        }

      outputFile << "Check point 1...\n\n";

      std::string column_header;
      std::getline(in,column_header);
      // output header to the file
      outputFile << "Column Header: \n";
      outputFile << column_header << "\n\n";

      // !!!!!!!
      // Split column_header into separate strings to attach as reference to data values.
      // Not certain what is the best construction yet that ensures values are not mixed up somewhere else
      const std::vector<std::string> x_mapped_field_entries = dealii::Utilities::split_string_list(column_header, ' ');
      outputFile << "The Column Header as a set of keys: \n";
      for (unsigned int i=0; i<x_mapped_field_entries.size(); ++i)
        {
          outputFile << x_mapped_field_entries[i] << " ";
        }
      outputFile << "\n\n";


      ///// Maybe don't mix stream and getline from above??? Need to make this end at eof
      ///// Need correct type for line
      //
      while (getline(in,temp)) //(in.good() && !in.eof())   // ????
        {
          std::istringstream line(temp);
          // read the two components of the coordinates
          Point<2> coordinates;
          line >> coordinates;
          coords.emplace_back(coordinates);

          outputFile << "Coordinates: " << coordinates << "\n";
          outputFile << "Data: ";
          // read the 6 data columns
          std::vector<double> values;
          for (unsigned int i=0; i<6; ++i)
            {
              double value;
              line >> value;

              values.emplace_back(value);

              outputFile << value << " ";
            }
          data.emplace_back (values);

          outputFile << "\n\nCheck point 2...\n\n";

        }

      //dealii::numerics::set_data(std::vector<Point<2> > coords, std::vector<std::array<double,6> > data);
      KDTree<querydim> kdtree;
      kdtree.set_points (coords);


      //get_closet_data ();

      Point<2> target_point(940, 76000); // !!! Temporary Values. Should Be Read In !!!
      int n_points = 3; // !!! Temporary Values. Should Be Read In !!!


      //const std::array<double,datadim> & get_closest_point_data (Point<querydim> point_in)
      std::vector<std::pair<unsigned int, double> > closest_point = kdtree.get_closest_points (target_point,n_points);
      //std::vector<std::pair<unsigned int, double> > KDTree< dim >::get_closest_points (const Point< dim > &target,
      //                                        const unsigned int n_points) const

      // Match the indices to the data
      outputFile << "The point needing to be assessed:\n";
      outputFile << target_point << "\n\nThe closest points and their distances:\n";
      for (unsigned int j=0; j<closest_point.size(); ++j)
        {
          outputFile << std::get<0>(closest_point[j]) << " " << std::get<1>(closest_point[j]) << "\n";
          outputFile << "Corresponding to: (T,P) --- " << coords[std::get<0>(closest_point[j])] << "\n";
          outputFile << "                  data  --- ";
          for (unsigned int k=0; k<data[std::get<0>(closest_point[j])].size(); ++k)
            {
              outputFile << data[std::get<0>(closest_point[j])][k] << " ";
            }
          outputFile << "\n\n";
        }



      outputFile.close();


      // Match the closest points with the data.
      // Create a new vector of coordinates, distance to target_point, and associated data.
      // Return the new vector.

      for (unsigned int icoord=0; icoord<coords.size(); ++icoord)
        {
          //look for a match and record the element number
          // use that element number to paste data vector onto the new vector
        }

      //return subtable_data;

//private: dealii::KDTree<2>
    }
    /*
    template <int querydim, int datadim>
    void
    UnstructuredLookup<querydim,datadim>::data_interpolation (*VECTOR-OF-CLOSEST-POINTS*, *INTERPOLATION-CHOICE*)
    {

      // Interpolate data based on choice. New choices can be added as desired.

      // Closest Points:
      // Simply return the closest points obtained from the KDTree
      if (interp_choice == "closest_points")
      {
        cout << "closest_points\n";
        //return closest_point_data;
      }

      // Radial Basis Functions:
      // Interpolate the data using Radial Basis Functions
      else if (interp_choice == "radial_basis_functions")
      {
        cout << "radial_basis_functions\n";

        // Do stuff....


        // return data;
      }

      // Linear Interpolation using Weighted Average:
      // Interpolate the data by making a simple weighted average between data points.
      // The weights are determined from the distances of the specified number of
      // closest coordinates returned from KDTree.
      else if (interp_choice == "linear_w_weigtedavg")
      {
        cout << "linear_w_weightedavg\n";

        // Do stuff....


        // return data;
      }

      // Polynomial Inerpolation:
      // Interpolate data based on a higher order polynomial (see Numerical Recipes).
      else if (interp_choice == "polynomial_interp")
      {
        cout << "polynomial_interp\n";

        // Do stuff....

        //return data;
      }
    }
    */
//template class UnstructuredLookup<2,6>;
    TEST_CASE("Utilities::KDTree")
    {
      UnstructuredLookup<3,6> lookup;
      datadirectory = "/home/paulbremner/data/ASPECT/development/devel_hackathon_2018/aspect/data/material-model/latent-heat-enthalpy-test/";
      material_file_names = "testdata.txt"
//lookup.load_ascii_data(get_mpi_communicator());
                            lookup.load_ascii_data(MPI_COMM_WORLD);
    }
  }
}



