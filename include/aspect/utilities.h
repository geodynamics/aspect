/*
  Copyright (C) 2014 by the authors of the ASPECT code.

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


#ifndef __aspect__utilities_h
#define __aspect__utilities_h

#include <aspect/global.h>

#include <deal.II/base/std_cxx1x/array.h>
#include <deal.II/base/point.h>
#include <deal.II/base/conditional_ostream.h>

#include <deal.II/base/table_indices.h>
#include <deal.II/base/function_lib.h>

#include <aspect/geometry_model/interface.h>
#include <aspect/simulator_access.h>



namespace aspect
{
  /**
   * A namespace for utility functions that might be used in many different
   * places to prevent code duplication.
   */
  namespace Utilities
  {
    using namespace dealii;

    template <int dim>
    std_cxx1x::array<double,dim>
    spherical_coordinates(const Point<dim> &position)
    {
      std_cxx1x::array<double,dim> scoord;

      scoord[0] = position.norm(); // R
      scoord[1] = std::atan2(position(1),position(0)); // Phi
      if (scoord[1] < 0.0) scoord[1] = 2*numbers::PI + scoord[1]; // correct phi to [0,2*pi]
      if (dim==3)
        {
          if (scoord[0] > std::numeric_limits<double>::min())
            scoord[2] = std::acos(position(2)/scoord[0]);
          else
            scoord[2] = 0.0;
        }
      return scoord;
    }

    template <int dim>
    Point<dim>
    cartesian_coordinates(const std_cxx1x::array<double,dim> &scoord)
    {
      Point<dim> ccoord;

      switch (dim)
        {
          case 2:
          {
            ccoord[0] = scoord[0] * std::cos(scoord[1]); // X
            ccoord[1] = scoord[0] * std::sin(scoord[1]); // Y
            break;
          }
          case 3:
          {
            ccoord[0] = scoord[0] * std::sin(scoord[2]) * std::cos(scoord[1]); // X
            ccoord[1] = scoord[0] * std::sin(scoord[2]) * std::sin(scoord[1]); // Y
            ccoord[2] = scoord[0] * std::cos(scoord[2]); // Z
            break;
          }
          default:
            Assert (false, ExcNotImplemented());
            break;
        }

      return ccoord;
    }

    /**
     * Returns spherical coordinates of a cartesian point. The returned array
     * is filled with radius, phi and theta (polar angle). If the dimension is
     * set to 2 theta is omitted. Phi is always normalized to [0,2*pi].
     *
     */
//    template <int dim>
//    std_cxx1x::array<double,dim>
//    spherical_coordinates(const Point<dim> &position);

    /**
     * Return the cartesian point of a spherical position defined by radius,
     * phi and theta (polar angle). If the dimension is set to 2 theta is
     * omitted.
     */
//    template <int dim>
//    Point<dim>
//    cartesian_coordinates(const std_cxx1x::array<double,dim> &scoord);

    /**
     * Checks whether a file named filename exists.
     *
     * @param filename File to check existence
     */
    bool fexists(const std::string &filename);



    template <int dim>
    class AsciiDataLookup
    {
      public:
        AsciiDataLookup(const unsigned int components,
                        const double scale_factor);

        /**
         * Loads a data text file. Throws an exception if the file does not exist,
         * if the data file format is incorrect or if the file grid changes over model runtime.
         */
        void
        load_file(const std::string &filename);

        /**
         * Returns the computed data (velocity, temperature, etc. - according to the used plugin)
         * in cartesian coordinates.
         *
         * @param position The current position to compute the data (velocity, temperature, etc.)
         * @param time_weight A weighting between the two current timesteps n and n+1
         */
        double
        get_data(const Point<dim> &position,
                 const unsigned int component) const;

      private:
        /**
         * The number of data components read in (=columns in the data file).
         */
        const unsigned int components;

        /**
         * Interpolation functions to access the data.
         */
        std::vector<Functions::InterpolatedUniformGridData<dim> *> data;

        /**
         * Model size
         */
        std_cxx11::array<std::pair<double,double>,dim> grid_extent;

        /**
         * Number of points in the data grid.
         */
        TableIndices<dim> table_points;

        /**
         * Scales the data boundary condition by a scalar factor. Can be
         * used to transform the unit of the data.
         */
        const double scale_factor;

        /**
         * Computes the table indices of each entry in the input data file.
         * The index depends on dim, grid_dim and the number of components.
         */
        TableIndices<dim>
        compute_table_indices(const unsigned int i) const;

    };

    template <int dim>
    class AsciiDataBase : public SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor
         */
        AsciiDataBase();

        /**
         * Declare the parameters all derived classes take from input files.
         */
        static
        void
        declare_parameters (ParameterHandler  &prm,
                            const std::string &default_directory,
                            const std::string &default_filename);

        /**
         * Read the parameters from the parameter file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

        /**
         * Directory in which the data files are present.
         */
        std::string data_directory;

        /**
         * Filename of data file. The file names can contain
         * the specifiers %s and/or %c (in this order), meaning the name of the
         * boundary and the number of the data file time step.
         */
        std::string data_file_name;

        /**
         * Scale the data by a scalar factor. Can be
         * used to transform the unit of the data (if they are not
         * specified in SI units (m/s or m/yr depending on the
         * "Use years in output instead of seconds" parameter).
         */
        double scale_factor;
    };

  }
}

#endif
