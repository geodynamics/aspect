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

#include <deal.II/base/table.h>
#include <deal.II/base/table_indices.h>
#include <deal.II/base/function_lib.h>

#include <aspect/geometry_model/box.h>
#include <aspect/geometry_model/spherical_shell.h>

#include <fstream>

namespace aspect
{
  /**
   * A namespace for utility functions that might be used in many different
   * places to prevent code duplication.
   */
  namespace Utilities
  {
    using namespace dealii;

    /**
     * Returns spherical coordinates of a cartesian point. The returned array
     * is filled with radius, phi and theta (polar angle). If the dimension is
     * set to 2 theta is omitted. Phi is always normalized to [0,2*pi].
     *
     */
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

    /**
     * Return the cartesian point of a spherical position defined by radius,
     * phi and theta (polar angle). If the dimension is set to 2 theta is
     * omitted.
     */
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

    namespace
    {
      bool
      fexists(const std::string &filename)
      {
        std::ifstream ifile(filename.c_str());
        return !(!ifile); // only in c++11 you can convert to bool directly
      }
    }

    template <int dim, int grid_dim>
    class AsciiDataLookup
    {
      public:
        AsciiDataLookup(const GeometryModel::Interface<dim> &geometry_model,
                        const unsigned int components,
                        const double scale_factor,
                        const types::boundary_id boundary_id = numbers::invalid_boundary_id)
      :
        components(components),
        geometry_model(geometry_model),
        data(components),
        old_data(components),
        scale_factor(scale_factor)
      {
          Assert (grid_dim <= dim,
              ExcMessage("You can only create an AsciiDataLookup object for a grid of at "
                  "most the number of dimensions your model uses."));

          const std_cxx11::array<unsigned int,grid_dim> dimensions = get_boundary_dimensions(boundary_id);
          const std_cxx11::array<std::pair<double,double>,dim> model_extent = get_model_extent(geometry_model);

          for (unsigned int i = 0; i < grid_dim; i++)
            {
              boundary_dimensions[i] = dimensions[i];
              grid_extent[i] = model_extent[boundary_dimensions[i]];
            }
      }

        void screen_output(const ConditionalOStream &pcout) const
        {
          std::ostringstream output;

          output << std::setprecision (3) << std::setw(3) << std::fixed << std::endl
              << "   Set up Ascii Data boundary module."  << std::endl
              << "   The grid extents to:";
          for (unsigned int i = 0; i < grid_dim; i++)
            output << " " << grid_extent[i].second;

          output << std::endl;

          pcout << output.str();
        }

        void
        load_file(const std::string &filename)
        {
          // Check whether file exists, we do not want to throw
          // an exception in case it does not, because it could be by purpose
          // (i.e. the end of the boundary condition is reached)
          AssertThrow (fexists(filename),
              ExcMessage (std::string("Ascii data file <")
          +
          filename
          +
          "> not found!"));

          std::ifstream in(filename.c_str(), std::ios::in);
          AssertThrow (in,
                       ExcMessage (std::string("Couldn't find data. Is file format correct?")));

          // Read header lines and if necessary reinit tables
          while (in.peek() == '#')
            {
              std::string line;
              getline(in,line);
              std::stringstream linestream(line);
              std::string word;
              while (linestream >> word)
                if (word == "POINTS:")
                  for (unsigned int i = 0; i < grid_dim; i++)
                    {
                      unsigned int temp_index;
                      linestream >> temp_index;

                      if (table_points[i] == 0)
                        table_points[i] = temp_index;
                      else
                        AssertThrow (table_points[i] == temp_index,
                            ExcMessage("The file grid must not change over model runtime. "
                                "Either you prescribed a conflicting number of points in "
                                "the input file, or the POINTS comment in your data files "
                                "is changing between following files."));
                    }
            }

          /**
           * Table for the new data. This peculiar reinit is necessary, because
           * there is no constructor for Table, which takes TableIndices as
           * argument.
           */
          Table<grid_dim,double> data_table;
          data_table.TableBase<grid_dim,double>::reinit(table_points);
          std::vector<Table<grid_dim,double> > data_tables(components+grid_dim,data_table);

          // Read data lines
          unsigned int i = 0;
          double temp_data;

          while (in >> temp_data)
            {
              data_tables[i%(components+grid_dim)](compute_table_indices(i)) = temp_data * scale_factor;
              i++;

              // TODO: add checks for coordinate ordering in data files
            }

          AssertThrow(i == (components + grid_dim) * data_table.n_elements(),
                      ExcMessage (std::string("Number of read in points does not match number of expected points. File corrupted?")));

          std_cxx11::array<unsigned int,grid_dim> table_intervals;
          for (unsigned int i = 0; i < grid_dim; i++)
            table_intervals[i] = table_points[i] - 1;

          for (unsigned int i = 0; i < components; i++)
            {
              if (old_data[i])
                delete old_data[i];
              old_data[i] = data[i];
              data[i] = new Functions::InterpolatedUniformGridData<grid_dim> (grid_extent,
                                                                              table_intervals,
                                                                              data_tables[grid_dim+i]);
            }
        }


        double
        get_data(const Point<dim> &position,
                 const unsigned int component,
                 const double time_weight) const
        {
          Point<dim> internal_position = position;

          if (dynamic_cast<const GeometryModel::SphericalShell<dim>*> (&geometry_model) != 0)
            {
              const std_cxx11::array<double,dim> spherical_position =
                  ::aspect::Utilities::spherical_coordinates(position);

              const GeometryModel::SphericalShell<dim> *geometry_model_spherical_shell =
                            dynamic_cast<const GeometryModel::SphericalShell<dim> *>(&geometry_model);

              internal_position[0] = spherical_position[0] - geometry_model_spherical_shell->inner_radius();

              for (unsigned int i = 1; i < dim; i++)
                internal_position[i] = spherical_position[i];
            }


          Point<grid_dim> data_position;
          for (unsigned int i = 0; i < grid_dim; i++)
            data_position[i] = internal_position[boundary_dimensions[i]];

          return time_weight * data[component]->value(data_position)
              + (1 - time_weight) * old_data[component]->value(data_position);
        }

      private:
        /**
         * The number of data components read in (=columns in the data file).
         */
        const unsigned int components;

        /**
         * A reference to the geometry model. Is needed to convert the
         * position into spherical coordinates if necessary.
         */
        const GeometryModel::Interface<dim> &geometry_model;

        /**
         * Interpolation functions to access the data.
         */
        std::vector<Functions::InterpolatedUniformGridData<grid_dim> *> data;
        std::vector<Functions::InterpolatedUniformGridData<grid_dim> *> old_data;

        /**
         * Model size
         */
        std_cxx11::array<std::pair<double,double>,grid_dim> grid_extent;

        /**
         * Number of points in the data grid.
         */
        TableIndices<grid_dim> table_points;

        /**
         * Dimensions of the boundary plane
         */
        unsigned int boundary_dimensions[grid_dim];

        /**
         * Scale the data boundary condition by a scalar factor. Can be
         * used to transform the unit of the data.
         */
        const double scale_factor;


        std_cxx11::array<std::pair<double,double>,dim>
        get_model_extent (const GeometryModel::Interface<dim> &geometry_model) const
        {
          std_cxx11::array<std::pair<double,double>,dim> data_extent;

          if (dynamic_cast<const GeometryModel::Box<dim>*> (&geometry_model) != 0)
            {
              const GeometryModel::Box<dim> *geometry_model_box = dynamic_cast<const GeometryModel::Box<dim> *>(&geometry_model);
              const Point<dim> box_origin = geometry_model_box->get_origin();
              const Point<dim> box_upper_bound = box_origin + geometry_model_box->get_extents();

              for (unsigned int i=0;i<dim;i++)
                data_extent[i] = std::make_pair<double,double> (box_origin[i],box_upper_bound[i]);
            }
          else if (dynamic_cast<const GeometryModel::SphericalShell<dim>*> (&geometry_model) != 0)
            {
              const GeometryModel::SphericalShell<dim> *geometry_model_spherical_shell =
                  dynamic_cast<const GeometryModel::SphericalShell<dim> *>(&geometry_model);

              data_extent[0] = std::make_pair(0,geometry_model_spherical_shell->maximal_depth());
              data_extent[1] = std::make_pair(0,geometry_model_spherical_shell->opening_angle() * numbers::PI / 180);
              if (dim==3)
                data_extent[2] = std::make_pair<double,double> (0,geometry_model_spherical_shell->opening_angle() * numbers::PI / 180);
            }
          else Assert (false, ExcNotImplemented());

          return data_extent;
        }

        std_cxx11::array<unsigned int,grid_dim>
        get_boundary_dimensions (const types::boundary_id boundary_id) const
        {
          std_cxx11::array<unsigned int,grid_dim> boundary_dimensions;

          if (boundary_id != numbers::invalid_boundary_id)
            {
          switch (dim)
          {
          case 2:
            if ((boundary_id == 2) || (boundary_id == 3))
              {
                boundary_dimensions[0] = 0;
              }
            else if ((boundary_id == 0) || (boundary_id == 1))
              {
                boundary_dimensions[0] = 1;
              }
            else
              AssertThrow(false,ExcNotImplemented());

            break;

          case 3:
            if ((boundary_id == 4) || (boundary_id == 5))
              {
                boundary_dimensions[0] = 0;
                boundary_dimensions[1] = 1;
              }
            else if ((boundary_id == 0) || (boundary_id == 1))
              {
                boundary_dimensions[0] = 1;
                boundary_dimensions[1] = 2;
              }
            else if ((boundary_id == 2) || (boundary_id == 3))
              {
                boundary_dimensions[0] = 0;
                boundary_dimensions[1] = 2;
              }
            else
              AssertThrow(false,ExcNotImplemented());

            break;

          default:
            AssertThrow(false,ExcNotImplemented());
          }
            }
          else
            for (unsigned int i = 0; i < grid_dim; i++)
              boundary_dimensions[i] = i;
          return boundary_dimensions;
        }


        TableIndices<grid_dim>
        compute_table_indices(const unsigned int i) const
        {
          TableIndices<grid_dim> idx;
          idx[0] = (i / (components+grid_dim)) % table_points[0];
          if (grid_dim >= 2)
            idx[1] = ((i / (components+grid_dim)) / table_points[0]) % table_points[1];
          if (grid_dim == 3)
            idx[2] = (i / (components+grid_dim)) / (table_points[0] * table_points[1]);

          return idx;
        }
    };


  }
}

#endif
