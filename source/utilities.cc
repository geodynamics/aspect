/*
  Copyright (C) 2014 - 2015 by the authors of the ASPECT code.

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

#include <aspect/global.h>
#include <aspect/utilities.h>

#include <deal.II/base/std_cxx1x/array.h>
#include <deal.II/base/point.h>

#include <deal.II/base/table.h>
#include <deal.II/base/table_indices.h>
#include <deal.II/base/function_lib.h>

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




    bool
    fexists(const std::string &filename)
    {
      std::ifstream ifile(filename.c_str());
      return static_cast<bool>(ifile); // only in c++11 you can convert to bool directly
    }

    template <int dim>
    AsciiDataLookup<dim>::AsciiDataLookup(const unsigned int components,
                                          const double scale_factor)
        :
        components(components),
        data(components),
        scale_factor(scale_factor)
      {}

    template <int dim>
    void
    AsciiDataLookup<dim>::load_file(const std::string &filename)
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
                       ExcMessage (std::string("Couldn't open data file <"
                           +
                           filename
                           +
                           ">.")));

          // Read header lines and if necessary reinit tables
          while (in.peek() == '#')
            {
              std::string line;
              getline(in,line);
              std::stringstream linestream(line);
              std::string word;
              while (linestream >> word)
                if (word == "POINTS:")
                  for (unsigned int i = 0; i < dim; i++)
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
          Table<dim,double> data_table;
          data_table.TableBase<dim,double>::reinit(table_points);
          std::vector<Table<dim,double> > data_tables(components+dim,data_table);

          // Read data lines
          unsigned int i = 0;
          double temp_data;

          while (in >> temp_data)
            {
              const unsigned int column_num = i%(components+dim);

              if (column_num >= dim)
                temp_data *= scale_factor;

              data_tables[column_num](compute_table_indices(i)) = temp_data;

              i++;

              // TODO: add checks for coordinate ordering in data files
            }

          AssertThrow(i == (components + dim) * data_table.n_elements(),
                      ExcMessage (std::string("Number of read in points does not match number of expected points. File corrupted?")));

          std_cxx11::array<unsigned int,dim> table_intervals;

          for (unsigned int i = 0; i < dim; i++)
            {
              table_intervals[i] = table_points[i] - 1;

              TableIndices<dim> idx;
              grid_extent[i].first = data_tables[i](idx);
              idx[i] = table_points[i] - 1;
              grid_extent[i].second = data_tables[i](idx);
            }

          for (unsigned int i = 0; i < components; i++)
            {
              if (data[i])
                delete data[i];
              data[i] = new Functions::InterpolatedUniformGridData<dim> (grid_extent,
                                                                         table_intervals,
                                                                         data_tables[dim+i]);
            }
        }


    template <int dim>
    double
    AsciiDataLookup<dim>::get_data(const Point<dim> &position,
                                   const unsigned int component) const
        {
          return data[component]->value(position);
        }


    template <int dim>
    TableIndices<dim>
    AsciiDataLookup<dim>::compute_table_indices(const unsigned int i) const
    {
      TableIndices<dim> idx;
      idx[0] = (i / (components+dim)) % table_points[0];
      if (dim >= 2)
        idx[1] = ((i / (components+dim)) / table_points[0]) % table_points[1];
      if (dim == 3)
        idx[2] = (i / (components+dim)) / (table_points[0] * table_points[1]);

      return idx;
    }



        template <int dim>
        AsciiDataBase<dim>::AsciiDataBase ()
        {}


        template <int dim>
        void
        AsciiDataBase<dim>::declare_parameters (ParameterHandler  &prm,
                                                const std::string &default_directory,
                                                const std::string &default_filename)
        {
          prm.enter_subsection ("Ascii data model");
          {
            prm.declare_entry ("Data directory",
                               default_directory,
                               Patterns::DirectoryName (),
                               "The name of a directory that contains the model data. This path "
                               "may either be absolute (if starting with a '/') or relative to "
                               "the current directory. The path may also include the special "
                               "text '$ASPECT_SOURCE_DIR' which will be interpreted as the path "
                               "in which the ASPECT source files were located when ASPECT was "
                               "compiled. This interpretation allows, for example, to reference "
                               "files located in the 'data/' subdirectory of ASPECT. ");
            prm.declare_entry ("Data file name",
                               default_filename,
                               Patterns::Anything (),
                               "The file name of the material data. Provide file in format: "
                               "(Velocity file name).\\%s%d where \\\\%s is a string specifying "
                               "the boundary of the model according to the names of the boundary "
                               "indicators (of a box or a spherical shell).%d is any sprintf integer "
                               "qualifier, specifying the format of the current file number. ");
            prm.declare_entry ("Scale factor", "1",
                               Patterns::Double (0),
                               "Scalar factor, which is applied to the boundary velocity. "
                               "You might want to use this to scale the velocities to a "
                               "reference model (e.g. with free-slip boundary) or another "
                               "plate reconstruction. Another way to use this factor is to "
                               "convert units of the input files. The unit is assumed to be"
                               "m/s or m/yr depending on the 'Use years in output instead of "
                               "seconds' flag. If you provide velocities in cm/yr set this "
                               "factor to 0.01.");
          }
          prm.leave_subsection();
        }


        template <int dim>
        void
        AsciiDataBase<dim>::parse_parameters (ParameterHandler &prm)
        {
          prm.enter_subsection("Ascii data model");
          {
            // Get the path to the data files. If it contains a reference
            // to $ASPECT_SOURCE_DIR, replace it by what CMake has given us
            // as a #define
            data_directory    = prm.get ("Data directory");
            {
              const std::string      subst_text = "$ASPECT_SOURCE_DIR";
              std::string::size_type position;
              while (position = data_directory.find (subst_text),  position!=std::string::npos)
                data_directory.replace (data_directory.begin()+position,
                                        data_directory.begin()+position+subst_text.size(),
                                        ASPECT_SOURCE_DIR);
            }

            data_file_name    = prm.get ("Data file name");
            scale_factor      = prm.get_double ("Scale factor");
          }
          prm.leave_subsection();
        }

        // Explicit instantiations
        template class AsciiDataLookup<1>;
        template class AsciiDataLookup<2>;
        template class AsciiDataLookup<3>;
        template class AsciiDataBase<2>;
        template class AsciiDataBase<3>;


    }
  }
