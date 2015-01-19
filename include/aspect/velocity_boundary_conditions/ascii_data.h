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


#ifndef __aspect__velocity_boundary_conditions_ascii_data_h
#define __aspect__velocity_boundary_conditions_ascii_data_h

#include <aspect/velocity_boundary_conditions/interface.h>
#include <aspect/simulator_access.h>

#include <aspect/geometry_model/box.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/utilities.h>

#include <deal.II/base/table.h>
#include <deal.II/base/table_indices.h>
#include <deal.II/base/function_lib.h>

#include <fstream>


namespace aspect
{
  namespace VelocityBoundaryConditions
  {
    using namespace dealii;

    namespace internal
     {
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
           AsciiDataLookup(const std_cxx11::array<unsigned int,3> &npoints,
                           const GeometryModel::Interface<dim> &geometry_model,
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
                 table_points[i] = npoints[boundary_dimensions[i]];
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


    /**
     * A class that implements prescribed velocity boundary conditions
     * determined from a AsciiData input file.
     *
     * @ingroup VelocityBoundaryConditionsModels
     */
    template <int dim>
    class AsciiData : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Empty Constructor.
         */
        AsciiData ();

        /**
         * Initialization function. This function is called once at the
         * beginning of the program. Checks preconditions.
         */
        virtual
        void
        initialize ();

        /**
         * A function that is called at the beginning of each time step. For
         * the current plugin, this function loads the next data files if
         * necessary and outputs a warning if the end of the set of data
         * files is reached.
         */
        virtual
        void
        update ();

        /**
         * Return the boundary velocity as a function of position. For the
         * current class, this function returns value from the text files.
         */
        Tensor<1,dim>
        boundary_velocity (const Point<dim> &position) const;


        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm);

      private:
        /**
         * A variable that stores the currently used data file of a
         * series. It gets updated if necessary by set_current_time.
         */
        int  current_file_number;

        /**
         * Time from which on the data file with number 'First data
         * file number' is used as boundary condition. Previous to this
         * time, a no-slip boundary condition is assumed. Depending on the setting
         * of the global 'Use years in output instead of seconds' flag
         * in the input file, this number is either interpreted as seconds or as years."
         */
        double first_data_file_model_time;

        /**
         * Number of the first data file to be loaded when the model time
         * is larger than 'First data file model time'.
         */
        int first_data_file_number;

        /**
         * In some cases the boundary files are not numbered in increasing
         * but in decreasing order (e.g. 'Ma BP'). If this flag is set to
         * 'True' the plugin will first load the file with the number
         * 'First data file number' and decrease the file number during
         * the model run.
         */
        bool decreasing_file_order;

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
         * Boundary id for which this plugin is created. Is initialized in
         * initialize() and used to determine the appropriate file name.
         */
        types::boundary_id boundary_id;

        /**
         * Time in model units (depends on other model inputs) between two
         * data files.
         */
        double data_file_time_step;

        /**
         * Number of grid points in data file
         */
        std_cxx1x::array<unsigned int,3> data_points;

        /**
         * Weight between data file n and n+1 while the current time is
         * between the two values t(n) and t(n+1).
         */
        double time_weight;

        /**
         * State whether we have time_dependent boundary conditions. Switched
         * off after finding no more data files to suppress attempts to
         * read in new files.
         */
        bool time_dependent;

        /**
         * Scale the boundary condition by a scalar factor. Can be
         * used to transform the unit of the velocities (if they are not
         * specified in the default unit (m/s or m/yr depending on the
         * "Use years in output instead of seconds" parameter).
         */
        double scale_factor;

        /**
         * Pointer to an object that reads and processes data we get from
         * text files.
         */
        std_cxx11::shared_ptr<internal::AsciiDataLookup<dim,dim-1> > lookup;

        /**
         * Handles the update of the data in lookup.
         */
        void
        update_data ();

        /**
         * Handles settings and user notification in case the time-dependent
         * part of the boundary condition is over.
         */
        void
        end_time_dependence (const int timestep);

        /**
         * Create a filename out of the name template.
         */
        std::string
        create_filename (const int timestep) const;
    };
  }
}


#endif
