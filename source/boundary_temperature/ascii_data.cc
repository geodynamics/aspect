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


#include <aspect/global.h>
#include <aspect/boundary_temperature/ascii_data.h>

#include <deal.II/base/parameter_handler.h>

#include <aspect/geometry_model/box.h>
#include <aspect/geometry_model/spherical_shell.h>


namespace aspect
{
  namespace BoundaryTemperature
  {
    namespace internal
    {
      template <int dim>
      AsciiDataBoundary<dim>::AsciiDataBoundary ()
      :
      current_file_number(0),
      first_data_file_model_time(0.0),
      first_data_file_number(0),
      decreasing_file_order(false),
      data_file_time_step(0.0),
      time_weight(0.0),
      time_dependent(true),
      lookups(),
      old_lookups()
      {}

      template <int dim>
      void
      AsciiDataBoundary<dim>::initialize(const std::set<types::boundary_id> &boundary_ids,
          const unsigned int components)
          {
        for (typename std::set<types::boundary_id>::const_iterator
            boundary_id = boundary_ids.begin();
            boundary_id != boundary_ids.end(); ++boundary_id)
          {
            const std_cxx11::array<unsigned int,dim-1> dimensions = get_boundary_dimensions(*boundary_id);

            std_cxx11::shared_ptr<Utilities::AsciiDataLookup<dim-1> > lookup;
            lookup.reset(new Utilities::AsciiDataLookup<dim-1> (components,
                Utilities::AsciiDataBase<dim>::scale_factor));
            lookups.insert(std::make_pair(*boundary_id,lookup));

            lookup.reset(new Utilities::AsciiDataLookup<dim-1> (components,
                Utilities::AsciiDataBase<dim>::scale_factor));
            old_lookups.insert(std::make_pair(*boundary_id,lookup));


            // Set the first file number and load the first files
            current_file_number = first_data_file_number;

            const int next_file_number =
                (decreasing_file_order) ?
                    current_file_number - 1
                    :
                    current_file_number + 1;

            this->get_pcout() << std::endl << "   Loading Ascii data boundary file "
                << create_filename (current_file_number,*boundary_id) << "." << std::endl << std::endl;
            lookups.find(*boundary_id)->second->load_file(create_filename (current_file_number,*boundary_id));

            // If the boundary condition is constant, switch
            // off time_dependence immediately. This also sets time_weight to 1.0.
            // If not, also load the second file for interpolation.
            if (create_filename (current_file_number,*boundary_id) == create_filename (current_file_number+1,*boundary_id))
              {
                lookups.find(*boundary_id)->second.swap(old_lookups.find(*boundary_id)->second);
                end_time_dependence (current_file_number, *boundary_id);
              }
            else
              {
                try
                {
                    this->get_pcout() << std::endl << "   Loading Ascii data boundary file "
                        << create_filename (next_file_number,*boundary_id) << "." << std::endl << std::endl;
                    lookups.find(*boundary_id)->second.swap(old_lookups.find(*boundary_id)->second);
                    lookups.find(*boundary_id)->second->load_file(create_filename (next_file_number,*boundary_id));
                }
                catch (...)
                {
                    end_time_dependence (current_file_number, *boundary_id);
                }
              }
          }
          }


      template <int dim>
      std_cxx11::array<unsigned int,dim-1>
      AsciiDataBoundary<dim>::get_boundary_dimensions (const types::boundary_id boundary_id) const
      {
        std_cxx11::array<unsigned int,dim-1> boundary_dimensions;

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
        return boundary_dimensions;
      }

      template <int dim>
      std::string
      AsciiDataBoundary<dim>::create_filename (const int filenumber,
          const types::boundary_id boundary_id) const
          {
        std::string templ = Utilities::AsciiDataBase<dim>::data_directory + Utilities::AsciiDataBase<dim>::data_file_name;
        const int size = templ.length();
        const std::string boundary_name = this->get_geometry_model().translate_id_to_symbol_name(boundary_id);
        char *filename = (char *) (malloc ((size + 10) * sizeof(char)));
        snprintf (filename, size + 10, templ.c_str (), boundary_name.c_str(),filenumber);
        std::string str_filename (filename);
        free (filename);
        return str_filename;
          }


      template <int dim>
      void
      AsciiDataBoundary<dim>::update ()
      {
        if (time_dependent && (this->get_time() - first_data_file_model_time >= 0.0))
          {
            // whether we need to update our data files. This looks so complicated
            // because we need to catch increasing and decreasing file orders and all
            // possible first_data_file_model_times and first_data_file_numbers.
            const bool need_update =
                static_cast<int> ((this->get_time() - first_data_file_model_time) / data_file_time_step)
                > std::abs(current_file_number - first_data_file_number);

            if (need_update)
              for (typename std::map<types::boundary_id,
                  std_cxx11::shared_ptr<Utilities::AsciiDataLookup<dim-1> > >::iterator
                  boundary_id = lookups.begin();
                  boundary_id != lookups.end(); ++boundary_id)
                update_data(boundary_id->first);

            time_weight = (this->get_time() - first_data_file_model_time) / data_file_time_step
                - std::abs(current_file_number - first_data_file_number);

            Assert ((0 <= time_weight) && (time_weight <= 1),
                    ExcMessage (
                        "Error in set_current_time. Time_weight has to be in [0,1]"));
          }
      }

      template <int dim>
      void
      AsciiDataBoundary<dim>::update_data (const types::boundary_id boundary_id)
      {
        // The last file, which was tried to be loaded was
        // number current_file_number +/- 1, because current_file_number
        // is the file older than the current model time
        const int old_file_number =
            (decreasing_file_order) ?
                current_file_number - 1
                :
                current_file_number + 1;

        //Calculate new file_number
        current_file_number =
            (decreasing_file_order) ?
                first_data_file_number
                - static_cast<unsigned int> ((this->get_time() - first_data_file_model_time) / data_file_time_step)
                :
                first_data_file_number
                + static_cast<unsigned int> ((this->get_time() - first_data_file_model_time) / data_file_time_step);

        const int next_file_number =
            (decreasing_file_order) ?
                current_file_number - 1
                :
                current_file_number + 1;

        // If the time step was large enough to move forward more
        // then one data file we need to load both current files
        // to stay accurate in interpolation
        if (std::abs(current_file_number - old_file_number) >= 1)
          try
        {
              this->get_pcout() << std::endl << "   Loading Ascii data boundary file "
                  << create_filename (current_file_number,boundary_id) << "." << std::endl << std::endl;
              lookups.find(boundary_id)->second.swap(old_lookups.find(boundary_id)->second);
              lookups.find(boundary_id)->second->load_file(create_filename (current_file_number,boundary_id));
        }
        catch (...)
        // If loading current_time_step failed, end time dependent part with old_file_number.
        {
            try
            {
                end_time_dependence (old_file_number,boundary_id);
                return;
            }
            catch (...)
            {
                // If loading the old file fails (e.g. there was no old file), cancel the model run.
                // We might get here, if the model time step is so large that step t is before the
                // whole boundary condition while step t+1 is already behind all files in time.
                AssertThrow (false,
                    ExcMessage (
                        "Loading new and old data file did not succeed. "
                        "Maybe the time step was so large we jumped over all files "
                        "or the files were removed during the model run. "
                        "Another possible way here is to restart a model with "
                        "previously time-dependent boundary condition after the "
                        "last file was already read. Aspect has no way to find the "
                        "last readable file from the current model time. Please "
                        "prescribe the last data file manually in such a case. "
                        "Cancelling calculation."));
            }
        }

        // Now load the data file. This part is the main purpose of this function.
        try
        {
            this->get_pcout() << std::endl << "   Loading Ascii data boundary file "
                << create_filename (next_file_number,boundary_id) << "." << std::endl << std::endl;
            lookups.find(boundary_id)->second.swap(old_lookups.find(boundary_id)->second);
            lookups.find(boundary_id)->second->load_file(create_filename (next_file_number,boundary_id));
        }

        // If loading current_time_step + 1 failed, end time dependent part with current_time_step.
        // We do not need to check for success here, because current_file_number was guaranteed to be
        // at least tried to be loaded before, and if it fails, it should have done before (except from
        // hard drive errors, in which case the exception is the right thing to be thrown).

        catch (...)
        {
            end_time_dependence (current_file_number,boundary_id);
        }
      }

      template <int dim>
      void
      AsciiDataBoundary<dim>::end_time_dependence (const int file_number,
          const types::boundary_id boundary_id)
          {
        // no longer consider the problem time dependent from here on out
        // this cancels all attempts to read files at the next time steps
        time_dependent = false;
        // Give warning if first processor
        this->get_pcout() << std::endl
            << "   Loading new data file did not succeed." << std::endl
            << "   Assuming constant boundary conditions for rest of model run."
            << std::endl << std::endl;
          }

      template <int dim>
      double
      AsciiDataBoundary<dim>::
      get_data_component (const types::boundary_id             boundary_indicator,
          const Point<dim>                    &position,
          const unsigned int                   component) const
          {
        if (this->get_time() - first_data_file_model_time >= 0.0)
          {
            Point<dim> internal_position = position;

            if (dynamic_cast<const GeometryModel::SphericalShell<dim>*> (&this->get_geometry_model()) != 0)
              {
                const std_cxx11::array<double,dim> spherical_position =
                    ::aspect::Utilities::spherical_coordinates(position);

                for (unsigned int i = 0; i < dim; i++)
                  internal_position[i] = spherical_position[i];
              }

            const std_cxx11::array<unsigned int,dim-1> boundary_dimensions =
                get_boundary_dimensions(boundary_indicator);

            Point<dim-1> data_position;
            for (unsigned int i = 0; i < dim-1; i++)
              data_position[i] = internal_position[boundary_dimensions[i]];

            const double old_data = old_lookups.find(boundary_indicator)->second->get_data(data_position,component);

            if (!time_dependent)
              return old_data;

            const double data = lookups.find(boundary_indicator)->second->get_data(data_position,component);

            return time_weight * data + (1 - time_weight) * old_data;
          }
        else
          return 0.0;
          }


      template <int dim>
      void
      AsciiDataBoundary<dim>::declare_parameters (ParameterHandler  &prm,
          const std::string &default_directory,
          const std::string &default_filename)
          {
        Utilities::AsciiDataBase<dim>::declare_parameters(prm,
            default_directory,
            default_filename);

        prm.enter_subsection ("Ascii data model");
        {
          prm.declare_entry ("Data file time step", "1e6",
              Patterns::Double (0),
              "Time step between following velocity files. "
              "Depending on the setting of the global 'Use years in output instead of seconds' flag "
              "in the input file, this number is either interpreted as seconds or as years. "
              "The default is one million, i.e., either one million seconds or one million years.");
          prm.declare_entry ("First data file model time", "0",
              Patterns::Double (0),
              "Time from which on the velocity file with number 'First velocity "
              "file number' is used as boundary condition. Previous to this "
              "time, a no-slip boundary condition is assumed. Depending on the setting "
              "of the global 'Use years in output instead of seconds' flag "
              "in the input file, this number is either interpreted as seconds or as years.");
          prm.declare_entry ("First data file number", "0",
              Patterns::Integer (),
              "Number of the first velocity file to be loaded when the model time "
              "is larger than 'First velocity file model time'.");
          prm.declare_entry ("Decreasing file order", "false",
              Patterns::Bool (),
              "In some cases the boundary files are not numbered in increasing "
              "but in decreasing order (e.g. 'Ma BP'). If this flag is set to "
              "'True' the plugin will first load the file with the number "
              "'First velocity file number' and decrease the file number during "
              "the model run.");
        }
        prm.leave_subsection();
          }


      template <int dim>
      void
      AsciiDataBoundary<dim>::parse_parameters (ParameterHandler &prm)
      {
        Utilities::AsciiDataBase<dim>::parse_parameters(prm);

        prm.enter_subsection("Ascii data model");
        {
          data_file_time_step             = prm.get_double ("Data file time step");
          first_data_file_model_time      = prm.get_double ("First data file model time");
          first_data_file_number          = prm.get_double ("First data file number");
          decreasing_file_order           = prm.get_bool   ("Decreasing file order");

          if (this->convert_output_to_years() == true)
            {
              data_file_time_step        *= year_in_seconds;
              first_data_file_model_time *= year_in_seconds;
            }
        }
        prm.leave_subsection();
      }
    }


    template <int dim>
    AsciiData<dim>::AsciiData ()
    {}


    template <int dim>
    void
    AsciiData<dim>::initialize ()
    {
      const std::set<types::boundary_id> boundary_ids = this->get_fixed_temperature_boundary_indicators();
      internal::AsciiDataBoundary<dim>::initialize(boundary_ids,
                                                    1);
    }


    template <int dim>
    void
    AsciiData<dim>::update ()
    {
      Interface<dim>::update ();
      internal::AsciiDataBoundary<dim>::update();
    }


    template <int dim>
    double
    AsciiData<dim>::temperature (const GeometryModel::Interface<dim> &geometry_model,
                                 const unsigned int                   boundary_indicator,
                                 const Point<dim>                    &position) const
    {
      const types::boundary_id boundary_id(boundary_indicator);
      return internal::AsciiDataBoundary<dim>::get_data_component(boundary_id,
                                                                   position,
                                                                   0);
    }


    template <int dim>
    double
    AsciiData<dim>::minimal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids) const
    {
      return 0;
    }


    template <int dim>
    double
    AsciiData<dim>::maximal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids) const
    {
      return 0;
    }


    template <int dim>
     void
     AsciiData<dim>::declare_parameters (ParameterHandler &prm)
     {
       prm.enter_subsection("Boundary temperature model");
       {
         internal::AsciiDataBoundary<dim>::declare_parameters(prm,
                                                               "$ASPECT_SOURCE_DIR/data/boundary-temperature/ascii-data/test/",
                                                               "box_2d_%s.%d.txt");
       }
       prm.leave_subsection();
     }


     template <int dim>
     void
     AsciiData<dim>::parse_parameters (ParameterHandler &prm)
     {
       prm.enter_subsection("Boundary temperature model");
       {
         internal::AsciiDataBoundary<dim>::parse_parameters(prm);
       }
       prm.leave_subsection();
     }
  }
}

// explicit instantiations
namespace aspect
{
  namespace BoundaryTemperature
  {
    namespace internal
    {
    template class AsciiDataBoundary<2>;
    template class AsciiDataBoundary<3>;
    }

    ASPECT_REGISTER_BOUNDARY_TEMPERATURE_MODEL(AsciiData,
                                               "ascii data",
                                               "Implementation of a model in which the boundary "
                                               "data is derived from files containing data "
                                               "in ascii format. Note the required format of the "
                                               "input data: The first lines may contain any number of comments"
                                               "if they begin with '#', but one of these lines needs to"
                                               "contain the number of grid points in each dimension as"
                                               "for example '# POINTS: 3 3'."
                                               "The order of the data columns "
                                               "has to be 'x', 'Temperature [K]' in a 2d model and "
                                               " 'x', 'y', 'Temperature [K]' in a 3d model, which means that "
                                               "there has to be a single column "
                                               "containing the temperature. "
                                               "Note that the data in the input "
                                               "files need to be sorted in a specific order: "
                                               "the first coordinate needs to ascend first, "
                                               "followed by the second in order to "
                                               "assign the correct data to the prescribed coordinates."
                                               "If you use a spherical model, "
                                               "then the data will still be handled as cartesian,"
                                               "however the assumed grid changes. 'x' will be replaced by "
                                               "the radial distance of the point to the bottom of the model, "
                                               "'y' by the azimuth angle and 'z' by the polar angle measured "
                                               "positive from the north pole. The grid will be assumed to be "
                                               "a latitude-longitude grid. Note that the order "
                                               "of spherical coordinates is 'r', 'phi', 'theta' "
                                               "and not 'r', 'theta', 'phi', since this allows "
                                               "for dimension independent expressions. ")
  }
}
